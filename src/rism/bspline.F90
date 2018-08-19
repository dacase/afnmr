!<compile=optimized>

! The 3D-RISM-KH software found here is copyright (c) 2010-2012 by 
! Andriy Kovalenko, Tyler Luchko and David A. Case.
! 
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
! 
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
! 
! You should have received a copy of the GNU General Public License in the
! ../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
! 
! Users of the 3D-RISM capability found here are requested to acknowledge
! use of the software in reports and publications.  Such acknowledgement
! should include the following citations:
! 
! 1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
! ibid. 112:10391-10417 (2000).   
! 
! 2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
! F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
! 
! 3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
! and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"

!> Assorted b-spline functions.
!!TODO: Document all methods.
module bspline

  interface cardinal_bspline
     module procedure cardinal_bspline_weights
     module procedure cardinal_bspline_weights_derivs
  end interface
  
contains

  subroutine cardinal_bspline_recursion(offset, order, weights)
    implicit none  
    _REAL_, intent(in) :: offset
    integer, intent(in) :: order
    _REAL_, intent(out) :: weights(order)

    _REAL_ :: invOrderM1
    integer :: j

    invOrderM1 = 1d0 / (order - 1)

    weights(order) = offset * weights(order - 1) * invOrderM1
    do j = 1, order - 2
       weights(order - j) = ((offset + j) * weights(order - j - 1) &
            + (order - j - offset) * weights(order - j)) * invOrderM1
    end do
    weights(1) = (1 - offset) * weights(1) * invOrderM1
  end subroutine cardinal_bspline_recursion


  subroutine cardinal_bspline_weights(offset, order, weights)
    implicit none
    _REAL_, intent(in) :: offset
    integer, intent(in) :: order
    _REAL_, intent(out) :: weights(order)
    integer :: order_n
    _REAL_, parameter :: div = 1. / 3.

    weights(2) = offset
    weights(1) = 1 - offset
    if (order > 2) then
       ! One pass to order 3:
       weights(3) = 0.5 * offset * weights(2)
       weights(2) = 0.5 * ((offset + 1.) * weights(1) + (2. - offset) * weights(2))
       weights(1) = 0.5 * (1. - offset) * weights(1)
       if (order > 3) then
          ! One pass to order 4:
          weights(4) = div * offset * weights(3)
          weights(3) = div * ((offset + 1.) * weights(2) + (3. - offset) * weights(3))
          weights(2) = div * ((offset + 2.) * weights(1) + (2. - offset) * weights(2))
          weights(1) = div * (1. - offset) * weights(1)
          ! And the rest.
          do order_n = 5, order
             call cardinal_bspline_recursion(offset, order_n, weights)
          end do
       end if
    end if
  end subroutine cardinal_bspline_weights

  
  !> Compute the derivative of the provided b-spline weights. Multiple
  !! calls result in higher order derivatives.
  !! The weights should be calculated up to order - 1, NOT order.
  subroutine differentiate_cardinal_bspline(order, weights, deriv_weights)
    implicit none
    integer, intent(in) :: order
    _REAL_, intent(in) :: weights(:)
    _REAL_, intent(out) :: deriv_weights(:)

    integer :: i

    deriv_weights(1) = -weights(1)

    do i = 2, order - 1
       deriv_weights(i) = weights(i - 1) - weights(i)
    end do

    deriv_weights(order) = weights(order - 1)
  end subroutine differentiate_cardinal_bspline

  
  !> Calculate the b-spline weights and their derivatives.
  subroutine cardinal_bspline_weights_derivs( &
       offset, order, weights, weight_derivs)
    implicit none
    _REAL_, intent(in) :: offset
    integer, intent(in) :: order
    _REAL_, intent(out) :: weights(order)
    _REAL_, intent(out) :: weight_derivs(order)

    ! Calculate weights up to order - 1.
    call cardinal_bspline(offset, order - 1, weights)

    ! Calculate the derivatives of the b-spline weights with respect
    ! to reciprocal space position.
    call differentiate_cardinal_bspline(order, weights, weight_derivs)

    ! One more recursion to obtain weights at desired order.
    call cardinal_bspline_recursion(offset, order, weights)
  end  subroutine cardinal_bspline_weights_derivs

  
  !> Use chain rule to convert bspline derivative from reciprocal
  !! space to real space.
  !! TODO: More efficient to do this outside any loops, so this
  !! subroutine should be replaced by manual multiply by user.
  subroutine bspline_deriv_k_to_r(boxLen, gridDim, order, weight_derivs)
    _REAL_, intent(in) :: boxLen
    integer, intent(in) :: gridDim
    integer, intent(in) :: order
    _REAL_, intent(inout) :: weight_derivs(order)
    integer :: i
    ! Chain rule: (dw / dp) * (dp / dx) where w is the b-spline
    ! weight, p is the reciprocal space position used to evaluate the
    ! b-spline, and x is the real space position.
    do i = 1, order
       weight_derivs(i) = weight_derivs(i) * gridDim / boxLen
    end do
  end subroutine bspline_deriv_k_to_r


  !> Calculate b-spline weights on a periodic grid, their grid point
  !! locations, and, optionally, their derivatives with respect to
  !! reciprocal space position.
  subroutine cardinal_bspline_periodic_grid( &
       x, boxLen, gridDim, order, nearestGridPoints, weights, weight_derivs)
    implicit none
    _REAL_, intent(in) :: x
    _REAL_, intent(in) :: boxLen
    integer, intent(in) :: gridDim
    integer, intent(in) :: order
    integer, intent(out) :: nearestGridPoints(:)
    _REAL_, intent(out) :: weights(:)
    _REAL_, intent(out), optional :: weight_derivs(:)

    _REAL_ :: distance

    call nearby_grid_points_periodic(x, boxLen, gridDim, order, &
         nearestGridPoints, distance)

    if (present(weight_derivs)) then
       call cardinal_bspline(distance, order, weights, weight_derivs)
    else
       call cardinal_bspline(distance, order, weights)
    end if
  end subroutine cardinal_bspline_periodic_grid

  
  subroutine nearby_grid_points_periodic( &
       x, boxLen, gridDim, order, nearestGridPoints, distance)
    implicit none
    _REAL_, intent(in) :: x
    _REAL_, intent(in) :: boxLen
    integer, intent(in) :: gridDim
    integer, intent(in) :: order
    integer, intent(out) :: nearestGridPoints(:)
    _REAL_, intent(out), optional :: distance !< Distance from nearest
                                              !! gridpoint to the left.
    _REAL_ :: xOnGrid
    integer :: nearbyGridPoint
    integer :: i

    ! Find the distance (in terms of fraction grid points) to the
    ! closest grid point left of x.
    xOnGrid = x / (boxLen / gridDim) - mod(order, 2) * 0.5
    nearbyGridPoint = floor(xOnGrid)
    distance = xOnGrid - nearbyGridPoint
    ! Identify the grid points where the b-spline will be evaluated,
    ! starting from the leftmost and considering periodicity such that
    ! all grid points are kept within the reference box.
    nearbyGridPoint = mod(mod(nearbyGridPoint - (order / 2 - 1), gridDim) + gridDim, gridDim)

    do i = 1, order
       nearestGridPoints(i) = mod(nearbyGridPoint, gridDim)
       nearbyGridPoint = nearbyGridPoint + 1
    end do
  end subroutine nearby_grid_points_periodic

  
  !> Calculate the Fourier coefficients of the cardinal b-spline.
  subroutine cardinal_bspline_Fourier_coefficients( &
       maxIndex, order, weights, fourierCoefficients, half)
    use constants, only : pi
    implicit none
    integer, intent(in) :: maxIndex !< Largest Fourier coefficient index to compute.
    integer, intent(in) :: order !< 
    _REAL_, intent(in) :: weights(order) !< 
    _REAL_, intent(out) :: fourierCoefficients(:) !< 
    logical, intent(in) :: half !< Calculate Fourier coefficients for only
                                !! half the range.

    integer :: m, n, i
    _REAL_ :: k
    integer :: o
    integer :: Nh
    integer :: Np
    _REAL_ :: tpion
    
    fourierCoefficients(:) = 0

    ! o = order / 2 - 1 for even order, order / 2 for odd order.
    o = (order - 1) / 2
    Nh = maxIndex - o
    Np = merge(maxIndex / 2 + 1, maxIndex, half)

    tpion = 2 * pi / maxIndex
    do m = 1, Np
       k = (m - 1) * tpion
       ! Iterates from weights((order - 1) / 2) to weights(order).
       do n = 0, order / 2
          fourierCoefficients(m) = fourierCoefficients(m) &
               + weights(n + o + 1) * cos(k * n)
       end do
       ! Iterates from weights(1) to weights((order - 1) / 2).
       do i = 0, o - 1
          fourierCoefficients(m) = fourierCoefficients(m) &
               + weights(i + 1) * cos(k * (Nh + i))
       end do
    end do
  end subroutine cardinal_bspline_Fourier_coefficients

end module bspline
