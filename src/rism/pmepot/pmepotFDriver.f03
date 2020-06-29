program pmepotFDriver
  implicit none
  integer, parameter :: numAtoms = 3
  real (kind = 8) :: positions(3, numAtoms)
  real (kind = 8) :: charges(numAtoms) = (/ 1, -1, 1 /)
  integer :: gridDim(3) = (/ 10, 10, 10 /)
  real (kind = 8) :: boxLen(3) = (/ 5, 5, 5 /)
  real (kind = 8) :: potential(10*10*10)
  
  integer :: atom
  real (kind = 8) :: positionOffset = 0

  real (kind = 8) :: numbas(3,3) = 0
  numbas(:,1) = (/ 1, 0, 0 /)
  numbas(:,2) = (/ 0, 1, 0 /)
  numbas(:,3) = (/ 0, 0, 1 /)

  do atom = 1, numAtoms
     positions(:, atom) = positionOffset
     positionOffset = positionOffset + 0.5
  end do
  
  print *, "Hello world!"

  call pmepot(numAtoms, positions, charges, &
       gridDim, boxLen, &
       25d0, numbas, numbas, &
       1d0, 3, 100d0, &
       potential)

  ! open(unit = 1, file = 'potential.dat', )

  print *, potential
  
end program pmepotFDriver
  
