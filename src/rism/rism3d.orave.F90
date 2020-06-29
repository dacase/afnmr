#include "../include/dprec.fh"

!The 3D-RISM-KH software found here is copyright (c) 2012 by
!Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

!> This program takes 3D solvent distribution functions from
!! 3D-RISM calculations and calculates orientational averages.  
module orave_m
  use safemem
  use getopts_c
  use rism_report_c
  use rism3d_grid_c
  use quaternion
  implicit none

  type dimension
     !len : length in [A]
     !d   : delta in [A]
     !min : smallest value in range [A]
     !max : largest value in range [A]
     !twist : deform the coordinates in perpendicular plane by
     !        applying a continuous rotation along this dimension.
     !        The default is for unwinding B-form DNA along the
     !        cylindrical axis.
     !        phi = phi - twist*i*dr
     _REAL_ :: len=0d0, d=0d0, max=0d0, min=0d0, twist=0.18587d0
     !n   : number of point (length/delta)
     !i   : the current index in this dimension
     integer :: n = 0, i
     !name : name of dimension
     character(len=16) :: name
     !applytwist :: apply twist to this dimension
     logical :: applytwist=.false.
  end type dimension

  integer, parameter :: CLEN=1024
  !possible coordinate types
  integer, parameter :: SPHERICAL=0, CYLINDRICAL=1, CARTESIAN=2
  !possible coordinate names
  integer, parameter :: XDIM=32,YDIM=16,ZDIM=8,RDIM=4,THETADIM=2,PHIDIM=1

  !system : type of coodinate system
  !dim    : sum of dimension.  E.g. X+Y
  integer :: system, dim

  !each possible dimension
  type(dimension),target,save :: x,y,z,radius,theta,phi
  !dummy dimension for output. a dimension with only one point.  Used for second
  !dimension for 1D output
  type(dimension),target,save :: dummy_dim

  !grid  :: grid object
  type(rism3d_grid),save :: grid

  !guvfile : filename for the input guv files
  character(len=CLEN),pointer :: guvfile(:)=>NULL()
  
  !method  : integration method for averaging selected by user
  !outfile : name of the output file, if it exists
  !pdbfile : pdbfile that contains all sample points
  character(len=CLEN) :: method, outfile='', pdbfile='test.pdb'

  !guv        :: (product(nr)) solvent distribution function
  _REAL_,pointer :: guv(:)=>NULL()

  !dxOrigin :: origin from the DX file.  Tells us how to translate the solute and
  !            how to interpret volume specifications
  !origin :: origin for integration
  !vec    :: orientation vector for cylindrical averaging
  !orient_quat :: quaterion that can be used to rotate data points to
  !               new coordinate system defined by vec
  _REAL_ :: dxOrigin(3), origin(3),vec(3), orient_quat(4)

  !abscissa :: list of abscissa for the particular method. Spherical
  !            averaging only.
  !grid_out :: averaged data for 1D or 2D averaging for each site (last index).
  !            If there is only one output dimension, the second dimension
  !            of grid is set to 1
  _REAL_, pointer :: abscissa(:,:)=>NULL(), grid_out(:,:,:)=>NULL()

  !dim_out1 :: output dimension1 
  !dim_out2 :: output dimension2 
  type(dimension), pointer :: dim_out1, dim_out2

  !crd1     :: first dimension in this coordinate system
  !crd2     :: second dimension in this coordinate system
  !crd3     :: third dimension in this coordinate system
  type(dimension), pointer :: crd1, crd2, crd3

  logical ::  writePDB=.false., verbose=.false.

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!get command line options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getOptions()
    use constants, only : PI
    implicit none
    integer :: err,i, ntoken
    character(len=CLEN),pointer :: extra(:)=>NULL()
    character(len=10*CLEN), pointer :: tempC(:)=>NULL()
    character(len=10*CLEN) :: tempOpt, strpair(2)
    _REAL_ :: twist
    integer :: iarg, narg
    type(dimension),pointer :: tempdim
    logical :: help, true 
    call getopts_add("system",'',min=1,max=1)
    call getopts_add("zvec",'0,0,1',"v",1)
    call getopts_add("origin",'default',"o",1)
    call getopts_add("x",.false.,"x")
    call getopts_add("y",.false.,"y")
    call getopts_add("z",.false.,"z")
    call getopts_add("r",.false.,"radius")
    call getopts_add("theta",.false.,"t")
    call getopts_add("phi",.false.,"p")
    call getopts_add("lenx",0d0,max=1)
    call getopts_add("leny",0d0,max=1)
    call getopts_add("lenz",0d0,max=1)
    call getopts_add("lenr",0d0,max=1)
    call getopts_add("dr",0d0,max=1)
    call getopts_add("dx",0d0,max=1)
    call getopts_add("dy",0d0,max=1)
    call getopts_add("dz",0d0,max=1)
    call getopts_add("dtheta",0d0,max=1)
    call getopts_add("dphi",0d0,max=1)
    call getopts_add("nx",0,max=1)
    call getopts_add("ny",0,max=1)
    call getopts_add("nz",0,max=1)
    call getopts_add("nr",0,max=1)
    call getopts_add("ntheta",0,max=1)
    call getopts_add("nphi",0,max=1)
    call getopts_add("twist",'')
    call getopts_add("guv",'',"g",min=1)
    call getopts_add("method",'uniform',max=1)
    call getopts_add("pdbout",'',max=1)
    call getopts_add("verbose",.false.)
    call getopts_add("help",.false.,"h")
    ntoken = getopts_process()
!    call getopts_summary(0)

    call getopts_get("help",help)
    if(help .or. ntoken == 0)then
       call usage()
       stop
    elseif(ntoken <0)then
       call usage()
       call rism_report_error("command line read failed.")
    end if
    
    !get input and output files
    guvfile=> getopts_getAll("guv",guvfile,len(guvfile))
    extra=> getopts_unread()
    if(size(extra) == 1) outfile = extra(1)
    if( getopts_narg("guv") == 0) then
       call usage()
       call rism_report_error("'--guv' must be specified")
    else if(size(extra) >1)then
       call usage()
       call flush(0)
       !concatenate all extra options into temp
       tempOpt=""
       do i=2, size(extra)
          tempOpt(len_trim(tempOpt)+1:) =  " "//trim(extra(i))
       end do
       call rism_report_error("unknown options:"//trim(tempOpt))
    end if

    !check if the user wants the debug output
    if(getopts_narg("pdbout") == 1)then
       call getopts_get("pdbout",1,pdbfile)
       writePDB=.true.
    end if
    
    !set up the averaging type.
    !We have three options for cooridinate systems:
    !  spherical, cylindrical and cartesian
    call getopts_get("system",1,tempopt)
    if(tempopt .eq. 'spherical')then
       system = SPHERICAL
       call rism_report_error("spherical coodinates not supported at this time.")
    elseif(tempopt .eq. 'cylindrical')then
       system = CYLINDRICAL
    elseif(tempopt .eq. 'cartesian')then
       system = CARTESIAN
    else
       call usage()
       call rism_report_error("unknown coordinate system: "//trim(tempopt))
    end if

    !get system origin and orientation

    !for comma separated arguments, read in arguments as strings, then parse
    !arguments into the arrays
    call getopts_get("origin",1,tempopt)
    if(trim(tempopt) .eq. 'default')then
       origin = huge(1d0)
    else
       read(tempopt,*,iostat=err) origin
       if(err/=0) call rism_report_error("Invalid cooridinates for --origin: "&
            //trim(tempopt))
    end if
    call getopts_get("zvec",1,tempopt)
    read(tempopt,*,iostat=err) vec
    if(err/=0) call rism_report_error("Invalid cooridinates for --zvec: "&
         //trim(tempopt))
    !make it a unit vector
    vec = vec/sqrt(sum(vec**2))

    !Find the dimensions that the user wants to keep.  We can check
    !for illegal values later
    dim=0
    call getopts_get("x",true)
    if(true)dim=dim+XDIM
    call getopts_get("y",true)
    if(true)dim=dim+YDIM
    call getopts_get("z",true)
    if(true)dim=dim+ZDIM
    call getopts_get("r",true)
    if(true)dim=dim+RDIM
    call getopts_get("theta",true)
    if(true)dim=dim+THETADIM
    call getopts_get("phi",true)
    if(true)dim=dim+PHIDIM
    !quick check
    if(dim==0) &
         call rism_report_error("atleast one output dimension must be specified")
    
    !collect box dimensions. Extra values can be ignored
    !X
    x%name="x"
    if( getopts_narg("lenx") == 1)&
         call getopts_get("lenx",1,x%len)
    if( getopts_narg("dx") == 1)&
         call getopts_get("dx",1,x%d)
    if( getopts_narg("nx") == 1)&
         call getopts_get("nx",1,x%n)
    !Y
    y%name="y"
    if( getopts_narg("leny") == 1)&
         call getopts_get("leny",1,y%len)
    if( getopts_narg("dy") == 1)&
         call getopts_get("dy",1,y%d)
    if( getopts_narg("ny") == 1)&
         call getopts_get("ny",1,y%n)
    !Z
    z%name="z"
    if( getopts_narg("lenz") == 1)&
         call getopts_get("lenz",1,z%len)
    if( getopts_narg("dz") == 1)&
         call getopts_get("dz",1,z%d)
    if( getopts_narg("nz") == 1)&
         call getopts_get("nz",1,z%n)
    !R
    radius%name="r"
    if( getopts_narg("lenr") == 1)&
         call getopts_get("lenr",1,radius%len)
    if( getopts_narg("dr") == 1)&
         call getopts_get("dr",1,radius%d)
    if( getopts_narg("nr") == 1)&
         call getopts_get("nr",1,radius%n)
    !THETA
    theta%name="theta"
    theta%len=2d0*PI
    if( getopts_narg("dtheta") == 1)&
         call getopts_get("dtheta",1,theta%d)
    if( getopts_narg("ntheta") == 1)&
         call getopts_get("ntheta",1,theta%n)
    !PHI
    phi%name="phi"
    phi%len=PI
    if( getopts_narg("dphi") == 1)&
         call getopts_get("dphi",1,phi%d)
    if( getopts_narg("nphi") == 1)&
         call getopts_get("nphi",1,phi%n)

    !get twist argument for unwinding along a dimension
    do iarg=1,getopts_narg("twist")
       call getopts_get("twist",iarg,tempopt)
       strpair='default'
       !list-directed input is whitespace/comma/slash delimited
       read(tempopt,*,iostat=err) strpair
       !check which dimension
       if(trim(strpair(1)) .eq. "x")then
          tempdim=>x
       elseif(trim(strpair(1)) .eq. "y")then
          tempdim=>y
       elseif(trim(strpair(1)) .eq. "z")then
          tempdim=>z
       elseif(trim(strpair(1)) .eq. "r" .or. trim(strpair(1)) .eq. "radius")then
          tempdim=>radius
       elseif(trim(strpair(1)) .eq. "t" .or. trim(strpair(1)) .eq. "theta")then
          tempdim=>theta
       elseif(trim(strpair(1)) .eq. "p" .or. trim(strpair(1)) .eq. "phi")then
          tempdim=>phi
       else
          call rism_report_error("unknown dimension to --twist: "//trim(strpair(1)))
       end if
       !check if the optional argument was passed.  Currently not checking bad input
       tempdim%applytwist=.true.
       read(strpair(2),*,iostat=err) twist
       if(err==0) tempdim%twist=twist
    end do
    !integration method
    call getopts_get("method",1,method)


    if(safemem_dealloc(extra)/=0) &
         call rism_report_error("GETOPTONS: failed to dealloc 'extra'")
    if(safemem_dealloc(tempC)/=0) &
         call rism_report_error("GETOPTONS: failed to dealloc 'tempC'")


    call getopts_get("verbose",verbose)
    if(verbose) call print_parameters()
    call getopts_cleanup()
  end subroutine getOptions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!read in input 3D distribution files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDistributionFiles
    use rism3d_opendx
    implicit none
    integer :: i
    _REAL_, pointer :: tmp1(:)
    _REAL_ :: delta(3)
    integer, parameter :: zeroR3(3)=0d0
    integer :: npos(3), nkpos(3), tmppos(3)

    !
    !User supplied distribution functions
    !

    !initialize objects
    call rism3d_grid_new(grid)

    !get grid size
    if(ubound(guvfile,1) >0)then
       call readDXHeader(guvfile(1),dxOrigin,delta,npos)
    end if

    !set grid size
    nkpos = npos
    nkpos(1) = nkpos(1)+2
    call rism3d_grid_resize(grid,delta,npos,nkpos,npos,nkpos,zeroR3,zeroR3)
    guv => safemem_realloc(guv,grid%totalLocalPointsR,.false.)
    !Guv
    if(ubound(guvfile,1) > 0)then
       do i = 1, ubound(guvfile,1)
          !check file size
          call readDXHeader(guvfile(i),dxOrigin,delta,tmppos)
          if(sum(abs(tmppos- npos))/=0) call rism_report_error(trim(guvfile(i))//" is the wrong size")
       end do
    end if
    if(origin(1) == huge(1d0)) origin = dxOrigin + grid%boxLength/2d0
  end subroutine readDistributionFiles


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!allocates memory and sets values for abscissa and weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setAbscissa()
    use constants, only : pi
    implicit none
    integer :: ix,iy,iz,itheta
    _REAL_ :: cross(3), angle, point(3)
    _REAL_ :: xaxis(3)=(/1,0,0/), yaxis(3)=(/0,1,0/), zaxis(3)=(/0,0,1/)
    !check that the origin is in the datagrid
    if(.not.inGrid(origin))&
         call rism_report_error("Origin must be located in solvent grid.")
    !set up rotation quaterion
    cross(1) = zaxis(2)*vec(3) - zaxis(3)*vec(2)
    cross(2) = zaxis(3)*vec(1) - zaxis(1)*vec(3)
    cross(3) = zaxis(1)*vec(2) - zaxis(2)*vec(1)
    angle = acos(dot_product(zaxis, vec))
    if(angle == 0d0) cross = (/0,0,1/)
    orient_quat = quaternionFromEulerAxis(angle,cross)

    !initialize the dummy dimension in case it is needed
    dummy_dim%n=1
    dummy_dim%i=0
    dummy_dim%name="dummy"

    if(system==CARTESIAN)then
       !1) complete gridding information for each dimension
       call complete_dimension(x)
       call complete_dimension(y)
       call complete_dimension(z)
       crd1=>x
       crd2=>y
       crd3=>z
       !2) check corners are in data grid
       do iz=-1,1,2
          do iy=-1,1,2
             do ix=-1,1,2
                point = ix*xaxis*x%len/2d0 + iy*yaxis*y%len/2d0 + iz*zaxis*z%len/2d0
                call rotate_quat(point,orient_quat)
                point=point+origin
                if(.not.inGrid(point))&
                     call rism_report_error("defined cartesian box does not "&
                     //"fit in data with given origin and orientation")
             end do
          end do
       end do

       !3) check that we have a valid combination of cartesian
       !dimensions and allocate grid
       if(dim==XDIM)then
          grid_out => safemem_realloc(grid_out,x%n,1,ubound(guvfile,1))
          dim_out1=>x
          dim_out2=>dummy_dim
       elseif(dim==YDIM)then
          grid_out => safemem_realloc(grid_out,y%n,1,ubound(guvfile,1))
          dim_out1=>y
          dim_out2=>dummy_dim
       elseif(dim==ZDIM)then
          grid_out => safemem_realloc(grid_out,z%n,1,ubound(guvfile,1))
          dim_out1=>z
          dim_out2=>dummy_dim
       elseif(dim==XDIM+YDIM)then
          grid_out => safemem_realloc(grid_out,x%n,y%n,ubound(guvfile,1))
          dim_out1=>x
          dim_out2=>y
       elseif(dim==YDIM+ZDIM)then
          grid_out => safemem_realloc(grid_out,y%n,z%n,ubound(guvfile,1))
          dim_out1=>y
          dim_out2=>z
       elseif(dim==ZDIM+XDIM)then
          grid_out => safemem_realloc(grid_out,x%n,z%n,ubound(guvfile,1))
          dim_out1=>x
          dim_out2=>z
       else
          call rism_report_error("illegal combination of output axes")
       end if
    elseif(system==CYLINDRICAL)then
       !1) complete gridding information for each dimension
       call complete_dimension(radius)
       call complete_dimension(theta)
       call complete_dimension(z)
       crd1=>radius
       crd2=>theta
       crd3=>z
       !2) check that ends are in data grid
       !cycle through each point in the end caps of the cylinder and check if they are in the grid
       do iz=-1,1,2
          do itheta = 0, theta%n-1
             point = radius%len*cos(itheta*theta%d)*xaxis &
                  + radius%len*sin(itheta*theta%d)*yaxis &
                  + iz*z%len/2d0*zaxis
             call rotate_quat(point,orient_quat)
             point=point+origin
             if(.not.inGrid(point))&
                  call rism_report_error("defined cylinder does not "&
                  //"fit in data with given origin and orientation")
          end do
       end do
       !3) check that we have a valid combination of cylindrical
       !dimensions and allocate grid
       if(dim==RDIM)then
          grid_out => safemem_realloc(grid_out,radius%n,1,ubound(guvfile,1))
          dim_out1=>radius
          dim_out2=>dummy_dim
       elseif(dim==THETADIM)then
          grid_out => safemem_realloc(grid_out,theta%n,1,ubound(guvfile,1))
          dim_out1=>theta
          dim_out2=>dummy_dim
       elseif(dim==ZDIM)then
          grid_out => safemem_realloc(grid_out,z%n,1,ubound(guvfile,1))
          dim_out1=>z
          dim_out2=>dummy_dim
       elseif(dim==RDIM+THETADIM)then
          grid_out => safemem_realloc(grid_out,radius%n,theta%n,ubound(guvfile,1))
          dim_out1=>radius
          dim_out2=>theta
       elseif(dim==RDIM+ZDIM)then
          grid_out => safemem_realloc(grid_out,radius%n,z%n,ubound(guvfile,1))
          dim_out1=>radius
          dim_out2=>z
       elseif(dim==THETADIM+ZDIM)then
          grid_out => safemem_realloc(grid_out,theta%n,z%n,ubound(guvfile,1))
          dim_out1=>theta
          dim_out2=>z
       else
          call rism_report_error("illegal combination of output axes")
       end if
    else
       call rism_report_error("Spherical averaging not yet supported.")
    end if
    if(verbose) call print_gridspecs()
  end subroutine setAbscissa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!fills in/asserts consistency of length, spacing and point number
!!!for dimension d.  Also handles twist.  If applytwist=.false. then twist=0
!!!IN:
!!!   d - the dimension to complete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine complete_dimension(d)
    implicit none
    type(dimension), intent(inout) :: d
    _REAL_ :: diff
    !d%d
    if(d%d == 0d0)then
       if(d%len == 0d0 .or. d%n==0)&
            call rism_report_error("dimension '"//trim(d%name)&
            //"' requires two of --len"//trim(d%name)//", --d"//trim(d%name)&
            //" or n"//trim(d%name)//".")
       d%d = d%len/(d%n-1)
    end if
    !d%len
    if(d%len == 0d0)then
       if(d%d == 0d0 .or. d%n==0)&
            call rism_report_error("dimension '"//trim(d%name)&
            //"' requires two of --len"//trim(d%name)//", --d"//trim(d%name)&
            //" or n"//trim(d%name)//".")
       d%len = d%d*(d%n-1)
    end if
    !d%n
    if(d%n == 0)then
       if(d%d == 0d0 .or. d%len==0d0)&
            call rism_report_error("dimension '"//trim(d%name)&
            //"' requires two of --len"//trim(d%name)//", --d"//trim(d%name)&
            //" or n"//trim(d%name)//".")
       d%n = d%len/d%d + 1
    end if
    !special case for a single data point
    if(d%n == 1)then
       d%d = 0d0
       d%len=0d0
    end if

    !check consistency
    if(abs(d%len - d%d*(d%n-1)) > EPSILON(1d0)*10d6)&
         call rism_report_error('(a,e16.8,a,e16.8,a,e16.8,a,e16.8)',"Inconsisent values for '"//trim(d%name)&
         //"' length != delta*(N-1): ",d%len," != ",  d%d*(d%n-1), " = ", d%d, " * ", dble(d%n-1))
    !set bounds
    if(d%name .eq. "r" .or. d%name .eq. "theta" .or. d%name .eq. "phi")then
       d%min=0d0
       d%max=d%len
    else
       d%min=-d%len/2d0
       d%max=d%len/2d0
    end if
    !check if twist is set
    if(.not.d%applytwist)&
         d%twist=0d0
  end subroutine complete_dimension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates orientational averages
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine average()
    use rism3d_opendx
    implicit none
    integer :: iguv
    _REAL_ :: tmporigin(3), grdspc(3)
    tmporigin=-1d0
    grdspc=-1d0
    grid_out=0d0
    do iguv = 1, ubound(guvfile,1)
       call readDX(guvfile(iguv),guv,grid%localDimsR,tmporigin,grdspc)
       if(system==CARTESIAN)then
          call averageCartesian(iguv)
       elseif(system==CYLINDRICAL)then
          call averageCylindrical(iguv)
       else
       end if
    end do
    grid_out=grid_out*(dim_out1%n*dim_out2%n)/(crd1%n*crd2%n*crd3%n)
  end subroutine average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates cartesian coordinate orientational averages
!!!IN:
!!!   iguv : the sequential number of the guv file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine averageCartesian(iguv)
    use rism_util, only : freeUnit
    use constants, only : pi
    implicit none
    integer, intent(in) :: iguv
    !point : point before rotation
    !rotatepoint : point after rotation
    !qtwist : an array of quaterions, one for each dimension, that are
    !         used to hold the twist rotation for this coordinate
    _REAL_ :: point(3), rotatepoint(3), qtwist(4,3), qfinal(4)
    integer :: i,j,k
    !the gfortran won't let me use %i as a counter so we use i,j and k
    !and then set %i
    do k=0, crd3%n-1
       crd3%i=k
       point(3) = crd3%i*crd3%d + crd3%min
       qtwist(:,3) = quaternionFromEulerAxis(point(3)*crd3%twist,(/0d0,0d0,1d0/))
       do j=0, crd2%n-1
          crd2%i=j
          point(2) = crd2%i*crd2%d + crd2%min
          call quat_mult(quaternionFromEulerAxis(point(2)*crd2%twist,(/0d0,1d0,0d0/)),qtwist(:,3),qtwist(:,2))
          do i=0, crd1%n-1
             crd1%i=i
             point(1) = crd1%i*crd1%d + crd1%min
             rotatepoint = point
             call quat_mult(quaternionFromEulerAxis(point(1)*crd1%twist,(/1d0,0d0,0d0/)),qtwist(:,2),qtwist(:,1))
             call quat_mult(orient_quat,qtwist(:,1),qfinal)
             call rotate_quat(rotatepoint,qfinal)
             rotatepoint=rotatepoint+origin
             grid_out(dim_out1%i+1,dim_out2%i+1,iguv) = &
                  grid_out(dim_out1%i+1,dim_out2%i+1,iguv)&
                  + interpVal(rotatepoint,guv)
          end do
       end do
    end do
  end subroutine averageCartesian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!calculates cylindrical coordinate orientational averages
!!!IN:
!!!   iguv : the sequential number of the guv file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine averageCylindrical(iguv)
    use rism_util, only : freeUnit
    use constants, only : pi
    implicit none
    integer, intent(in) :: iguv
    !point : point before rotation
    !rotatepoint : point after rotation
    !qtwist : an array of quaterions, one for each dimension, that are
    !         used to hold the twist rotation for this coordinate
    _REAL_ :: point(3), rotatepoint(3), qtwist(4,3), qfinal(4)
    integer :: i,j,k
    !the gfortran won't let me use %i as a counter so we use i,j and k
    !and then set %i
    do k = 0, crd3%n - 1
       crd3%i = k
       point(3) = crd3%i * crd3%d + crd3%min
       qtwist(:, 3) = quaternionFromEulerAxis(point(3) * crd3%twist, (/ 0d0, 0d0, 1d0 /))
       do j = 0, crd2%n - 1
          crd2%i = j
          do i = 0, crd1%n - 1
             crd1%i = i
             point(2) = (crd1%i * crd1%d + crd1%min) * sin(crd2%i * crd2%d + crd2%min)
             point(1) = (crd1%i * crd1%d + crd1%min) * cos(crd2%i * crd2%d + crd2%min)
             call quat_mult(&
                  quaternionFromEulerAxis(point(2) * crd2%twist, (/ 0d0, 1d0, 0d0 /)), &
                  qtwist(:, 3), qtwist(:, 2))
             call quat_mult(&
                  quaternionFromEulerAxis(point(1) * crd1%twist, (/ 1d0, 0d0, 0d0 /)), &
                  qtwist(:, 2), qtwist(:, 1))
             rotatepoint = point
             call quat_mult(orient_quat, qtwist(:, 1), qfinal)
             call rotate_quat(rotatepoint, qfinal)
             rotatepoint = rotatepoint + origin
             grid_out(dim_out1%i + 1, dim_out2%i + 1, iguv) = &
                  grid_out(dim_out1%i + 1, dim_out2%i + 1, iguv) &
                  + interpVal(rotatepoint, guv)
          end do
       end do
    end do
  end subroutine averageCylindrical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Determines if a given point is inside the grid
!!!IN:
!!!   point :: xyz point
!!!OUT:
!!!    .true. if it is in the grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function inGrid(point)
    implicit none
    _REAL_, intent(in) :: point(3)
    logical :: inGrid
    _REAL_ :: maxExtent(3)
    integer :: id
    inGrid = .true.

    maxExtent = dxOrigin + grid%boxLength - grid%spacing

    do id=1,3
       if(point(id) < dxOrigin(id) .or. point(id) > maxExtent(id))then
          inGrid = .false.
          return
       end if
    end do
    
  end function inGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Interpolates the value of point from the grid
!!!IN:
!!!   point :: xyz point
!!!   dist  :: distribution to interpolate
!!!OUT:
!!!    interpolated value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function interpVal(point,dist) result(val)
    implicit none
    _REAL_, intent(in) :: point(3), dist(:)
    _REAL_ :: val
    _REAL_ :: findex(3)
    !nb :: neighbour list
    integer :: nb(8,3), i,j, k
    val = 1d0
    !1) get fraction grid index
    findex = (point - dxOrigin)/grid%spacing+1d0
    !2) get eight nearest neighbours
    nb(1,:) = floor(findex) !x000
!    if(nb(1,3) == 1) write(0,*) "NB", NB(1,3), point(3), dxorigin(3), findex(3)
    nb(8,:) = nb(1,:)+1 !x111
!    if(nb(8,3) == grid%localDimsR(3)) write(0,*) "NB", NB(8,3), point(3), dxorigin(3), findex(3)
    !catches the case that the point is exactly on the upper boundary
    do i = 1,3
       if(nb(1,i) == grid%localDimsR(i) .and. floor(findex(i)) == ceiling(findex(i)))then
          nb(1,i) = nb(1,i)-1
          nb(8,i) = nb(8,i)-1
       end if
    end do
    !check that both are inbounds
    if(any(nb(1,:) < 1))then
       call rism_report_error("(a,3(f8.3))", "Interpolate lower out-of-bounds error for: ",point)
    end if
    if(nb(8,1) > grid%localDimsR(1) .or. nb(8,2) > grid%localDimsR(2) .or. nb(8,3) > grid%localDimsR(3))then
!       write(0,*) nb(8,:), grid%localDimsR, point(3), (floor(findex(3)) == ceiling(findex(3)))
       call rism_report_error("(a,3(f8.3))", "Interpolate upper out-of-bounds error for: ",point)
    end if

    nb(2,:) = (/nb(1,1),nb(1,2),nb(8,3)/) !x001
    nb(3,:) = (/nb(1,1),nb(8,2),nb(1,3)/) !x010
    nb(4,:) = (/nb(1,1),nb(8,2),nb(8,3)/) !x011
    nb(5,:) = (/nb(8,1),nb(1,2),nb(1,3)/) !x100
    nb(6,:) = (/nb(8,1),nb(1,2),nb(8,3)/) !x101
    nb(7,:) = (/nb(8,1),nb(8,2),nb(1,3)/) !x110
    do i = 1,8
       do j = i+1,8
          if(nb(i,1) == nb(j,1) .and.&
               nb(i,2) == nb(j,2) .and.&
               nb(i,3) == nb(j,3))then
             write(0,*) "neighbours ", i, " and ", j, " are the same:"
             write(0,*) nb(i,:)
             write(0,*) nb(j,:)
             write(0,*) point
             write(0,*) findex
          end if
       end do
    end do
    !turn fractional index to fractional position in voxel
    findex = findex-nb(1,:)
    if(any(findex > 1d0) .or. any(findex < 0d0))then
       write(0,*) "bad findex "
       write(0,*) findex
       write(0,*) point
       
    end if
    !interpolate
    call blend_103(findex(1),findex(2),findex(3),&
         dist(nb(1,1) + (nb(1,2)-1)*grid%localDimsR(1) + (nb(1,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(2,1) + (nb(2,2)-1)*grid%localDimsR(1) + (nb(2,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(3,1) + (nb(3,2)-1)*grid%localDimsR(1) + (nb(3,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(4,1) + (nb(4,2)-1)*grid%localDimsR(1) + (nb(4,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(5,1) + (nb(5,2)-1)*grid%localDimsR(1) + (nb(5,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(6,1) + (nb(6,2)-1)*grid%localDimsR(1) + (nb(6,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(7,1) + (nb(7,2)-1)*grid%localDimsR(1) + (nb(7,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         dist(nb(8,1) + (nb(8,2)-1)*grid%localDimsR(1) + (nb(8,3)-1)*grid%localDimsR(1)*grid%localDimsR(2)),&
         val)

  end function interpVal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!write a line to a pdbfile
!!!IN:
!!!   unit : unit to write to
!!!   serial : atom number
!!!   atmname : atom name
!!!   resname : residue name
!!!   chain   : chain name
!!!   resid   : residue number
!!!   xyz     : coordinates
!!!   beta_o  : (optional) beta value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writePDBLine(unit, serial,atmname,resname,chain,resid,xyz,beta_o)
    implicit none
    integer, intent(in) :: unit,serial,resid
    character(len=*), intent(in) :: atmname, resname,chain
    _REAL_, intent(in) :: xyz(3)
    _REAL_, optional, intent(in) :: beta_o
    _REAL_ :: beta
    beta = 0
    if(present(beta_o)) beta = beta_o
    write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,2A2)')&
         "ATOM  ",serial,atmname,"",resname,chain,resid,'',xyz,beta
  end subroutine writePDBLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeOutput
    use rism_util, only : freeUnit
    implicit none
    integer :: unit, iostat
    if(len_trim(outfile) > 0)then
       unit = freeUnit()
       open(unit,file=outfile,status='replace',iostat=iostat)
       if(iostat/=0)&
            call rism_report_error("failed to open :"//trim(outfile))
    else
       unit = 6
    end if
    if(ubound(grid_out,2)==1)then
       call write1DOutput(unit)
    else
       call write2DOutput(unit)
    end if
    if(len_trim(outfile) > 0)then
       close(unit)
    end if
  end subroutine writeOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print 1D radial results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write1DOutput(unit)
    use rism_util, only : rmExPrec
    implicit none
    integer, intent(in) :: unit
    character(len=80) :: fmt
    integer :: i, iv
    write(fmt,'(a,i4,a)') '(1p,',(ubound(grid_out,3)),'(e16.7E3,1x))'
    write(unit,'(a16,1x,a16)') "#"//trim(dim_out1%name), "value"
    do i=1,dim_out1%n
       write(unit,fmt,advance='no') (i-1)*dim_out1%d +dim_out1%min
       do iv = 1, ubound(grid_out,3)
          write(unit,fmt,advance='no') &
               rmExPrec(grid_out(i,1,iv))
       end do
       write(unit,'(a)')
    end do
  end subroutine write1DOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!print 2D radial results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write2DOutput(unit)
    use rism_util, only : rmExPrec
    implicit none
    integer, intent(in) :: unit
    character(len=80) :: fmt
    integer :: i, j
    write(fmt,'(a,i4,a)') '(1p,',(2+ubound(grid_out,3)),'(e16.7E3,1x))'
    write(unit,'(a16,1x,a16,1x,a16)') "#"//trim(dim_out1%name),&
         trim(dim_out2%name), "value"
    do i=1,dim_out1%n
       do j=1,dim_out2%n
          if(grid_out(1,j,1) /= huge(1d0))&
               write(unit,fmt) (i-1)*dim_out1%d + dim_out1%min,&
               (j-1)*dim_out2%d + dim_out2%min, rmExPrec(grid_out(i,j,:))
       end do
       write(unit,'(a)')
    end do
  end subroutine write2DOutput



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Prints input paramters to standard out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_parameters()
    implicit none
    integer :: unit, i
    integer,parameter :: titlewidth=23
    character(len=CLEN) :: whtspc, fmttxt, fmttxtvec
    unit = rism_report_getMUnit()
    write(whtspc,'(a)') " "
    write(fmttxt,'(a,i2,a)') '(a',titlewidth,',a)'
    write(fmttxtvec,'(a,i2,a)') '(1p,a',titlewidth,',3(e16.8))'
    write(unit,'(a)') "Calculation parameters:"
    write(unit,'(a)') "-----------------------"
    write(unit,fmttxt) "  Input density files:"//whtspc,trim(guvfile(1))
    do i=2,ubound(guvfile,1)
       write(unit,fmttxt) ""//whtspc,trim(guvfile(i))
    end do
    write(unit,fmttxt) "  Coordinate system:"//whtspc,trim(systemname(system))
    if(origin(1) == huge(1d0))then
       write(unit,fmttxt) "  Origin:"//whtspc,'<center>'
    else
       write(unit,fmttxtvec) "  Origin:"//whtspc,origin
    end if
    write(unit,fmttxtvec) "  z-vector:"//whtspc,vec
    write(unit,fmttxt) "  Output dimensions:"//whtspc,trim(dimensionnames(dim))
    write(unit,fmttxt) "  Abscissa:"//whtspc,trim(method)
    write(unit,fmttxt,advance='no') "  Output:"//whtspc
    if(len_trim(outfile) > 0)then
       write(unit,'(a)') trim(outfile)
    else
       write(unit,'(a)') '<stdout>'
    end if
    
  end subroutine print_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Prints grid specs to standard out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_gridspecs()
    implicit none
    integer :: unit, i
    integer,parameter :: titlewidth=23
    character(len=CLEN) :: whtspc, fmttxt, fmttxtreal, fmttxtint
    character(len=1) :: xyz(3)
    xyz(1)="x"
    xyz(2)="y"
    xyz(3)="z"
    unit = rism_report_getMUnit()
    write(whtspc,'(a)') ' '
    write(fmttxt,'(a,i2,a)') '(a',titlewidth,',a)'
    write(fmttxtint,'(a,i2,a)') '(a',titlewidth,',a5,a6,i8)'
    write(fmttxtreal,'(a,i2,a)') '(1p,a',titlewidth,',a5,a6,e16.8)'
    write(unit,'(a)') "Grid parameters:"
    write(unit,'(a)') "-----------------------"
    write(unit,fmttxt) "  Input Grid:"//whtspc
    do i=1,3
       write(unit,fmttxtreal) whtspc,xyz(i)//whtspc,"length"//whtspc,grid%boxLength(i)
       write(unit,fmttxtreal) whtspc,xyz(i)//whtspc,"d"//xyz(i)//whtspc,grid%spacing(i)
       write(unit,fmttxtreal) whtspc,xyz(i)//whtspc,"min"//whtspc,dxorigin(i)
       write(unit,fmttxtreal) whtspc,xyz(i)//whtspc,"max"//whtspc,dxorigin(i)+grid%boxLength(i)
       write(unit,fmttxtint) whtspc,xyz(i)//whtspc,"N"//whtspc,grid%localDimsR(i)
       write(unit,*)
    end do
    write(unit,fmttxt) "  Box:"//whtspc
    write(unit,'(a80)') dimstring(crd1,o_indent=titlewidth)
    write(unit,*)
    write(unit,'(a80)') dimstring(crd2,o_indent=titlewidth)
    write(unit,*)
    write(unit,'(a80)') dimstring(crd3,o_indent=titlewidth)
    write(unit,fmttxt) "  Output grid:"//whtspc
    write(unit,'(a80)') dimstring(dim_out1,o_indent=titlewidth)
    
    if(trim(dim_out2%name) .ne. 'dummy') then
       write(unit,*)
       write(unit,'(a80)') dimstring(dim_out2,o_indent=titlewidth)
    end if
  end subroutine print_gridspecs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a text expression for the coordinate system name
!!!IN:
!!!   system : system type that correspond to one of the defined constants
!!!OUT:
!!!   string representation of the coordinate system name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function systemname(system)
    implicit none
    integer, intent(in) :: system
    character(len=CLEN) :: systemname
    if(system == CARTESIAN) systemname='cartesian'
    if(system == CYLINDRICAL) systemname='cylindrical'
    if(system == SPHERICAL) systemname='SPHERICAL'
  end function systemname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a text expression for output dimensions.  Expects the
!!!output dimensions to be represented as a sum of dimension constants
!!!IN:
!!!   dim : dimension type that is the sum of one or two defined
!!!         constants. It is assumed that any dimension is only listed
!!!         once
!!!OUT:
!!!    string representation of the combination of dimension names
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dimensionnames(dim)
    implicit none
    integer, intent(in) :: dim
    character(len=CLEN) :: dimensionnames
    integer :: tdim,count
    dimensionnames=''
    count = 0
    tdim=dim
    if(tdim - XDIM >= 0)then
       tdim = tdim-XDIM
       dimensionnames = 'x'
       count = count+1
    end if
    if(tdim - YDIM >= 0)then
       tdim = tdim-YDIM
       if(count>0) dimensionnames = trim(dimensionnames)//'-'
       dimensionnames = trim(dimensionnames)//'y'
       count = count+1
    end if
    if(tdim - RDIM >= 0)then
       tdim = tdim-RDIM
       if(count>0) dimensionnames = trim(dimensionnames)//'-'
       dimensionnames = trim(dimensionnames)//'r'
       count = count+1
    end if
    if(tdim - THETADIM >= 0)then
       tdim = tdim-THETADIM
       if(count>0) dimensionnames = trim(dimensionnames)//'-'
       dimensionnames = trim(dimensionnames)//'theta'
       count = count+1
    end if
    if(tdim - PHIDIM >= 0)then
       tdim = tdim-PHIDIM
       if(count>0) dimensionnames = trim(dimensionnames)//'-'
       dimensionnames = trim(dimensionnames)//'phi'
       count = count+1
    end if
  end function dimensionnames

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a multiline string representation of the dimension object.
!!!Each line is optionally idented with whitespace
!!!IN:
!!!   dim : dimension object
!!!   o_indent :: (optional) number of characters of whitespace to indent
!!!OUT:
!!!    array string representation of the combination of dimension names
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function dimstring(dim,o_indent)
    implicit none
    type(dimension), intent(in) :: dim
    integer, optional, intent(in) :: o_indent
    character(len=CLEN) :: dimstring(6)
    character(len=CLEN) :: whtspc, fmttxt, fmttxtreal, fmttxtint
    integer :: indent
    dimstring=''
    indent=0
    write(whtspc,'(a)') " "
    if(present(o_indent)) indent=o_indent
    write(fmttxt,'(a,i2,a)') '(a',indent,',a)'
    write(fmttxtint,'(a,i2,a)') '(a',indent,',a6,a6,i8)'
    write(fmttxtreal,'(a,i2,a)') '(1p,a',indent,',a6,a6,e16.8)'
    write(dimstring(1),fmttxtreal) whtspc,dim%name//whtspc,"length"//whtspc,dim%len
    write(dimstring(2),fmttxtreal) whtspc,dim%name//whtspc,"d"//dim%name//whtspc,dim%d
    write(dimstring(3),fmttxtreal) whtspc,dim%name//whtspc,"min"//whtspc,dim%min
    write(dimstring(4),fmttxtreal) whtspc,dim%name//whtspc,"max"//whtspc,dim%max
    write(dimstring(5),fmttxtint) whtspc,dim%name//whtspc,"N"//whtspc,dim%n
    if(dim%applytwist)then
       write(dimstring(6),fmttxtreal) whtspc,dim%name//whtspc,"twist"//whtspc,dim%twist
    else
       write(dimstring(6),fmttxtreal) whtspc,dim%name//whtspc,"twist"//whtspc,0d0
    end if
    dimstring(1) = trim(dimstring(1))
  end function dimstring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Clean up memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine cleanup
    implicit none
    integer*8 :: memstats(10)
    integer :: unit , i
    unit = rism_report_getMUnit()
    call rism3d_grid_destroy(grid)
    if(safemem_dealloc(guv)/=0) call rism_report_error("Failed to deallocate guv")
    if(safemem_dealloc(guvfile)/=0) call rism_report_error("Failed to deallocate guvfile")
    if(safemem_dealloc(abscissa)/=0) call rism_report_error("Failed to deallocate abscissa")
    if(safemem_dealloc(grid_out)/=0) call rism_report_error("Failed to deallocate grid_out")
    if(verbose)then
       memstats = memStatus()
       write(unit,'(a)') "rism3d.orave memory allocation summary"
       write(unit,'(a)') "Type         Current         Maximum"
       write(unit,'(a,i12,a,f12.5,a)') "Integer  ",memstats(1)," B ",&
            dble(memstats(6))/BYTES_PER_GB," GB"
       write(unit,'(a,i12,a,f12.5,a)') "Real     ",memstats(2)," B ",&
            dble(memstats(7))/BYTES_PER_GB," GB"
       write(unit,'(a,i12,a,f12.5,a)') "Logical  ",memstats(3)," B ",&
            dble(memstats(8))/BYTES_PER_GB," GB"
       write(unit,'(a,i12,a,f12.5,a)') "Character",memstats(4)," B ",&
            dble(memstats(9))/BYTES_PER_GB," GB"
       write(unit,'(a)') "---------------------------------------"
       write(unit,'(a,i12,a,f12.5,a)') "Total    ",memstats(5)," B ",&
            dble(memstats(10))/BYTES_PER_GB," GB"
    end if
  end subroutine cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Usage description
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine usage
    implicit none
    integer :: munit
    munit = rism_report_getMUnit()
    call rism_report_setMUnit(rism_report_getEUnit())
    call rism_report_message("USAGE: rism3d.orave (-g|--guv) guvfiles ")
    call rism_report_message("                     --system (spherical|cylindrical|cartesian)")
    call rism_report_message("                     [-o|--origin x,y,z] [-v|--zvec x,y,z]")
    call rism_report_message("                     [-x] [-y] [-z] [-r] [-t|--theta] [-p|--phi]")
    call rism_report_message("                     [--lenx len] [--dx del] [--nx num]")
    call rism_report_message("                     [--leny len] [--dy del] [--ny num]")
    call rism_report_message("                     [--lenz len] [--dz del] [--nz num]")
    call rism_report_message("                     [--lenr len] [--dr del] [--nr num]")
    call rism_report_message("                     [--dtheta del] [--ntheta num]")
    call rism_report_message("                     [--dphi del] [--nphi num]")
    call rism_report_message("                     [--twist dim[,angle]")
    call rism_report_message("                     [--pdbout file] [--method name]")
    call rism_report_message("                     [-h|--help] [output file]")

    call rism_report_message("")
    call rism_report_message("Averages the distributions (--guv) given on the commandline in ")
    call rism_report_message("requested coordinate system (--system) about a selected center ")
    call rism_report_message("(--origin) with the z-axis rotated to a new position (--zvec). ")
    call rism_report_message("One or two dimensions for the coordinate system may be selected ")
    call rism_report_message("for output (-x -y -z -r -t -p). The the subvolume is defined by ")
    call rism_report_message("providing at least two of the length (len), spacing (del) or ")
    call rism_report_message("number of points (num) for each dimension. Theta always has ")
    call rism_report_message("length 2PI and phi has length PI. For spherical averaging, a ")
    call rism_report_message("an alternate sampling method for the angles may be selected ")
    call rism_report_message("(--method). Output is to standard out unless a file name is ")
    call rism_report_message("specified.")
    call rism_report_message("")
    call rism_report_message("--guv     : list of guv files to average.")
    call rism_report_message("--system  : spherical, cylindrical or cartesian coordinate system.")
    call rism_report_message("--origin  : (Default: center of guv) center for averaging ")
    call rism_report_message("            coodidinate system.")
    call rism_report_message("--zvec    : (Default: 0,0,1) rotate the guv z-axis to this ")
    call rism_report_message("            z-axis. The magnitude of this vector is not used.")
    call rism_report_message("-x        : output the x-dimension.")
    call rism_report_message("-y        : output the y-dimension.")
    call rism_report_message("-z        : output the z-dimension.")
    call rism_report_message("-r        : output the r-dimension.")
    call rism_report_message("--theta   : output the theta-dimension.")
    call rism_report_message("--phi     : output the phi-dimension.")
    call rism_report_message("--lenx    : length along x-axis.")
    call rism_report_message("--leny    : length along y-axis.")
    call rism_report_message("--lenz    : length along z-axis.")
    call rism_report_message("--lenr    : size of radius.")
    call rism_report_message("--dx      : x-dimension spacing.")
    call rism_report_message("--dy      : y-dimension spacing.")
    call rism_report_message("--dz      : z-dimension spacing.")
    call rism_report_message("--dr      : radial-dimension spacing")
    call rism_report_message("--dtheta  : theta-dimension spacing.")
    call rism_report_message("--dphi    : phi-dimension spacing.")
    call rism_report_message("--nx      : Number of sampling points for x.")
    call rism_report_message("--ny      : Number of sampling points for y.")
    call rism_report_message("--nz      : Number of sampling points for z.")
    call rism_report_message("--nr      : Number of sampling points for r.")
    call rism_report_message("--ntheta  : Number of sampling points for theta.")
    call rism_report_message("--nphi    : Number of sampling points for phi.")
    call rism_report_message("--twist   : Continuously rotate the data about the given ")
    call rism_report_message("            dimension (dim) by a given angle (ang0). Default angle")
    call rism_report_message("            is 0.18587; sufficient to unwind B DNA.")
    call rism_report_message("            ang(dim) = -ang0*dim")
    call rism_report_message("--pdbout  : Debugging. Considerably slows calculation. Outputs a")
    call rism_report_message("            PDB format file of all the points sampled. Points in")
    call rism_report_message("            the solvent box are atom type 'N'; those outside are")
    call rism_report_message("            of type 'O'.  Sampled values are placed in the occupancy")
    call rism_report_message("            column.")
    call rism_report_message("--method  : Use abscissa and weights from method")
    call rism_report_message("              LEB (Lebedev) (spherical)")
    call rism_report_message("              REP (REPULSION) (spherical)")
    call rism_report_message("              ZCW (spherical)")
    call rism_report_message("              SHREWD_REP (SHREWD + REPULSION) (spherical)")
    call rism_report_message("              SHREWD_ZCW (SHREWD + ZCW) (spherical)")
    call rism_report_message("              uniform (spherical or cylindrical)")
    call rism_report_message("            There are only a few sets of valid NP for each method")
    call rism_report_message("            except 'uniform'.")
    call rism_report_setMUnit(munit)
  end subroutine usage
end module orave_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program orave
  use orave_m
  implicit none

  call getOptions()
  call readDistributionFiles()
  call setAbscissa()
  call average()
  call writeOutput()
  call cleanup()

end program orave

