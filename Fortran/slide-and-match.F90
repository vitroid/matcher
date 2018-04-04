!Modified points from slide-and-overlay:
!  It reads two files. First is for the template, and the second is for the lattice to be slided.

subroutine readAR3A(FILE, x,y,z, nmol)
  implicit none
  integer :: FILE
  real(kind=8), intent(INOUT) :: x(*), y(*), z(*)
  integer, intent(OUT) :: nmol
  !   local
  integer  :: i
  read(FILE,*) nmol
  do i=1, nmol
     read(FILE,*) x(i),y(i),z(i)
  enddo
end subroutine readAR3A


!trim the list of molecules within the radius
subroutine trimall(rad,x,y,z,nin,nout)
  implicit none
  real(kind=8) :: rad
  real(kind=8), intent(INOUT) :: x(*), y(*), z(*)
  integer, intent(IN) :: nin
  integer, intent(OUT) :: nout
  !local
  integer :: head
  real(kind=8) :: rr
  head = 1
  nout = 0
  do while ( head <= nin )
     rr = x(head)**2 + y(head)**2 + z(head)**2
     if ( rr < rad**2 ) then
        nout = nout + 1
        x(nout) = x(head)
        y(nout) = y(head)
        z(nout) = z(head)
     endif
     head = head + 1
  enddo
end subroutine trimall



!put all in the image cell
subroutine wrapall(nmol,x,y,z,bx,by,bz)
  implicit none
  integer, intent(IN) :: nmol
  real(kind=8), intent(INOUT) :: x(*), y(*), z(*)
  real(kind=8), intent(IN)    :: bx, by, bz
  !local
  integer :: i
  do i=1, nmol
     x(i) = x(i) - dnint( x(i) / bx ) * bx
     y(i) = y(i) - dnint( y(i) / by ) * by
     z(i) = z(i) - dnint( z(i) / bz ) * bz
  enddo
end subroutine wrapall



!copy all
subroutine copyall(nmol,x,y,z,tx,ty,tz)
  implicit none
  integer, intent(IN) :: nmol
  real(kind=8), intent(IN) :: x(*), y(*), z(*)
  real(kind=8), intent(OUT) :: tx(*), ty(*), tz(*)
  !local
  integer :: i
  do i=1, nmol
     tx(i) = x(i)
     ty(i) = y(i)
     tz(i) = z(i)
  enddo
end subroutine copyall


!give offset to all
subroutine shiftall(nmol,x,y,z,dx,dy,dz)
  implicit none
  integer, intent(IN) :: nmol
  real(kind=8), intent(INOUT) :: x(*), y(*), z(*)
  real(kind=8), intent(IN) :: dx,dy,dz
  !local
  integer :: i
  do i=1, nmol
     x(i) = x(i) + dx
     y(i) = y(i) + dy
     z(i) = z(i) + dz
  enddo
end subroutine shiftall



!look up the nearest neighbor
function trial0(nr,xr,yr,zr,ns,xs,ys,zs)
  implicit none
  real(kind=8) :: trial0
  integer, intent(IN) :: nr,ns
  real(kind=8), intent(IN) :: xr(*), yr(*), zr(*)
  real(kind=8), intent(IN) :: xs(*), ys(*), zs(*)
  !local
  integer :: ir,is
  real(kind=8) :: dmin,dx,dy,dz,d,score
  score = 0.0
  do ir=1,nr
     dmin = 99999.9d0
     do is=1,ns
        dx = xr(ir) - xs(is)
        dy = yr(ir) - ys(is)
        dz = zr(ir) - zs(is)
        d = dx**2 + dy**2 + dz**2
        if ( d < dmin ) then
           dmin = d
        endif
     enddo
     score = score + dmin
  enddo
  trial0 = score
end function trial0



!look up the nearest neighbor
function trial(nr,xr,yr,zr,ns,xs,ys,zs,vx,vy,vz)
  implicit none
  real(kind=8) :: trial
  integer, intent(IN) :: nr,ns
  real(kind=8), intent(IN) :: xr(*), yr(*), zr(*)
  real(kind=8), intent(IN) :: xs(*), ys(*), zs(*)
  real(kind=8), intent(OUT):: vx,vy,vz
  !local
  integer :: ir,is
  real(kind=8) :: dmin,dx,dy,dz,d,score,xmin,ymin,zmin
  integer :: nearest(nr), smin
  score = 0.0
  vx = 0d0
  vy = 0d0
  vz = 0d0  
  do ir=1,nr
     dmin = 99999.9d0
    do is=1,ns
        dx = xr(ir) - xs(is)
        dy = yr(ir) - ys(is)
        dz = zr(ir) - zs(is)
        d = dx**2 + dy**2 + dz**2
        if ( d < dmin ) then
           dmin = d
           xmin = dx
           ymin = dy
           zmin = dz
           smin = is
        endif
     enddo
     score = score + dmin
     vx = vx + xmin
     vy = vy + ymin
     vz = vz + zmin
     nearest(ir) = smin
  enddo
  !write(0,*) score, "before"
  vx = vx / nr
  vy = vy / nr
  vz = vz / nr
  score = 0.0
  do ir=1,nr
     dmin = 99999.9d0
    do is=1,ns
        dx = xr(ir) - xs(is) - vx
        dy = yr(ir) - ys(is) - vy
        dz = zr(ir) - zs(is) - vz
        d = dx**2 + dy**2 + dz**2
        if ( d < dmin ) then
           dmin = d
        endif
     enddo
     score = score + dmin
  enddo  
  !write(0,*) score, "after"
  trial = score
end function trial




program main
  use mt95, only: genrand_init, genrand_real2
  implicit none
  !declaration
  real(kind=8) :: trial
  !constants
  integer, parameter :: MAXMOL = 50000
  integer, parameter :: STDIN=5, STDOUT=6, STDERR=0
  !variables
  integer :: FILE
  character(len=5) :: tag
  !reference lattice
  integer  :: nr
  real(kind=8) :: xr(MAXMOL), yr(MAXMOL), zr(MAXMOL)
  !original lattice
  integer  :: no
  real(kind=8) :: xo(MAXMOL), yo(MAXMOL), zo(MAXMOL)
  !shifted lattice
  integer  :: ns
  real(kind=8) :: xs(MAXMOL), ys(MAXMOL), zs(MAXMOL)
  !offsets
  real(kind=8) :: shiftx, shifty, shiftz
  real(kind=8) :: bestx, besty, bestz
  real(kind=8) :: width, bestscore
  real(kind=8) :: dx,dy,dz
  real(kind=8) :: vx,vy,vz
  integer      :: loop,l
  !box size
  real(kind=8) :: bx, by, bz
  real(kind=8) :: score
  !file names
  character(len=1024) :: template,lattice,loopstr
  !specify two files at the command line
  call getarg(1,loopstr)
  call getarg(2,template)
  call getarg(3,lattice)
  read(loopstr,'(i10)') loop
  FILE=10
99   format(a5)
  open(FILE,file=template)
  do
     read(FILE,99) tag
     if ( tag .eq. "@AR3A" ) then
        call readAR3A(FILE, xr,yr,zr, nr)
        exit
     else
        print tag
     endif
  enddo
  close(10)
  open(FILE,file=lattice)
  do
     read(FILE,99) tag
     if ( tag .eq. "@BOX3" ) then
        read(FILE,*) bx,by,bz
        !write(STDERR,*) bx,by,bz, "BOX"
     else if ( tag .eq. "@AR3A" ) then
        call readAR3A(FILE, xo,yo,zo, no)
        exit
     else
        print tag
     endif
  enddo

  ! trial
  call genrand_init

  !infinite loop
  do l=1,loop
     width = 2.0
     bestscore = 300.0
     call genrand_real2(dx)
     call genrand_real2(dy)
     call genrand_real2(dz)
     bestx = (dx-0.5)*bx
     besty = (dy-0.5)*by
     bestz = (dz-0.5)*bz
     do while ( width > 0.1 )
        call genrand_real2(dx)
        call genrand_real2(dy)
        call genrand_real2(dz)
        shiftx = bestx + (dx*2-1)*width
        shifty = besty + (dy*2-1)*width
        shiftz = bestz + (dz*2-1)*width
        call copyall(no, xo,yo,zo, xs,ys,zs)
        call shiftall(no, xs,ys,zs, shiftx, shifty, shiftz)
        call wrapall(no,xs,ys,zs,bx,by,bz)
        call trimall(8.0d0+3.0, xs,ys,zs, no, ns)
        !write(STDERR,*) nmol,nr,ns
        score = trial(nr,xr,yr,zr,ns,xs,ys,zs,vx,vy,vz)
        if ( score < bestscore ) then
           !write(STDERR,*) shiftx, shifty, shiftz, score, width
           bestscore = score
           shiftx = shiftx + vx
           shifty = shifty + vy
           shiftz = shiftz + vz
           bestx = shiftx
           besty = shifty
           bestz = shiftz
        endif
        width = width * 0.99
     end do
     write(STDOUT,*) shiftx, shifty, shiftz, score !, width
  enddo
end program main
  
