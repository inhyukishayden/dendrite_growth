!=============================================================
! mc.f90  (Fortran 90/2008 + OpenMP)
!
! - implicit none + explicit kinds (real64)
! - allocatable arrays sized to what is actually used
! - standard RNG via random_number/random_seed
! - OpenMP parallelization on major kernels (field solve + gradients + init)
! - fixes critical bugs in original code:
!     (1) box dimensions set before vanode scaling
!     (2) j1/j2 properly defined in gradient loop
!     (3) verr is OpenMP-reduction (no race)
!
! Notes:
! - The Monte Carlo particle-growth loop is kept SERIAL for correctness
!   (it mutates shared cell lists and has early-exit logic).
! - The Laplace relaxation is still Jacobi (as original), but OpenMP-safe.
!
! Compile (gfortran):
!   gfortran -O3 -fopenmp -march=native -ffast-math mc.f90 -o mc
!
! Compile (ifx):
!   ifx -O3 -qopenmp -xHost mc.f90 -o mc
!=============================================================

program mc
  use iso_fortran_env, only : dp => real64
  use omp_lib
  implicit none

  !----------------------------
  ! Problem/control parameters
  !----------------------------
  integer, parameter :: necellx = 100, necelly = 100, necellz = 300
  integer, parameter :: npcellx = 10,  npcelly = 10,  npcellz = 30
  integer, parameter :: maxN    = 2000000
  integer, parameter :: maxCell = 3001   ! max particles per particle-cell (1st dim of xcell/ycell/zcell)

  integer, parameter :: nzV  = 3*necellz + 2   ! v has boundaries at 1 and necellz+2, plus extended region
  integer, parameter :: nzE  = 3*necellz       ! vex/vey/vez stored for k=1..3*necellz

  !----------------------------
  ! I/O units
  !----------------------------
  integer, parameter :: utraj=111, ucoord=301, uef=302, upot=303, umax=305, ucond=501

  !----------------------------
  ! Scalars
  !----------------------------
  integer  :: nsd, maxnf
  real(dp) :: vanode, dt, rdrate, vaprob

  real(dp) :: boxx, boxy, boxz, pi
  real(dp) :: cellx, celly, cellz, csqx, csqy, csqz
  real(dp) :: sig, datt, rmove

  integer  :: nf, nfold, ncount
  real(dp) :: xnew, ynew, znew
  real(dp) :: xmove, ymove, zmove
  real(dp) :: exmove, eymove, ezmove
  integer  :: nx, ny, nz
  integer  :: npx, npy, npz
  integer  :: ncx, ncy, ncz
  real(dp) :: prob
  real(dp) :: start_time, end_time

  ! Laplace/Jacobi
  integer  :: icount
  real(dp) :: verr, denom

  ! temps in collision checks
  integer  :: ix, iy, iz, ix1, iy1
  integer  :: ipar
  real(dp) :: dx, dy, dz, dist

  !----------------------------
  ! Big arrays
  !----------------------------
  real(dp), allocatable :: x(:), y(:), z(:)
  real(dp), allocatable :: v(:,:,:), vnew(:,:,:)
  real(dp), allocatable :: vex(:,:,:), vey(:,:,:), vez(:,:,:)
  integer,  allocatable :: ncell(:,:,:)
  real(dp), allocatable :: xcell(:,:,:,:), ycell(:,:,:,:), zcell(:,:,:,:)

  !----------------------------
  ! Open files and read input
  !----------------------------
  open(unit=utraj,  file='traj_3d.gro',      status='unknown')
  open(unit=ucoord, file='coord_3d.dat',     status='unknown')
  open(unit=uef,    file='efield_3d.dat',    status='unknown')
  open(unit=upot,   file='epot_3d.dat',      status='unknown')
  open(unit=umax,   file='maxheight_3d.dat', status='unknown')
  open(unit=ucond,  file='cond.dat',         status='old')

  read(ucond,*) nsd
  read(ucond,*) maxnf
  read(ucond,*) vanode
  read(ucond,*) nuintv
  read(ucond,*) dt
  read(ucond,*) rdrate
  read(ucond,*) vaprob

  call cpu_time(start_time)

  !----------------------------
  ! Define constants (IMPORTANT: before using vanode scaling)
  !----------------------------
  boxx = 166.7_dp
  boxy = 166.7_dp
  boxz = 500.1_dp
  pi   = acos(-1.0_dp)

  ! Your original had: vanode = vanode*dble(boxz/boxx) THEN later again uses vanode*(boxz/boxx).
  ! That double-scales. I preserve your *stated* scaling ONCE here.
  vanode = vanode * (boxz/boxx)

  cellx = 2.0_dp*boxx / real(necellx,dp)
  celly = 2.0_dp*boxy / real(necelly,dp)
  cellz = 2.0_dp*boxz / real(necellz,dp)
  csqx  = cellx*cellx
  csqy  = celly*celly
  csqz  = cellz*cellz

  sig   = 1.2_dp
  datt  = 0.1_dp
  rmove = 1.67_dp*sqrt(rdrate)

  denom = (2.0_dp/csqx + 2.0_dp/csqy + 2.0_dp/csqz)

  !----------------------------
  ! Allocate
  !----------------------------
  allocate(x(maxN), y(maxN), z(maxN))
  allocate(v(necellx,necelly,nzV), vnew(necellx,necelly,nzV))
  allocate(vex(necellx,necelly,nzE), vey(necellx,necelly,nzE), vez(necellx,necelly,nzE))
  allocate(ncell(npcellx,npcelly,npcellz))
  allocate(xcell(maxCell,npcellx,npcelly,npcellz))
  allocate(ycell(maxCell,npcellx,npcelly,npcellz))
  allocate(zcell(maxCell,npcellx,npcelly,npcellz))

  call init_rng(nsd)

  call init_arrays(x,y,z, v,vnew, vex,vey,vez, ncell, xcell,ycell,zcell)

  write(*,*) 'INIT DONE'

  !----------------------------
  ! Initialize electrode potential profile in v
  ! v(:,:,1) = 0, v(:,:,necellz+2)=vanode, linear ramp in between (k=2..necellz+1)
  !----------------------------
  call init_potential(v, vanode, necellz)

  !----------------------------
  ! Initial walker position
  !----------------------------
  xnew   = urand()*boxx
  ynew   = urand()*boxy
  znew   = boxz
  nf     = 1
  nfold  = 0
  ncount = 1

  !============================
  ! Main growth loop (SERIAL)
  !============================
  do while (nf <= maxnf)

    !-------------------------------------------------------
    ! Electric field recompute occasionally (your mod(nf,10)==1)
    !-------------------------------------------------------
    if (nfold /= nf) then
      nfold = nfold + 1

      if (mod(nf,nuintv) == 1) then
        icount = 1

        ! reset vnew
        !$omp parallel do collapse(3) schedule(static)
        do nx=1,necellx
          do ny=1,necelly
            do nz=1,necellz+2
              vnew(nx,ny,nz) = 0.0_dp
            enddo
          enddo
        enddo
        !$omp end parallel do

        ! re-apply boundary conditions and initial ramp
        call init_potential(v, vanode, necellz)

        ! Jacobi relaxation loop
        do
          verr = 0.0_dp

          call jacobi_step(v, vnew, csqx, csqy, csqz, denom, necellx, necelly, necellz)

          ! swap (copy) and compute residual (OpenMP reduction)
          !$omp parallel do collapse(3) reduction(+:verr) schedule(static)
          do nx=1,necellx
            do ny=1,necelly
              do nz=2,necellz+1
                verr = verr + (v(nx,ny,nz) - vnew(nx,ny,nz))**2
                v(nx,ny,nz) = vnew(nx,ny,nz)
              enddo
            enddo
          enddo
          !$omp end parallel do

          ! Conductor constraint: set potential to 0 where particles exist (serial)
          if (nf > 1) then
            do ipar = 1, nf-1
              nx = int(floor(x(ipar)/(boxx/real(necellx,dp)))) + 1
              ny = int(floor(y(ipar)/(boxy/real(necelly,dp)))) + 1
              nz = int(floor(z(ipar)/(boxz/real(necellz,dp)))) + 2   ! matches your original +2
              nx = clamp_int(nx, 1, necellx)
              ny = clamp_int(ny, 1, necelly)
              nz = clamp_int(nz, 1, necellz+2)
              v(nx,ny,nz) = 0.0_dp
            enddo
          endif

          if (verr < 1.0e-9_dp) exit
          icount = icount + 1
          if (mod(icount,10) == 0) write(*,'(I8,1X,ES14.6)') icount, verr
        enddo

        ! Compute gradients -> vex/vey/vez (OpenMP)
        call compute_efield(v, vex, vey, vez, cellx, celly, cellz, necellx, necelly, necellz)

        ! Extend field above top (k=necellz+1..3*necellz) with top value (OpenMP)
        !$omp parallel do collapse(3) schedule(static)
        do nx=1,necellx
          do ny=1,necelly
            do nz=necellz+1, 3*necellz
              vex(nx,ny,nz) = vex(nx,ny,necellz)
              vey(nx,ny,nz) = vey(nx,ny,necellz)
              vez(nx,ny,nz) = vez(nx,ny,necellz)
            enddo
          enddo
        enddo
        !$omp end parallel do

      endif
    endif

    !=======================================================
    !  Random move
    !=======================================================
101 continue
    call random_displacement(rmove, xmove, ymove, zmove)

    nx = int(floor(xnew/(boxx/real(necellx,dp)))) + 1
    ny = int(floor(ynew/(boxy/real(necelly,dp)))) + 1
    nz = int(floor(znew/(boxz/real(necellz,dp)))) + 1
    nx = clamp_int(nx, 1, necellx)
    ny = clamp_int(ny, 1, necelly)
    nz = clamp_int(nz, 1, 3*necellz)

    exmove = 56.0_dp * vex(nx,ny,nz)
    eymove = 56.0_dp * vey(nx,ny,nz)
    ezmove = 56.0_dp * vez(nx,ny,nz)

    xnew = xnew + xmove*sqrt(dt) + exmove*dt
    xnew = xnew - boxx*floor(xnew/boxx)   ! periodic
    ynew = ynew + ymove*sqrt(dt) + eymove*dt
    ynew = ynew - boxy*floor(ynew/boxy)   ! periodic
    znew = znew + zmove*sqrt(dt) + ezmove*dt

    if (znew < 0.0_dp) goto 101
    if (znew > boxz) then
      znew = boxz
      goto 101
    endif

    ! Particle-cell index (for neighbor search)
    npx = int(floor(xnew/(boxx/real(npcellx,dp)))) + 1
    npy = int(floor(ynew/(boxy/real(npcelly,dp)))) + 1
    npz = int(floor(znew/(boxz/real(npcellz,dp)))) + 1
    npx = clamp_int(npx, 1, npcellx)
    npy = clamp_int(npy, 1, npcelly)
    npz = clamp_int(npz, 1, npcellz)

    !=======================================================
    !  Collision rejection against existing particles
    !=======================================================
    if (nf > 1) then
      do ix = npx-1, npx+1
        ix1 = ix
        if (ix1 == 0) ix1 = npcellx
        if (ix1 == npcellx+1) ix1 = 1

        do iy = npy-1, npy+1
          iy1 = iy
          if (iy1 == 0) iy1 = npcelly
          if (iy1 == npcelly+1) iy1 = 1

          do iz = npz-1, npz+1
            if (iz < 1 .or. iz > npcellz) cycle   ! SAFETY (your original could access iz=0 or 31)

            do ipar = 1, ncell(ix1,iy1,iz)
              dx = xnew - xcell(ipar,ix1,iy1,iz)
              dx = dx - boxx*nint(dx/boxx)
              dy = ynew - ycell(ipar,ix1,iy1,iz)
              dy = dy - boxy*nint(dy/boxy)
              dz = znew - zcell(ipar,ix1,iy1,iz)
              dist = sqrt(dx*dx + dy*dy + dz*dz)
              if (dist < sig) goto 101
            enddo
          enddo
        enddo
      enddo
    endif

    !=======================================================
    !  Deposition logic (two cases)                        
    !=======================================================

    ! Case 1: near bottom
    if (znew < sig + 0.1_dp) then
      prob = urand()
      if (prob > vaprob) then
        goto 101
      else
        call deposit_particle(nf, xnew,ynew,znew, x,y,z, ncell,xcell,ycell,zcell, boxx,boxy,boxz, ucoord,umax, ncount)
        xnew = urand()*boxx
        ynew = urand()*boxy
        znew = boxz
        nf = nf + 1
        goto 102
      endif
    endif

    ! Case 2: within (sig, sig+datt) of an existing particle (adhesion)
    do ix = npx-1, npx+1
      ix1 = ix
      if (ix1 == 0) ix1 = npcellx
      if (ix1 == npcellx+1) ix1 = 1

      do iy = npy-1, npy+1
        iy1 = iy
        if (iy1 == 0) iy1 = npcelly
        if (iy1 == npcelly+1) iy1 = 1

        do iz = npz-1, npz+1
          if (iz < 1 .or. iz > npcellz) cycle

          do ipar = 1, ncell(ix1,iy1,iz)
            dx = xnew - xcell(ipar,ix1,iy1,iz)
            dx = dx - boxx*nint(dx/boxx)
            dy = ynew - ycell(ipar,ix1,iy1,iz)
            dy = dy - boxy*nint(dy/boxy)
            dz = znew - zcell(ipar,ix1,iy1,iz)
            dist = sqrt(dx*dx + dy*dy + dz*dz)

            if (dist > sig .and. dist < sig + datt) then
              prob = urand()
              if (prob > vaprob) then
                goto 101
              else
                call deposit_particle(nf, xnew,ynew,znew, x,y,z, ncell,xcell,ycell,zcell, boxx,boxy,boxz, ucoord,umax, ncount)
                xnew = urand()*boxx
                ynew = urand()*boxy
                znew = boxz
                nf = nf + 1
                goto 102
              endif
            endif
          enddo
        enddo
      enddo
    enddo

102 continue
    ncount = ncount + 1

  enddo  ! nf loop

  !----------------------------
  ! Output efield and potential
  !----------------------------
  do nx=1,necellx
    do ny=1,necelly
      do nz=1,necellz
        call write_field_cell(uef, upot, nx,ny,nz, vex,vey,vez, v)
      enddo
    enddo
  enddo

  call cpu_time(end_time)
  write(*,'(A,F12.3)') 'CPU time (s): ', end_time-start_time

contains

  subroutine init_rng(nsd)
    implicit none
    integer, intent(in) :: nsd
    integer :: n, i
    integer, allocatable :: seed(:)
    call random_seed(size=n)
    allocate(seed(n))
    do i=1,n
      seed(i) = ieor(1234567 + 104729*i, nsd)
    enddo
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine init_rng

  real(dp) function urand()
    implicit none
    call random_number(urand)
  end function urand

  integer function clamp_int(val, lo, hi)
    implicit none
    integer, intent(in) :: val, lo, hi
    clamp_int = val
    if (clamp_int < lo) clamp_int = lo
    if (clamp_int > hi) clamp_int = hi
  end function clamp_int

  subroutine init_arrays(x,y,z, v,vnew, vex,vey,vez, ncell, xcell,ycell,zcell)
    implicit none
    real(dp), intent(inout) :: x(:),y(:),z(:)
    real(dp), intent(inout) :: v(:,:,:), vnew(:,:,:)
    real(dp), intent(inout) :: vex(:,:,:), vey(:,:,:), vez(:,:,:)
    integer,  intent(inout) :: ncell(:,:,:)
    real(dp), intent(inout) :: xcell(:,:,:,:), ycell(:,:,:,:), zcell(:,:,:,:)
    integer :: i,j,k,l

    !$omp parallel do schedule(static)
    do i=1,size(x)
      x(i)=0.0_dp; y(i)=0.0_dp; z(i)=0.0_dp
    enddo
    !$omp end parallel do

    !$omp parallel do collapse(3) schedule(static)
    do i=1,size(v,1)
      do j=1,size(v,2)
        do k=1,size(v,3)
          v(i,j,k)=0.0_dp
        enddo
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do collapse(3) schedule(static)
    do i=1,size(vnew,1)
      do j=1,size(vnew,2)
        do k=1,size(vnew,3)
          vnew(i,j,k)=0.0_dp
        enddo
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do collapse(3) schedule(static)
    do i=1,size(vex,1)
      do j=1,size(vex,2)
        do k=1,size(vex,3)
          vex(i,j,k)=0.0_dp
          vey(i,j,k)=0.0_dp
          vez(i,j,k)=0.0_dp
        enddo
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do collapse(3) schedule(static)
    do i=1,size(ncell,1)
      do j=1,size(ncell,2)
        do k=1,size(ncell,3)
          ncell(i,j,k)=0
        enddo
      enddo
    enddo
    !$omp end parallel do

    !$omp parallel do collapse(4) schedule(static)
    do i=1,size(xcell,1)
      do j=1,size(xcell,2)
        do k=1,size(xcell,3)
          do l=1,size(xcell,4)
            xcell(i,j,k,l)=0.0_dp
            ycell(i,j,k,l)=0.0_dp
            zcell(i,j,k,l)=0.0_dp
          enddo
        enddo
      enddo
    enddo
    !$omp end parallel do
  end subroutine init_arrays

  subroutine init_potential(v, vanode, necz)
    implicit none
    real(dp), intent(inout) :: v(:,:,:)
    real(dp), intent(in)    :: vanode
    integer,  intent(in)    :: necz
    integer :: i,j,k

    ! bottom electrode
    !$omp parallel do collapse(2) schedule(static)
    do i=1,size(v,1)
      do j=1,size(v,2)
        v(i,j,1) = 0.0_dp
        v(i,j,necz+2) = vanode
      enddo
    enddo
    !$omp end parallel do

    ! linear ramp
    !$omp parallel do collapse(3) schedule(static)
    do i=1,size(v,1)
      do j=1,size(v,2)
        do k=2,necz+1
          v(i,j,k) = vanode/real(necz,dp) * real(k-1,dp)
        enddo
      enddo
    enddo
    !$omp end parallel do
  end subroutine init_potential

  subroutine jacobi_step(v, vnew, csqx, csqy, csqz, denom, nex, ney, nez)
    implicit none
    real(dp), intent(in)    :: v(:,:,:)
    real(dp), intent(inout) :: vnew(:,:,:)
    real(dp), intent(in)    :: csqx, csqy, csqz, denom
    integer,  intent(in)    :: nex, ney, nez
    integer :: i,j,k, i1,i2,j1,j2

    !$omp parallel do collapse(2) private(i1,i2,j1,j2,k) schedule(static)
    do i=1,nex
      i1 = i+1; if (i1 == nex+1) i1 = 1
      i2 = i-1; if (i2 == 0)     i2 = nex
      do j=1,ney
        j1 = j+1; if (j1 == ney+1) j1 = 1
        j2 = j-1; if (j2 == 0)     j2 = ney
        do k=2,nez+1
          vnew(i,j,k) = ( (v(i1,j,k)+v(i2,j,k))/csqx &
                        + (v(i,j1,k)+v(i,j2,k))/csqy &
                        + (v(i,j,k+1)+v(i,j,k-1))/csqz ) / denom
        enddo
      enddo
    enddo
    !$omp end parallel do
  end subroutine jacobi_step

  subroutine compute_efield(v, vex, vey, vez, cellx, celly, cellz, nex, ney, nez)
    implicit none
    real(dp), intent(in)    :: v(:,:,:)
    real(dp), intent(inout) :: vex(:,:,:), vey(:,:,:), vez(:,:,:)
    real(dp), intent(in)    :: cellx, celly, cellz
    integer,  intent(in)    :: nex, ney, nez
    integer :: i,j,k, i1,i2,j1,j2
    real(dp) :: ex,ey,ez

    ! compute for k=1..nez from potential indices (k=2..nez+1)
    !$omp parallel do collapse(2) private(i1,i2,j1,j2,k,ex,ey,ez) schedule(static)
    do i=1,nex
      i1 = i+1; if (i1 == nex+1) i1 = 1
      i2 = i-1; if (i2 == 0)     i2 = nex
      do j=1,ney
        j1 = j+1; if (j1 == ney+1) j1 = 1
        j2 = j-1; if (j2 == 0)     j2 = ney
        do k=2,nez+1
          ex = -(v(i1,j,k)-v(i2,j,k))/cellx
          ey = -(v(i,j1,k)-v(i,j2,k))/celly
          ez = -(v(i,j,k+1)-v(i,j,k-1))/cellz
          vex(i,j,k-1) = ex
          vey(i,j,k-1) = ey
          vez(i,j,k-1) = ez
        enddo
      enddo
    enddo
    !$omp end parallel do
  end subroutine compute_efield

  subroutine random_displacement(rmove, xmove, ymove, zmove)
    implicit none
    real(dp), intent(in)  :: rmove
    real(dp), intent(out) :: xmove, ymove, zmove
    real(dp) :: b1, b2, bsq, bh

    do
      b1  = 1.0_dp - 2.0_dp*urand()
      b2  = 1.0_dp - 2.0_dp*urand()
      bsq = b1*b1 + b2*b2
      if (bsq <= 1.0_dp) exit
    enddo
    bh    = sqrt(1.0_dp - bsq)
    xmove = rmove * 2.0_dp*b1*bh
    ymove = rmove * 2.0_dp*b2*bh
    zmove = rmove * (1.0_dp - 2.0_dp*bsq)
  end subroutine random_displacement

  subroutine deposit_particle(nf, xnew,ynew,znew, x,y,z, ncell,xcell,ycell,zcell, boxx,boxy,boxz, ucoord, umax, ncount)
    implicit none
    integer,  intent(in)    :: nf, ucoord, umax, ncount
    real(dp), intent(in)    :: xnew,ynew,znew, boxx,boxy,boxz
    real(dp), intent(inout) :: x(:),y(:),z(:)
    integer,  intent(inout) :: ncell(:,:,:)
    real(dp), intent(inout) :: xcell(:,:,:,:), ycell(:,:,:,:), zcell(:,:,:,:)

    integer :: ncx,ncy,ncz, idx
    real(dp) :: zmax

    x(nf) = xnew
    y(nf) = ynew
    z(nf) = znew

    write(ucoord,'(I10,3F14.6)') nf, x(nf), y(nf), z(nf)
    write(*,'(I10,1X,I12)') nf, ncount

    zmax = maxval(z(1:nf))
    write(umax,'(I10,1X,F14.6)') nf, zmax

    ncx = int(floor(x(nf)/(boxx/real(size(ncell,1),dp)))) + 1
    ncy = int(floor(y(nf)/(boxy/real(size(ncell,2),dp)))) + 1
    ncz = int(floor(z(nf)/(boxz/real(size(ncell,3),dp)))) + 1
    ncx = clamp_int(ncx, 1, size(ncell,1))
    ncy = clamp_int(ncy, 1, size(ncell,2))
    ncz = clamp_int(ncz, 1, size(ncell,3))

    idx = ncell(ncx,ncy,ncz) + 1
    if (idx > size(xcell,1)) then
      write(*,*) 'ERROR: particle-cell overflow at (',ncx,ncy,ncz,') idx=',idx
      stop 2
    endif
    ncell(ncx,ncy,ncz) = idx
    xcell(idx,ncx,ncy,ncz) = x(nf)
    ycell(idx,ncx,ncy,ncz) = y(nf)
    zcell(idx,ncx,ncy,ncz) = z(nf)
  end subroutine deposit_particle

  subroutine write_field_cell(uef, upot, i,j,k, vex,vey,vez, v)
    implicit none
    integer, intent(in) :: uef, upot, i,j,k
    real(dp), intent(in) :: vex(:,:,:), vey(:,:,:), vez(:,:,:), v(:,:,:)
    real(dp) :: eval, exn, eyn, ezn

    eval = sqrt(vex(i,j,k)*vex(i,j,k) + vey(i,j,k)*vey(i,j,k) + vez(i,j,k)*vez(i,j,k))
    if (eval > 0.0_dp) then
      exn = vex(i,j,k)/eval
      eyn = vey(i,j,k)/eval
      ezn = vez(i,j,k)/eval
    else
      exn = 0.0_dp; eyn = 0.0_dp; ezn = 0.0_dp
    endif

    write(uef,'(3I6,1X,4F14.6)') i,j,k, exn, eyn, ezn, eval
    write(upot,'(3I6,1X,F14.6)') i,j,k, v(i,j,k)
  end subroutine write_field_cell

end program mc

