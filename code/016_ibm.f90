!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module ibm
!-----------------------------------------------------------------------
! immersed boundary method-related parameters
!-----------------------------------------------------------------------
use types, only: rprec
use param
implicit none
! ---
real(rprec), dimension(:, :, :), allocatable :: phi_uv, phi_w, phi_txz      ! signed-distance function phi(x) at uvp-node and w-node and txz-node
real(rprec), dimension(:, :, :, :), allocatable :: norm_phiuv, norm_phiw    ! normal direction of phi(x) for uvp-node and w-node
real(rprec), dimension(:, :), allocatable :: z_tpg      ! z_tpg - the nondimensional z-coordinate where phi=0
real(rprec), dimension(:, :, :), allocatable :: fx, fy, fz
real(rprec), dimension(3) :: delta_i        ! the increment of index for increasing one normalized unit length (zi) in each direction

real(rprec), dimension(:), allocatable, private :: arr_ibm

contains

    subroutine init_ibm()
    !-----------------------------------------------------------------------
    !   Initialize the immersed boundary layer
    !-----------------------------------------------------------------------
    implicit none
    ! ---
    integer(MPI_OFFSET_KIND) :: disp0, disp3d
    integer :: ind
    ! ---
    allocate(phi_uv(nxt, nynpy, 0:nz), phi_w(nxt, nynpy, 0:nz),                 &
             norm_phiuv(nxt, nynpy, 0:nz, 3), norm_phiw(nxt, nynpy, 0:nz, 3),   &
             phi_txz(nxt, nynpy, 0:nz), z_tpg(nxt, nynpy))

    allocate(fx(nxt, nynpy, nz), fy(nxt, nynpy, nz), fz(nxt, nynpy, nz))
    ! --- Calculate the thickness of immerse boundary and interpolation area
    phic = ratio_phic2dz * dz
    dphi = ratio_dphi * dz
    ! --- Readin topography
    disp3d = np * sizeof(phi_uv(1:nxt, 1:nynpy, 1:nz-1))

    disp0 = 0
    call read_file_mpi_3d(phi_uv(1:nxt, 1:nynpy, 1:nz-1), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(norm_phiuv(1:nxt, 1:nynpy, 1:nz-1, 1), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(norm_phiuv(1:nxt, 1:nynpy, 1:nz-1, 2), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(norm_phiuv(1:nxt, 1:nynpy, 1:nz-1, 3), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(phi_w(1:nxt, 1:nynpy, 1:nz-1), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(norm_phiw(1:nxt, 1:nynpy, 1:nz-1, 1), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(norm_phiw(1:nxt, 1:nynpy, 1:nz-1, 2), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(norm_phiw(1:nxt, 1:nynpy, 1:nz-1, 3), fn_ibm, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_xy(z_tpg(1:nxt, 1:nynpy), fn_ibm, disp0)

    ! ---
    call mpi_sendrecv(phi_uv(1, 1, nz - 1), nxt*nynpy, MPI_RPREC, up,   11, &
                      phi_uv(1, 1, 0),      nxt*nynpy, MPI_RPREC, down, 11, &
                      comm, status, ierr)
    call mpi_sendrecv(phi_uv(1, 1, 1),  nxt*nynpy, MPI_RPREC, down, 12, &
                      phi_uv(1, 1, nz), nxt*nynpy, MPI_RPREC, up,   12, &
                      comm, status, ierr)
    call mpi_sendrecv(phi_w(1, 1, nz - 1), nxt*nynpy, MPI_RPREC, up,   13,  &
                      phi_w(1, 1, 0),      nxt*nynpy, MPI_RPREC, down, 13,  &
                      comm, status, ierr)
    call mpi_sendrecv(phi_w(1, 1, 1),  nxt*nynpy, MPI_RPREC, down, 14, &
                      phi_w(1, 1, nz), nxt*nynpy, MPI_RPREC, up,   14, &
                      comm, status, ierr)

    do ind = 1, 3
        call mpi_sendrecv(norm_phiuv(1, 1, nz - 1, ind), nxt*nynpy, MPI_RPREC, up,   15,    &
                          norm_phiuv(1, 1, 0, ind),      nxt*nynpy, MPI_RPREC, down, 15,    &
                          comm, status, ierr)
        call mpi_sendrecv(norm_phiuv(1, 1, 1, ind),  nxt*nynpy, MPI_RPREC, down, 16,    &
                          norm_phiuv(1, 1, nz, ind), nxt*nynpy, MPI_RPREC, up,   16,    &
                          comm, status, ierr)
        call mpi_sendrecv(norm_phiw(1, 1, nz - 1, ind), nxt*nynpy, MPI_RPREC, up,   17, &
                          norm_phiw(1, 1, 0, ind),      nxt*nynpy, MPI_RPREC, down, 17, &
                          comm, status, ierr)
        call mpi_sendrecv(norm_phiw(1, 1, 1, ind),  nxt*nynpy, MPI_RPREC, down, 18, &
                          norm_phiw(1, 1, nz, ind), nxt*nynpy, MPI_RPREC, up,   18, &
                          comm, status, ierr)
    end do
  
    call mpi_bcast(z_tpg, nxt*nynpy, mpi_rprec, 0, comm_col, ierr)
    ! ---
    phi_txz = phi_w
    ! --- Set up the boundary condition
    if ( coordz == 0 ) then
        phi_uv(:, :, 0) = phi_uv(:, :, 1)
        phi_w(:, :, 0) = phi_w(:, :, 1)
        norm_phiuv(:, :, 0, :) = norm_phiuv(:, :, 1, :)
        norm_phiw(:, :, 0, :) = norm_phiw(:, :, 1, :)
        phi_txz(:, :, 1) = phi_uv(:, :, 1)
        phi_txz(:, :, 0) = phi_txz(:, :, 1)
    end if
    
    if ( coordz == npz-1 ) then
        phi_uv(:, :, nz) = phi_uv(:, :, nz-1)
        phi_w(:, :, nz) = phi_w(:, :, nz-1)
        norm_phiuv(:, :, nz, :) = norm_phiuv(:, :, nz-1, :)
        norm_phiw(:, :, nz, :) = norm_phiw(:, :, nz-1, :)
        phi_txz(:, :, nz) = phi_txz(:, :, nz-1)
    end if

    ! --- Initialize other variable
    delta_i = [1/dx, 1/dy, 1/dz]
    z0_ibm = z0_ibm / z_i
    allocate(arr_ibm(max(nxt2, nyt2)))      !!!!!!!!!!!
    ! --- 
    end subroutine init_ibm
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine smoothVel_cubic(fun, dir, phi)
    !-----------------------------------------------------------------------
    !   Smooth the velocity component at the location inside the immersed
    !   boundary using cubic function
    !----------------------------------------------------------------------- 
    implicit none
    ! ---
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(inout) :: fun
    character(*), intent(in) :: dir
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(in) :: phi
    integer :: iz
    ! ---
    do iz = 0, nz
        call smooth2D_cubic(fun(:, :, iz), dir, phi(:, :, iz), nxt, nynpy)
    end do
    ! ---
    end subroutine smoothVel_cubic
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine smoothTau_cubic(fun, dir, phi)
    !-----------------------------------------------------------------------
    !   Smooth the SGS stress component at the location inside the immersed
    !   boundary using cubic function
    !   The smooth area is where phi<=-phic/2
    !-----------------------------------------------------------------------
    implicit none
    ! ---    
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(inout) :: fun
    character(*), intent(in) :: dir
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(in) :: phi
    integer :: iz
    ! ---
    do iz=0, nz
        call smooth2D_cubic(fun(:, :, iz), dir, phi(:, :, iz), nxt, nynpy, -phic/2)
    end do
    ! ---
    end subroutine smoothTau_cubic
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine smooth2D_cubic(fun, dir, phi, nx, ny, phic)
    !-----------------------------------------------------------------------
    !   Smooth the a 2D array at the location inside the immersed
    !   boundary using cubic function
    !-----------------------------------------------------------------------
    implicit none
    ! ---
    integer, intent(in) :: nx, ny
    real(rprec), dimension(nx, ny), intent(inout) :: fun
    character(*), intent(in) :: dir
    real(rprec), dimension(nx, ny), intent(in) :: phi
    real(rprec), optional :: phic
    ! ---
    real(rprec) :: phicVal
    integer :: is, ie ! the indices enclosed the range for smoothing area
    integer :: ix, iy
    ! --- Set the critical phi value
    if (.not. present(phic)) then
        phicVal = 0.
    else
        phicVal = phic
    end if
        
    ! --- Use cubic function to smooth the quantities
    select case (trim(dir))
    case ('x')
        do iy = 1, ny
        do ix = 1, nx
            if (phi(ix, iy)>phicVal .and. phi(modulo(ix, nx)+1, iy)<=phicVal) then
                !- Find indices start and end bound
                is = modulo(ix-2, nx) + 1
                ie = modulo(ix+1, nx) + 1
                do while(phi(ie, iy)<=phicVal)
                    ie = modulo(ie, nx) + 1
                end do
                ie = modulo(ie, nx) + 1
        
                !- Interpolation
                if (is <= ie) then !LV5
                    call interp1d(fun(is:ie, iy), ie-is+1)
                else if (is > ie) then !LV5
                    arr_ibm(1:nx-is+1) = fun(is:nx, iy)
                    arr_ibm(nx-is+2:nx-is+1+ie) = fun(1:ie, iy)
                    call interp1d(arr_ibm(1:nx-is+1+ie), nx-is+1+ie)
                    fun(is:nx, iy) = arr_ibm(1:nx-is+1)
                    fun(1:ie, iy) = arr_ibm(nx-is+2:nx-is+1+ie)
                end if
            end if
        end do
        end do
    case ('y')
        do ix = 1, nx
        do iy = 1, ny
            if (phi(ix, iy)>phicVal .and. phi(ix, modulo(iy, ny)+1)<=phicVal) then
                !4- Find indices start and end bound
                is = modulo(iy-2, ny) + 1
                ie = modulo(iy+1, ny) + 1
                do while(phi(ix, ie)<=phicVal) 
                    ie = modulo(ie, ny) + 1
                end do
                ie = modulo(ie, ny) + 1
        
                !4- Interpolation
                if (is <= ie) then
                    call interp1d(fun(ix, is:ie), ie-is+1)
                else if (is > ie) then
                    arr_ibm(1:ny-is+1) = fun(ix, is:ny)
                    arr_ibm(ny-is+2:ny-is+1+ie) = fun(ix, 1:ie)
                    call interp1d(arr_ibm(1:ny-is+1+ie), ny-is+1+ie)
                    fun(ix, is:ny) = arr_ibm(1:ny-is+1)
                    fun(ix, 1:ie) = arr_ibm(ny-is+2:ny-is+1+ie)
                end if
            end if
        end do
        end do
    case default
        write(*,*) "Part smooth2D_cubic: direction selection error. STOP !!!"
        stop
    end select 
    ! ---
    end subroutine smooth2D_cubic
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine interp1d(arr, narr)
    !-----------------------------------------------------------------------
    !   Interpolate data with two points at each side
    !-----------------------------------------------------------------------    
    implicit none
    ! ---
    real(rprec), dimension(:), intent(inout) :: arr
    integer, intent(in) :: narr
    ! ---      
    real(rprec) :: coef(0:3)    ! coefficient of cubic interpolation
    real(rprec) :: xbi(4)       ! index of boundary points (two points at each side)
    real(rprec) :: di(4)        ! denominator used in each coefficients
    real(rprec) :: xi(3:narr-2) ! index of each point
    integer :: ind
    ! --- Check if the size of array (much larger than 4)
    if ( narr <= 4 ) then
        write(*,*) "The size of array for cubic interpolation must be larger than 4. STOP!"
        stop
    end if
    ! --- The relative coordinate location
    xi(3:narr-2) = [(dble(ind), ind=3, narr-2)]
    xbi(1) = dble(1)
    xbi(2) = dble(2)
    xbi(3) = dble(narr-1)
    xbi(4) = dble(narr)
          
    ! --- Calculate the coefficient of cubic interpolation
    di(1) = (xbi(2)-xbi(1))*(xbi(3)-xbi(1))*(xbi(4)-xbi(1))
    di(2) = (xbi(1)-xbi(2))*(xbi(3)-xbi(2))*(xbi(4)-xbi(2))
    di(3) = (xbi(1)-xbi(3))*(xbi(2)-xbi(3))*(xbi(4)-xbi(3))
    di(4) = (xbi(1)-xbi(4))*(xbi(2)-xbi(4))*(xbi(3)-xbi(4))
        
    coef(0) = xbi(2)*xbi(3)*xbi(4) / di(1) * arr(1) &
            + xbi(1)*xbi(3)*xbi(4) / di(2) * arr(2) &
            + xbi(1)*xbi(2)*xbi(4) / di(3) * arr(narr-1) &
            + xbi(1)*xbi(2)*xbi(3) / di(4) * arr(narr)
    coef(1) = -(xbi(2)*xbi(3)+xbi(3)*xbi(4)+xbi(4)*xbi(2)) / di(1) * arr(1) &
            -(xbi(1)*xbi(3)+xbi(3)*xbi(4)+xbi(4)*xbi(1)) / di(2) * arr(2) &
            -(xbi(1)*xbi(2)+xbi(2)*xbi(4)+xbi(4)*xbi(1)) / di(3) * arr(narr-1) &
            -(xbi(1)*xbi(2)+xbi(2)*xbi(3)+xbi(3)*xbi(1)) / di(4) * arr(narr)
    coef(2) = (xbi(2)+xbi(3)+xbi(4)) / di(1) * arr(1) &
            + (xbi(1)+xbi(3)+xbi(4)) / di(2) * arr(2) &
            + (xbi(1)+xbi(2)+xbi(4)) / di(3) * arr(narr-1) &
            + (xbi(1)+xbi(2)+xbi(3)) / di(4) * arr(narr)
    coef(3) = -arr(1)/di(1) - arr(2)/di(2) - arr(narr-1)/di(3) - arr(narr)/di(4)
        
    ! --- Interpolate the values
    arr(3:narr-2) = coef(3)*xi(3:narr-2)**3 + coef(2)*xi(3:narr-2)**2 &
                  + coef(1)*xi(3:narr-2)    + coef(0)
    ! ---
    end subroutine interp1d
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine restoreval_vel(fun, phi)
    !-----------------------------------------------------------------------
    !   Recover the vel at the location inside the immersed boundary
    !-----------------------------------------------------------------------
    implicit none
    ! ---
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(inout) :: fun
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(in) :: phi
    ! --- Recover the value
    where(phi<=0) fun=0._rprec
    ! ---    
    end subroutine restoreval_vel
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine restoreval_tau(fun, phi)
    !-----------------------------------------------------------------------
    !   Recover the SGS stress at the location inside the immersed boundary
    !-----------------------------------------------------------------------
    implicit none
    ! ---  
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(inout) :: fun
    real(rprec), dimension(nxt, nynpy, 0:nz), intent(in) :: phi
    ! --- Recover the value
    where(phi<=-phic/2) fun=0._rprec
    ! ---   
    end subroutine restoreval_tau
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine wallstress_ibm()
    !-----------------------------------------------------------------------
    !   Set the wall stress at the immsersed boundary.
    !   txx, txy, tyy, tzz are at uvp-node, txz, tyz are at w-node
    !   Only assign stress at the layer 0<phi<dphi. For phi<=0, use cubic
    !   function to smooth the stress
    !-----------------------------------------------------------------------
    use sim_param, only: u, v, w, txx, txy, txz, tyy, tyz, tzz
    implicit none
    ! --- Variable Declaration ###
    real(rprec), dimension(3) :: ixv    ! ixv - the index vector at the point of velocity interpolation. ixv=[ix, iy, iz]
    real(rprec), dimension(3) :: vel    ! vel - the interpolated velocity
    real(rprec), dimension(3) :: vel_t  ! vel_t - the tangertial component of interpolated velocity
    real(rprec), dimension(3) :: vel_n  ! vel_n - the normal component of interpolated velocity
    real(rprec) :: c000, c001, c010, c011, c100, c101, c110, c111
    real(rprec) :: tau                  ! tau - the wall stress with respect to rotated axis
    real(rprec), dimension(3) :: e2p    ! e2p - the unit vector of second rotated axis
    real(rprec) :: ai1, aj3, ai3, aj1   ! direction cosine between global axis and rotated axis
    ! ---      
    integer :: ixs, ixe, iys, iye, izs, ize
    integer :: ix, iy, iz
    ! --- uvec_x, uvec_y, uvec_z - unit vector in each axis
    real(rprec), dimension(3), parameter ::         &
        uvec_x = [1._rprec, 0._rprec, 0._rprec],    &
        uvec_y = [0._rprec, 1._rprec, 0._rprec],    &
        uvec_z = [0._rprec, 0._rprec, 1._rprec]

    ! --- Stress at uvp-node (txx, txy, tyy, tzz)
    !! Important: Notice that the txz and tyz of 1st level of whole domain are
    !! >>at uvp-node, and calculate here too
    do iz = 1, nz-1
    do iy = 1, nynpy
    do ix = 1, nxt
        if (phi_uv(ix, iy, iz)>-phic/2 .and. phi_uv(ix, iy, iz)<phic) then !LV2
        !-- Step 1: Interpolate the tangertial velocity at the distance dphi from the immersed boundary
            ixv = [ix, iy, iz]&
                + (dphi-phi_uv(ix, iy, iz)) * norm_phiuv(ix, iy, iz, :) * delta_i
            ixv(1) = modulo(ixv(1), dble(nxt))
            ixv(2) = modulo(ixv(2), dble(nyt))        !!!!!!!!!!
            !-- Get u, v at the dphi away from the immersed boundary
            ixs = modulo(floor(ixv(1)-1), nxt)+1
            ixe = modulo(ixs, nxt)+1
            iys = modulo(floor(ixv(2)-1), nyt)+1
            iye = modulo(iys, nyt)+1
            izs = floor(ixv(3))
            ize = izs+1
            c000 = (ixe-ixv(1)) * (iye-ixv(2)) * (ize-ixv(3))
            c001 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ize-ixv(3))
            c010 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3))
            c011 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3))
            c100 = (ixe-ixv(1)) * (iye-ixv(2)) * (ixv(3)-izs)
            c101 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ixv(3)-izs)
            c110 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)-izs)
            c111 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)-izs)
            izs = min(izs, nz)
            ize = min(ize, nz)
            vel(1:2) = c000*[u(ixs, iys, izs), v(ixs, iys, izs)]    &
                     + c001*[u(ixe, iys, izs), v(ixe, iys, izs)]    &
                     + c010*[u(ixs, iye, izs), v(ixs, iye, izs)]    &
                     + c011*[u(ixe, iye, izs), v(ixe, iye, izs)]    &
                     + c100*[u(ixs, iys, ize), v(ixs, iys, ize)]    &
                     + c101*[u(ixe, iys, ize), v(ixe, iys, ize)]    &
                     + c110*[u(ixs, iye, ize), v(ixs, iye, ize)]    &
                     + c111*[u(ixe, iye, ize), v(ixe, iye, ize)]
        
            ! --- Get w at the dphi away from the immersed boundary
            izs = floor(ixv(3)+0.5)
            ize = izs+1
            c000 = (ixe-ixv(1)) * (iye-ixv(2)) * (ize-ixv(3)-.5)
            c001 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ize-ixv(3)-.5)
            c010 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3)-.5)
            c011 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3)-.5)
            c100 = (ixe-ixv(1)) * (iye-ixv(2)) * (ixv(3)+.5-izs)
            c101 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ixv(3)+.5-izs)
            c110 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)+.5-izs)
            c111 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)+.5-izs)
            izs = min(izs, nz)
            ize = min(ize, nz)
            vel(3) = c000*w(ixs, iys, izs) + c001*w(ixe, iys, izs)&
                   + c010*w(ixs, iye, izs) + c011*w(ixe, iye, izs)&
                   + c100*w(ixs, iys, ize) + c101*w(ixe, iys, ize)&
                   + c110*w(ixs, iye, ize) + c111*w(ixe, iye, ize)
        
        ! --- Step 2: Calculate the stress at the current location using log-law
            vel_n = dot_product(vel, norm_phiuv(ix, iy, iz, :))&
                  * norm_phiuv(ix, iy, iz, :)
            vel_t = vel - vel_n
            e2p = cross_prod(norm_phiuv(ix, iy, iz, :), vel_t) / norm(vel_t)
            tau = -(vonk/log(dphi/z0_ibm))**2 * dot_product(vel_t, vel_t)
        
            ai1 = dot_product(uvec_x, vel_t/norm(vel_t))
            aj3 = dot_product(uvec_x, norm_phiuv(ix, iy, iz, :))
            ai3 = aj3
            aj1 = ai1
            txx(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
            aj3 = dot_product(uvec_y, norm_phiuv(ix, iy, iz, :))
            aj1 = dot_product(uvec_y, vel_t/norm(vel_t))
            txy(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
            ai1 = aj1
            ai3 = aj3
            tyy(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
            ai1 = dot_product(uvec_z, vel_t/norm(vel_t))
            aj3 = dot_product(uvec_z, norm_phiuv(ix, iy, iz, :))
            ai3 = aj3
            aj1 = ai1
            tzz(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
        
            ! --- Calculate the txz and tyz at the 1st level of the whole domain
            if (coordz == 0 .and. iz==1 ) then
                ai1 = dot_product(uvec_x, vel_t/norm(vel_t))
                aj3 = dot_product(uvec_z, norm_phiw(ix, iy, iz, :))
                ai3 = dot_product(uvec_x, norm_phiw(ix, iy, iz, :))
                aj1 = dot_product(uvec_z, vel_t/norm(vel_t))
                txz(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
                ai1 = dot_product(uvec_y, vel_t/norm(vel_t))
                ai3 = dot_product(uvec_y, norm_phiw(ix, iy, iz, :))
                tyz(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
            end if
        end if
    end do
    end do
    end do
        
    !! Stress at w-node (txz, tyz)
    !! Important: Notice that the txz and tyz of 1st level of whole domain are
    !! >>at uvp-node, and calculate at the previous loop
    do iz = 1, nz-1
        ! Skip the txz and tyz at the 1st level of whole domain
        if (iz <= 1 .and. coordz == 0) cycle
        
        do iy = 1, nynpy
        do ix = 1, nxt
            if (phi_txz(ix, iy, iz)>-phic/2 .and. phi_txz(ix, iy, iz)<phic) then
            ! Step 1: Interpolate the tangertial velocity at the distance dphi from the immersed boundary
                ixv = [ix, iy, iz]&
                    + (dphi-phi_txz(ix, iy, iz)) * norm_phiw(ix, iy, iz, :) * delta_i
                ixv(1) = modulo(ixv(1), dble(nxt))
                ixv(2) = modulo(ixv(2), dble(nyt))
                ! Get u, v at the dphi away from the immersed boundary
                ixs = modulo(floor(ixv(1)-1), nxt)+1
                ixe = modulo(ixs, nxt)+1
                iys = modulo(floor(ixv(2)-1), nyt)+1
                iye = modulo(iys, nyt)+1
                izs = floor(ixv(3)-0.5)
                ize = izs+1
                c000 = (ixe-ixv(1)) * (iye-ixv(2)) * (ize-ixv(3)+.5)
                c001 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ize-ixv(3)+.5)
                c010 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3)+.5)
                c011 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3)+.5)
                c100 = (ixe-ixv(1)) * (iye-ixv(2)) * (ixv(3)-.5-izs)
                c101 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ixv(3)-.5-izs)
                c110 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)-.5-izs)
                c111 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)-.5-izs)
                izs = min(izs, nz)
                ize = min(ize, nz)
                vel(1:2) = c000*[u(ixs, iys, izs), v(ixs, iys, izs)]    &
                         + c001*[u(ixe, iys, izs), v(ixe, iys, izs)]    &
                         + c010*[u(ixs, iye, izs), v(ixs, iye, izs)]    &
                         + c011*[u(ixe, iye, izs), v(ixe, iye, izs)]    &
                         + c100*[u(ixs, iys, ize), v(ixs, iys, ize)]    &
                         + c101*[u(ixe, iys, ize), v(ixe, iys, ize)]    &
                         + c110*[u(ixs, iye, ize), v(ixs, iye, ize)]    &
                         + c111*[u(ixe, iye, ize), v(ixe, iye, ize)]
        
                ! Get w at the dphi away from the immersed boundary
                izs = floor(ixv(3))
                ize = izs+1
                c000 = (ixe-ixv(1)) * (iye-ixv(2)) * (ize-ixv(3))
                c001 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ize-ixv(3))
                c010 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3))
                c011 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ize-ixv(3))
                c100 = (ixe-ixv(1)) * (iye-ixv(2)) * (ixv(3)-izs)
                c101 = (ixv(1)-modulo(ixs,nxt)) * (iye-ixv(2)) * (ixv(3)-izs)
                c110 = (ixe-ixv(1)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)-izs)
                c111 = (ixv(1)-modulo(ixs,nxt)) * (ixv(2)-modulo(iys,nyt)) * (ixv(3)-izs)
                izs = min(izs, nz)
                ize = min(ize, nz)
                vel(3) = c000*w(ixs, iys, izs) + c001*w(ixe, iys, izs)  &
                       + c010*w(ixs, iye, izs) + c011*w(ixe, iye, izs)  &
                       + c100*w(ixs, iys, ize) + c101*w(ixe, iys, ize)  & 
                       + c110*w(ixs, iye, ize) + c111*w(ixe, iye, ize)
        
                ! Step 2: Calculate the stress at the current location using log-law
                vel_n = dot_product(vel, norm_phiw(ix, iy, iz, :))* norm_phiw(ix, iy, iz, :)
                vel_t = vel - vel_n
                e2p = cross_prod(norm_phiw(ix, iy, iz, :), vel_t) / norm(vel_t)
                tau = -(vonk/log(dphi/z0_ibm))**2 * dot_product(vel_t, vel_t)
        
                ai1 = dot_product(uvec_x, vel_t/norm(vel_t))
                aj3 = dot_product(uvec_z, norm_phiw(ix, iy, iz, :))
                ai3 = dot_product(uvec_x, norm_phiw(ix, iy, iz, :))
                aj1 = dot_product(uvec_z, vel_t/norm(vel_t))
                txz(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
                ai1 = dot_product(uvec_y, vel_t/norm(vel_t))
                ai3 = dot_product(uvec_y, norm_phiw(ix, iy, iz, :))
                tyz(ix, iy, iz) = (ai1*aj3+ai3*aj1) * tau
            end if
        end do
        end do
    end do
    ! ---    
    end subroutine wallstress_ibm
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine bodyforce_ibm()
    !-----------------------------------------------------------------------
    !   Set a bodyforce inside the immsersed boundary in order to force the
    !   velicity field inside the immsersed boundary to zero
    !-----------------------------------------------------------------------
    use sim_param, only: u, v, w, dpdx, dpdy, dpdz
    implicit none
    ! ---   
    where (phi_uv(:, :, 1:nz-1) <= 0)
        fx(:, :, 1:nz-1) = tadv1*dpdx(:, :, 1:nz-1) - u(:, :, 1:nz-1)/dt
        fy(:, :, 1:nz-1) = tadv1*dpdy(:, :, 1:nz-1) - v(:, :, 1:nz-1)/dt
    end where
          
    where (phi_w(:, :, 1:nz-1)<=0) !LV1
        fz(:, :, 1:nz-1) = tadv1*dpdz(:, :, 1:nz-1) - w(:, :, 1:nz-1)/dt
    end where
    ! ---
    end subroutine bodyforce_ibm
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    subroutine csMask_ibm()
    !-----------------------------------------------------------------------
    !   Mask Cs inside the immersed boundary when using lagrangian dynamic model.
    !   Cs is at w-node except the lowest one.
    !-----------------------------------------------------------------------
    use sgsmodule, only: Cs_opt2 
    implicit none
    ! ---    
    if (coordz == 0) then 
        where (phi_uv(:, :, 1) <= 0)
            Cs_opt2(:, :, 1) = real(1.e-24, kind=rprec)
        end where
            
        where (phi_w(:, :, 2:nz) <= 0)
            Cs_opt2(:, :, 2:nz) = real(1.e-24, kind=rprec)
        end where
    else 
        where (phi_w(:, :, 1:nz) <= 0)
            Cs_opt2(:, :, 1:nz) = real(1.e-24, kind=rprec)
        end where 
    end if 
    ! ---   
    end subroutine csMask_ibm
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------    
    function cross_prod(vec1, vec2)
    !-----------------------------------------------------------------------
    !   Compute the cross product of two vectors
    !-----------------------------------------------------------------------
    implicit none
    ! ---
    real(rprec), dimension(3) :: cross_prod
    real(rprec), dimension(3) :: vec1, vec2
    ! ---
    cross_prod(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_prod(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_prod(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    ! ---
    end function cross_prod
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------      
    function norm(vec)
    !-----------------------------------------------------------------------
    !   norm of a vector
    !-----------------------------------------------------------------------
    implicit none
    ! ---
    real(rprec) :: norm
    real(rprec), dimension(3) :: vec
          
    norm = sqrt(dot_product(vec, vec))
    ! ---    
    end function norm
! ---
end module 
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------