!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_coef_press()
!-----------------------------------------------------------------------
!   initialize velocity field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use fft
implicit none
! ---
integer :: jx, jy, jz
! ---
allocate(coef_ap(nxhnpy, nyt, 0:nz), coef_bp(nxhnpy, nyt, 0:nz), coef_cp(nxhnpy, nyt, 0:nz))
coef_ap = 0._rprec
coef_bp = 0._rprec
coef_cp = 0._rprec
! ---
if ( coordz == 0 ) then
    coef_ap(:, :, 0) = BOGUS !--was 0._rprec
    coef_bp(:, :, 0) = -1._rprec
    coef_cp(:, :, 0) = 1._rprec

    izs = 1
    ize = nz - 1
else if ( coordz == npz - 1 ) then
    coef_ap(:, :, nz) = -1._rprec
    coef_bp(:, :, nz) = 1._rprec
    coef_cp(:, :, nz) = BOGUS !--was 0._rprec

    izs = 0
    ize = nz - 1
else 
    izs = 0
    ize = nz - 1
end if
! ---
do jz = izs, ize
do jy = 1, nyt
    if ( jy == nyt/2 + 1 ) cycle

    do jx = 1, nxhnpy-1
        if (coordy == 0 .and. jx*jy == 1) cycle
        
        coef_ap(jx, jy, jz) = 1._rprec/(dz**2)
        coef_bp(jx, jy, jz) = -(kx_2d_mpi(jx, jy)**2 + ky_2d_mpi(jx, jy)**2 + 2._rprec/(dz**2))
        coef_cp(jx, jy, jz) = 1._rprec/(dz**2)
    end do
end do
end do
! ---
end subroutine init_coef_press
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine press_stag_array(func, dfdx, dfdy)
! p_hat contains the physical space pressure on exit
!--provides p_hat, dfdx, dfdy 0:nz-1
!-------------------
! Boundary Layer version with 4th order derivs in vertical.
!  04 December 1995
!        Mods.
!        12/6: added upper and lower boundary conditions.
!        12/8: Corrected sign on imag. x and y deriv of RHS.
!        12/17 Forcing pressure equal zero at wall
!                        Removed forcing of <P>=0 for all z.
!                    Will need to change BC when hetero. surface stress.
!        12/17 Added ficticious node below wall (Note solver on Nz+1 x Nz+1
!          prev. version of this subroutine saved as p_bl4.old.12.17
!    12/18 Revised 2st deriv stencil at wall (avg of deriv at -1 and at 1)
!    12/21 Redid FDD to 2nd order accurate.
!    12/22 Now value of P(wall) diagnosed from prev P, to match gradient BC
!....1/13: major changes.
!                Broke out mean pressure for separate solution
!    1/20 back to tridag for natrix solution (same sol'n as LUDCMP...A Keeper!)
!....1/23 Staggered solution
!.........Using Nz+1 levels for computing P, but tossing out level below ground
!....4/1 Changed sign on Div T_iz at wall and lid (five places)
!-------------------
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
use types, only: rprec
use param
use sim_param, only: u, v, w, divtz
use fft
use bottombc, only: T_s
use intermediate, only: rH_x, rH_y, rH_z, H_x, H_y, H_z, &
                        rtopw, rbottomw, topw, bottomw
implicit none
! ---
real(kind=rprec), dimension(nxt, nynpy, 0:nz), intent(out) :: func
real(kind=rprec), dimension(nxt, nynpy, nz), intent(out) :: dfdx, dfdy
! ---
complex(kind=rprec), dimension(nxhnpy, nyt, 0:nz) :: p_hat
complex(kind=rprec), dimension(nxhnpy, nyt, nz) :: dpdx, dpdy

complex(kind=rprec), dimension(nxhnpy, nyt, 0:nz) :: RHS_col

real(kind=rprec), dimension(nxt, nynpy) :: surf_temp_pr_modi
real(kind=rprec) :: const, T_s_mn
integer :: jx, jy, jz, jz_min, jz_max

! ---
if ( coordz == 0 ) then
    p_hat(:, :, 0) = (0._rprec, 0._rprec)
else
    p_hat(:, :, 0) = BOGUS
end if

!#####################################################
!if ( coordz == npz - 1 ) then
!    p_hat(:, :, nz) = (0._rprec, 0._rprec)
!else 
!    p_hat(:, :, nz) = BOGUS
!end if
!#####################################################

! --- Get the right hand side ready (Loop over levels)
const = 1._rprec/(nxt*nyt)
!--experiment: was nz here (see below experiments)
! temp storage for sum of RHS terms.  normalized for fft
! sc: recall that the old timestep guys already contain the pressure term
! no forces
rH_x(:, :, 1:nz-1) = const / tadv1 * (u(:, :, 1:nz-1)/dt)
rH_y(:, :, 1:nz-1) = const / tadv1 * (v(:, :, 1:nz-1)/dt)
rH_z(:, :, 1:nz-1) = const / tadv1 * (w(:, :, 1:nz-1)/dt)

do jz = 1, nz - 1
    call fftn_mpi(rH_x(:, :, jz), H_x(:, :, jz))
    call fftn_mpi(rH_y(:, :, jz), H_y(:, :, jz))
    call fftn_mpi(rH_z(:, :, jz), H_z(:, :, jz))
end do

H_x(:, :, 0) = BOGUS
H_y(:, :, 0) = BOGUS
H_z(:, :, 0) = BOGUS

!--experiment this causes blow-up
H_x(:, :, nz) = BOGUS
H_y(:, :, nz) = BOGUS
  
if ( coordz == npz - 1 ) then
    H_z(:, :, nz) = (0._rprec, 0._rprec)
else
    H_z(:, :, nz) = BOGUS
end if

if ( coordz == 0 ) then
    if (theta_flag .and. (jt .gt. theta_init_time) .and. (lbc .eq. 0)) then
        ! Pressure bc modification at the surface for surf temp bc
        call calc_hor_avg3(T_s, T_s_mn)
        surf_temp_pr_modi = 0._rprec
        surf_temp_pr_modi = g*(z_i/(u_scale**2))*(T_s - T_s_mn)
        rbottomw(:, :) = const*(divtz(:, :, 1) - surf_temp_pr_modi(:, :))
    else
        rbottomw(:, :) = const * divtz(:, :, 1)
    end if
    call fftn_mpi(rbottomw(:, :), bottomw(:, :))
end if

if ( coordz == npz - 1 ) then
    rtopw(:, :) = const * divtz(:, :, nz)
    call fftn_mpi(rtopw(:, :), topw(:, :))
end if

! --- Loop over (Kx,Ky) to solve for Pressure amplitudes
!--switch order of inner/outer loops here
if ( coordz == 0 ) then
    RHS_col(:, :, 0) = -dz * bottomw(:, :)
else if ( coordz == npz - 1 ) then
    RHS_col(:, :, nz) = -dz * topw(:, :)
end if

!--couldx maybe combine some of these to less communication is needed
!--fill H_x, H_y, H_z at jz=0 (from nz-1)
!--cant just change lbz above, since u,v,w (jz=0) are not in sync yet
call mpi_sendrecv(H_x(1, 1, nz - 1), nxhnpy*nyt, MPI_CPREC, up,   1,  &
                  H_x(1, 1, 0),      nxhnpy*nyt, MPI_CPREC, down, 1,  &
                  comm, status, ierr)
call mpi_sendrecv(H_y(1, 1, nz - 1), nxhnpy*nyt, MPI_CPREC, up,   2,  &
                  H_y(1, 1, 0),      nxhnpy*nyt, MPI_CPREC, down, 2,  &
                  comm, status, ierr)
call mpi_sendrecv(H_z(1, 1, nz - 1), nxhnpy*nyt, MPI_CPREC, up,   3,  &
                  H_z(1, 1, 0),      nxhnpy*nyt, MPI_CPREC, down, 3,  &
                  comm, status, ierr)
!--fill H_x, H_y, H_z at jz=nz (from 1)
! call mpi_sendrecv(H_x(1, 1, 1),  nxhnpy*nyt, MPI_CPREC, down, 4,    &
!                   H_x(1, 1, nz), nxhnpy*nyt, MPI_CPREC, up,   4,    &
!                   comm, status, ierr)
! call mpi_sendrecv(H_y(1, 1, 1),  nxhnpy*nyt, MPI_CPREC, down, 5,    &
!                   H_y(1, 1, nz), nxhnpy*nyt, MPI_CPREC, up,   5,    &
!                   comm, status, ierr)
call mpi_sendrecv(H_z(1, 1, 1),  nxhnpy*nyt, MPI_CPREC, down, 6,    &
                  H_z(1, 1, nz), nxhnpy*nyt, MPI_CPREC, up,   6,    &
                  comm, status, ierr)

do jz = izs, ize
do jy = 1, nyt
    if (jy == nyt/2 + 1) cycle

    do jx = 1, nxhnpy - 1
        if (coordy == 0 .and. jx*jy == 1) cycle

        RHS_col(jx, jy, jz) = eye*(kx_2d_mpi(jx, jy)*H_x(jx, jy, jz)  + &
                                   ky_2d_mpi(jx, jy)*H_y(jx, jy, jz)) + &
                                  (H_z(jx, jy, jz+1) - H_z(jx, jy, jz))/dz
    end do
end do
end do

!--this skips zero wavenumber solution, nyquist freqs
call tridag_array_pipelined(0, coef_ap, coef_bp, coef_cp, RHS_col, p_hat)
! call mpi_barrier(comm, ierr)
!--zero-wavenumber solution
!--wait for p_hat(1, 1, 1) from "down"
if ( coordy == 0 ) then
    call mpi_recv(p_hat(1, 1, 1), 1, MPI_CPREC, down, 8, comm, status, ierr)

    if ( coordz == 0 ) then
        p_hat(1, 1, 0) = 0._rprec
        p_hat(1, 1, 1) = p_hat(1, 1, 0) - dz*bottomw(1, 1)
    end if

    !###################################################################
    !if ( coordz == npz - 1 ) then
    !    p_hat(1, 1, nz) = 0._rprec
    !    p_hat(1, 1, nz-1) = p_hat(1, 1, nz) + dz*topw(1, 1)
    !    jz_max = nz - 2
    !else
    !    jz_max = nz
    !end if
    !###################################################################

    do jz = 2, nz
        p_hat(1, 1, jz) = p_hat(1, 1, jz - 1) + H_z(1, 1, jz)*dz
    end do

    !--send p_hat(1, 1, nz) to "up"
    call mpi_send(p_hat(1, 1, nz), 1, MPI_CPREC, up, 8, comm, ierr)

    !--make sure 0 <-> nz-1 are syncronized
    !-- 1 <-> nz shouldx be in sync already
    call mpi_sendrecv(p_hat(1, 1, nz - 1), nxhnpy*nyt, MPI_CPREC, up,   2, &
                      p_hat(1, 1, 0),      nxhnpy*nyt, MPI_CPREC, down, 2, &
                      comm, status, ierr)
end if

!--zero the nyquist freqs
if ( coordy == npy - 1 ) then
    p_hat(nxhnpy, :, :) = 0._rprec
end if
p_hat(:, nyt/2 + 1, :) = 0._rprec

!===========================================================================
!...Now need to get p_hat(wave,level) to physical p(jx,jy,jz)
!.....Loop over height levels
call ifftn_mpi(p_hat(:, :, 0), func(:, :, 0))
do jz = 1, nz - 1
    do jy = 1, nyt
    do jx = 1, nxhnpy
        dpdx(jx, jy, jz) = eye*kx_2d_mpi(jx, jy)*p_hat(jx, jy, jz)
        dpdy(jx, jy, jz) = eye*ky_2d_mpi(jx, jy)*p_hat(jx, jy, jz)
    end do
    end do
    call ifftn_mpi(dpdx(:, :, jz),  dfdx(:, :, jz))
    call ifftn_mpi(dpdy(:, :, jz),  dfdy(:, :, jz))
    call ifftn_mpi(p_hat(:, :, jz), func(:, :, jz))
end do

!--nz level is not needed elsewhere (although its valid)
dfdx(:, :, nz) = BOGUS
dfdy(:, :, nz) = BOGUS
func(:, :, nz) = BOGUS
! ---
end subroutine press_stag_array
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine tridag_array_pipelined(tag, a, b, c, r, u)
!-----------------------------------------------------------------------
!--this assumes MPI
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016) 
!-----------------------------------------------------------------------
use types, only:rprec
use param
implicit none
! ---
integer, intent(in) :: tag !--base tag
real(kind=rprec), dimension(nxhnpy, nyt, 0:nz), intent(in):: a, b, c
complex(kind=rprec), dimension(nxhnpy, nyt, 0:nz), intent(in) :: r
complex(kind=rprec), dimension(nxhnpy, nyt, 0:nz), intent(out):: u

integer :: n, nchunks
character(64) :: fmt

integer :: chunksize
integer :: cstart, cend
integer :: jx, jy, j, j_min, j_max
integer :: tag0
integer :: q

real(kind=rprec), dimension(nxhnpy, nyt) :: bet
real(kind=rprec), dimension(nxhnpy, nyt, 0:nz) :: gam

! --- want to skip nyt/2+1 and 1, 1
n = nz
nchunks = nyt/2
chunksize = nyt/nchunks !--make sure nchunks divides nyt evenly

if (coordz == 0) then
    do jy = 1, nyt
    do jx = 1, nxhnpy - 1
        if (b(jx, jy, 0) == 0._rprec) then
            write (*, *) 'tridag_array: rewrite eqs, jx, jy= ', jx, jy
            stop
        end if
    end do
    end do

    bet = b(:, :, 0)
    u(:, :, 0) = r(:, :, 0)/bet     ! cyan should be u(:,:,1)

    j_min = 0 !--this is only for backward pass
else
    j_min = 1 !--this is only for backward pass
end if

if (coordz == npz - 1) then
    j_max = n
else
    j_max = n - 1
end if


do q = 1, nchunks

    cstart = 1 + (q - 1)*chunksize
    cend = cstart + chunksize - 1

    tag0 = tag + 10*(q - 1)

    if (coordz /= 0) then
        !--wait for c(:,:,1), bet(:,:), u(:,:,1) from "down"
        !--may want irecv here with a wait at the end
        call mpi_recv(c(1, cstart, 0), nxhnpy*chunksize, MPI_RPREC, down, tag0 + 1, &
                      comm, status, ierr)
        call mpi_recv(bet(1, cstart), nxhnpy*chunksize, MPI_RPREC, down, tag0 + 2, &
                      comm, status, ierr)
        call mpi_recv(u(1, cstart, 0), nxhnpy*chunksize, MPI_CPREC, down, tag0 + 3, &
                      comm, status, ierr)
    end if

    
    do j = 1, j_max
    do jy = cstart, cend
        if (jy == nyt/2 + 1) cycle

        do jx = 1, nxhnpy - 1
            if (coordy == 0 .and. jx*jy == 1) cycle

            gam(jx, jy, j) = c(jx, jy, j - 1)/bet(jx, jy)
            bet(jx, jy) = b(jx, jy, j) - a(jx, jy, j)*gam(jx, jy, j)

            if (bet(jx, jy) == 0._rprec) then
                write (*, *) 'tridag_array failed at jx,jy,j=', jx, jy, j
                write (*, *) 'a,b,c,gam,bet=', a(jx, jy, j), b(jx, jy, j), &
                            c(jx, jy, j-1), gam(jx, jy, j), bet(jx, jy)
                stop
            end if

            u(jx, jy, j) = (r(jx, jy, j) - a(jx, jy, j)*u(jx, jy, j - 1))/ &
                            bet(jx, jy)
        end do
    end do
    end do

    if (coordz /= npz - 1) then
        !--send c(n-1), bet, u(n-1) to "up"
        !--may not want blocking sends here
        call mpi_send(c(1, cstart, n - 1), nxhnpy*chunksize, MPI_RPREC, up, tag0 + 1, &
                      comm, ierr)
        call mpi_send(bet(1, cstart), nxhnpy*chunksize, MPI_RPREC, up, tag0 + 2, &
                      comm, ierr)
        call mpi_send(u(1, cstart, n - 1), nxhnpy*chunksize, MPI_CPREC, up, tag0 + 3, &
                      comm, ierr)
    end if
end do

  
do q = 1, nchunks
    cstart = 1 + (q - 1)*chunksize
    cend = cstart + chunksize - 1

    tag0 = tag + 10*(q - 1)

    if (coordz /= npz - 1) then
      !--wait for u(n), gam(n) from "up"
        call mpi_recv(u(1, cstart, n), nxhnpy*chunksize, MPI_CPREC, up, tag0 + 4, &
                      comm, status, ierr)
        call mpi_recv(gam(1, cstart, n), nxhnpy*chunksize, MPI_RPREC, up, tag0 + 5, &
                      comm, status, ierr)
    end if

!if ( coordz == npz-1 ) then
!  j_max = n-1
!else
!  j_max = n
!end if

    do j = n - 1, j_min, -1
        !--intend on removing cycle statements/repl with something faster
        do jy = cstart, cend
            if (jy == nyt/2 + 1) cycle

            do jx = 1, nxhnpy - 1
                if (coordy == 0 .and. jx*jy == 1) cycle
                u(jx, jy, j) = u(jx, jy, j) - gam(jx, jy, j + 1)*u(jx, jy, j + 1)
            end do
        end do
    end do

    !--send u(2), gam(2) to "down"
    call mpi_send(u(1, cstart, 1), nxhnpy*chunksize, MPI_CPREC, down, tag0 + 4, &
                  comm, ierr)
    call mpi_send(gam(1, cstart, 1), nxhnpy*chunksize, MPI_RPREC, down, tag0 + 5, &
                  comm, ierr)
end do
! ---
end subroutine tridag_array_pipelined
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------