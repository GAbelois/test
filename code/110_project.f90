!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine project()
!-----------------------------------------------------------------------
!--provides u, v, w at 1:nz
!-----------------------------------------------------------------------
use param
use sim_param
! use scalars_module, only:inflow_pollen
use ibm
implicit none
! ---
integer :: jx, jy, jz, jz_min
real(rprec) :: RHS
! ---
if ( coordz == 0 ) then
    jz_min = 2
else
    jz_min = 1
end if

u(:, :, 1:nz-1) = u(:, :, 1:nz-1) - dt * tadv1 * dpdx(:, :, 1:nz-1)
v(:, :, 1:nz-1) = v(:, :, 1:nz-1) - dt * tadv1 * dpdy(:, :, 1:nz-1)
w(:, :, jz_min:nz-1) = w(:, :, jz_min:nz-1) - dt * tadv1 * dpdz(:, :, jz_min:nz-1)

if (flag_ibm) then
    u(:, :, 1:nz-1) = u(:, :, 1:nz-1) + dt*fx(:, :, 1:nz-1)
    v(:, :, 1:nz-1) = v(:, :, 1:nz-1) + dt*fy(:, :, 1:nz-1)
    w(:, :, jz_min:nz-1) = w(:, :, jz_min:nz-1) + dt*fz(:, :, jz_min:nz-1)
end if
! ---
if ( inflow ) then
    if ( read_inflow_file ) then
        call readin_inflow(nums)
        call update_fringe()
    end if
end if

!--enfore bc at top
if ( coordz == npz - 1 ) then
    if ( ubc_mom == 'stress free' ) then    ! cyan 2019-01-15
        u(:, :, nz-1) = u(:, :, nz-2)
        v(:, :, nz-1) = v(:, :, nz-2)    
    end if
    w(:, :, nz) = 0._rprec
end if

if ( coordz == 0 ) then
    w(:, :, 1) = 0._rprec
end if
! ---
end subroutine project
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine readin_inflow(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, p, RHSx, RHSy, RHSz, theta, PCon
use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX, F_XX,F_KX2,F_XX2
use scalars_module, only : RHS_T
use canopy, only:nutrient_inflow
implicit none
! ---
integer, intent(in) :: num

character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
disp3d = npy * npz * sizeof(u(jx_fringe+1:nxt, 1:nynpy, 1:nz-1))
! --- output velocity field at the outlet
write(fname, '(a,i8.8,a)') '../../../inflow/case4/output/inflow/vel_tt', num, '.out'
disp0 = 0
call read_file_mpi_fringe(u(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0 + disp3d
call read_file_mpi_fringe(v(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0 + disp3d
call read_file_mpi_fringe(w(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(RHSx(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(RHSy(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(RHSz(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(Cs_opt2(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(F_LM(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(F_MM(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(F_QN(1, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call read_file_mpi_yz(F_NN(1, 1:nynpy, 1:nz-1), fname, disp0)
! --- output the temperature field
if (theta_flag) then
    write(fname, '(a,i8.8,a)') '../../../inflow/case4/output/inflow/temp_tt', num, '.out'
    disp0 = 0
    call read_file_mpi_fringe(theta(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    ! disp0 = disp0 + disp2d
    ! call read_file_mpi_yz(RHS_T(1, 1:nynpy, 1:nz-1), fname, disp0)
end if
 
if (PCon_flag) then
    write(fname, '(a,i8.8,a)') '../../../inflow/case4/output/inflow/pcon_tt', num, '.out'
    disp0 = 0
    do ipcon = 1, npcon		
        call read_file_mpi_fringe(PCon(jx_fringe+1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)		
        disp0 = disp0 + disp3d		
    end do
end if 
! ---
end subroutine readin_inflow
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine update_fringe()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, theta, PCon
use scalars_module, only : RHS_T
implicit none
! ---
integer :: i, jx, ipcon

do i = 1, jx_relax
    jx = jx_fringe + i
    u(jx, :, :) = (1._rprec - wgt(i)) * u(jx_fringe, :, :) + wgt(i) * u(jx, :, :)
    v(jx, :, :) = (1._rprec - wgt(i)) * v(jx_fringe, :, :) + wgt(i) * v(jx, :, :)
    w(jx, :, :) = (1._rprec - wgt(i)) * w(jx_fringe, :, :) + wgt(i) * w(jx, :, :)
    theta(jx, :, :) = (1._rprec - wgt(i)) * theta(jx_fringe, :, :) + wgt(i) * theta(jx, :, :)
end do

if (PCon_flag) then		
    do ipcon = 1, npcon		
    do i = 1, jx_relax		
        jx = jx_fringe + i		
        PCon(jx, :, :, ipcon) = (1._rprec - wgt(i)) * PCon(jx_fringe, :, :, ipcon) +    &		
                                            wgt(i)  * PCon(jx, :, :, ipcon)		
    end do		
    end do		
end if
! ---
end subroutine update_fringe
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!

