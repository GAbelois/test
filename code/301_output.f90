
subroutine output_field(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, p, theta, PCon
use stokes_drift, only: ust, vst
use sgsmodule, only: dissip
implicit none
! ---
integer, intent(in) :: num

real(rprec), dimension(nxt, nynpy, 0:nz) :: p_org
character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
do jz = 1, nz-1
    p_org(:, :,jz) = p(:, :,jz) - 0.5*(u(:, :, jz)*u(:, :, jz) +   &
                                       v(:, :, jz)*v(:, :, jz) +   &
                                       w(:, :, jz)*w(:, :, jz))    &
                                -     (u(:, :, jz)*ust(jz)+v(:, :, jz)*vst(jz))   &
                                  - 0.5*(ust(jz)**2 + vst(jz)**2)
end do
! ---
disp3d = np * sizeof(real(u(1:nxt, 1:nynpy, 1:nz-1)))

write(fname, '(a,i8.8,a)') '../output/out/field_', num, '.out'
disp0 = 0
call write_file_mpi_3d_sp(real(u(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
disp0 = disp0+disp3d
call write_file_mpi_3d_sp(real(v(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
disp0 = disp0+disp3d
call write_file_mpi_3d_sp(real(w(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
disp0 = disp0+disp3d

if (flag_opt_pre) then
    call write_file_mpi_3d_sp(real(p_org(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
    disp0 = disp0+disp3d
end if

if (flag_opt_eps) then
    call write_file_mpi_3d_sp(real(dissip(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
    disp0 = disp0+disp3d
end if

if ( pcon_flag ) then
    do ipcon = 1, npcon 
        call write_file_mpi_3d_sp(real(PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon)), fname, disp0)
        disp0 = disp0+disp3d
    end do
end if

if ( theta_flag ) then
    call write_file_mpi_3d_sp(real(theta(1:nxt, 1:nynpy, 1:nz-1)), fname, disp0)
end if
! ---
end subroutine output_field
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
subroutine output_fringe(num)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
use types, only: rprec
use param
use sim_param, only: u, v, w, theta, PCon   ! , p, RHSx, RHSy, RHSz
! use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX, F_XX,F_KX2,F_XX2
! use scalars_module, only : RHS_T
implicit none
! ---
integer, intent(in) :: num

character(80) :: fname
integer(MPI_OFFSET_KIND) :: disp0, disp3d
integer :: jz, ipcon
! ---
disp3d = npy * npz * sizeof(u(jx_fringe+1:nxt, 1:nynpy, 1:nz-1))
! --- output velocity field at the outlet
write(fname, '(a,i8.8,a)') '../output/inflow/vel_tt', num, '.out'
disp0 = 0
call write_file_mpi_fringe(u(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0 + disp3d
call write_file_mpi_fringe(v(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
disp0 = disp0 + disp3d
call write_file_mpi_fringe(w(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(RHSx(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(RHSy(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(RHSz(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(Cs_opt2(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_LM(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_MM(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_QN(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! disp0 = disp0 + disp2d
! call write_file_mpi_yz(F_NN(nxt, 1:nynpy, 1:nz-1), fname, disp0)
! --- output the temperature field
if (theta_flag) then
    write(fname, '(a,i8.8,a)') '../output/inflow/temp_tt', num, '.out'
    disp0 = 0
    call write_file_mpi_fringe(theta(jx_fringe+1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    ! disp0 = disp0 + disp2d
    ! call write_file_mpi_yz(RHS_T(nxt, 1:nynpy, 1:nz-1), fname, disp0)
end if

if (PCon_flag) then		
    write(fname, '(a,i8.8,a)') '../output/inflow/pcon_tt', num, '.out'		
    disp0 = 0		
    do ipcon = 1, npcon		
        call write_file_mpi_fringe(PCon(jx_fringe+1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)		
        disp0 = disp0 + disp3d		
    end do		
end if
! ---
end subroutine output_fringe
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!
!subroutine output_loop
!!--------------------------------------------------------------------!
!!BC revised by Bicheng Chen to add output sample
!!-BC add variable flagVSpl, ttVSpl, flagCSpl, ttCSpl
!!-BC add new subroutine checkVSpl and checkCSpl                                               
!!--------------------------------------------------------------------! 
!use param,only:path, output, c_count, theta_flag, theta_init_time, jt, jt_total,  &
!               jan_diurnal_run, flagVSpl, ttVSpl, flagCSpl,             &
!               ttCSpl, flag_srfV, tt_srfV, PCon_FLAG, PCon_init,        &
!               use_avgslice, base_time
!use scalars_module2,only:scalar_slice, pollen_slice, budget_TKE_scalar
!use io, only:calc_mean, avgslice, checkpoint, checkVSpl, checkCSpl,     &
!             check_srfV, post_spec, io_spec, spec_write_start,          &
!             spec_write_end, spec_write_freqz
!implicit none
!! ---
!jt_total = jt_total + 1
!!call calc_mean()
!
!!cyan if (output) then
!!    if (mod(jt_total, base_time)==0) then
!!    !-- move all stuff into checkpoint (Bicheng Chen 06/26/2015)
!!        call checkpoint()
!!    end if
!!end if
!! --- BC added by Bicheng Chen for output the sample data
!if (flagVSpl) then
!    if (mod(jt_total, ttVSpl)==0) then
!        call checkVSpl()
!    end if
!end if
!
!if (flagCSpl) then
!    if (mod(jt_total, ttCSpl)==0) then
!        call checkCSpl()
!    end if
!end if
!!BC END
!! ---
!if (flag_srfV) then
!    if (mod(jt_total, tt_srfV)==0) then
!        call check_srfV()
!    end if
!end if
!
!! ---
!!cyanif ((io_spec) .and. (jt_total .gt. spec_write_start .and. jt_total .le. spec_write_end)) then
!!!  if ((io_spec) .and. mod(jt_total,spec_write_freqz)==0) then
!!    if (mod(jt_total,spec_write_freqz)==0) call post_spec(jt_total)
!!end if
!
!!cyan if (time_spec.gt.0) call timeseries_spec
!! ---
!end subroutine output_loop
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------