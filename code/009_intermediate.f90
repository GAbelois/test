!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
module intermediate
!-----------------------------------------------------------------------
!   save the intermediate variables
!-----------------------------------------------------------------------
use types, only: rprec
!use param, only: nxt, nyt, nz

implicit none
save
public
!### Output Setup ###
 !## scalar_modules2.f90, io.f90
 ! The variable average_dim_select_flag generates the following values based
 ! >>on the value of average_dim_num in param.f90 :-
 ! >>a) average_dim_num = 2 :
 ! >>average_dim_select_flag = 0 ==> average the 3d array over x and y and
 ! >>output the z profile
 ! >>b) average_dim_num = 1 :
 ! >>average_dim_select_flag = 1 ==> average the 3d array over y and
 ! >>output the (x,z) profile
integer :: average_dim_select_flag
integer :: dim1_size
integer :: dim2_size
integer :: dim1_global
integer :: dim2_global
! --- io.f90 and scalar_module2.f90, subroutine **slice
real(rprec), dimension(:), allocatable :: ax_1d
real(rprec), dimension(:), allocatable :: avg_out_1d
real(rprec), dimension(:, :), allocatable :: avg_out
! --- io.f90, subroutine avgslice
real(rprec), dimension(:, :), allocatable :: ap, au, av, aw, p2, u2, v2, w2,    &
                                             auw, avw, acs
real(rprec), dimension(:, :), allocatable :: adudz, advdz, aCs_Ssim, abeta_sgs, &
                                             abetaclip_sgs
real(rprec), dimension(:, :), allocatable :: atxx, atxz, atyy, atyz, atzz
real(rprec), dimension(:, :), allocatable :: u3, v3, w3
real(rprec), dimension(:, :), allocatable :: adudt, advdt, adwdt
! --- scalar_module2.f90, subroutine scalar_slice and pollen_slice
real(rprec), dimension(:, :), allocatable :: atheta, t2, q2, asgs_t3, awt
real(rprec), dimension(:, :), allocatable :: adTdz, anu_t, t3, var_t
real(rprec), dimension(:, :), allocatable :: aPCon, PCon2, asgs_PCon3, awPCon
real(rprec), dimension(:, :), allocatable :: adPCondz, var_PCon, aKc_t, aCs2Sc
real(rprec), dimension(:, :), allocatable :: asgs_PCon1, auPCon
! --- io.f90, subroutine avg_stats
real(rprec), dimension(:, :), allocatable :: ubar_avg, vbar_avg, thetabar_avg,  &
                                             Cs2bar_avg, Nutbar_avg, ubar_tot,  &
                                             vbar_tot, thetabar_tot, Cs2bar_tot, Nutbar_tot
real(rprec), dimension(:, :), allocatable :: upr2bar_avg, stressbar_avg, &
                                             upr2bar_tot, stressbar_tot
real(rprec), dimension(:, :, :), allocatable :: Eozbar_avg, Eozbar_tot
! --- scalar_module2.f90, subroutine budget_TKE_scalar
real(rprec), dimension(:, :), allocatable :: awu2, awv2, awT2
real(rprec), dimension(:, :), allocatable :: auv, awp
real(rprec), dimension(:, :), allocatable :: autau13, avtau23, adissip

!### Flow Parameters ###
! --- convec.f90
real(rprec), dimension(:, :, :), allocatable :: cc_big
real(rprec), dimension(:, :, :), allocatable :: u_big, v_big, w_big
real(rprec), dimension(:, :, :), allocatable :: cross_x_big, cross_y_big, cross_z_big

real(rprec), dimension(:, :, :), allocatable :: cross_x, cross_y, cross_z
real(rprec), dimension(:, :, :), allocatable :: vort1, vort2, vort3
real(rprec), dimension(:, :, :), allocatable :: vort1_big, vort2_big, vort3_big
! --- scalar_module.f90
real(rprec), dimension(:, :, :), allocatable :: u_cell, v_cell
real(rprec), dimension(:, :, :, :), allocatable :: pu0_cell, pv0_cell
real(rprec), dimension(:, :, :, :), allocatable :: pu_cell, pv_cell
real(rprec), dimension(:, :), allocatable :: decay_cell
! --- press_stag_array.f90
real(rprec), dimension(:, :, :), allocatable, target :: rH_x, rH_y, rH_z
complex(rprec), dimension(:, :, :), allocatable :: H_x, H_y, H_z
real(rprec), dimension(:, :), allocatable :: rtopw, rbottomw
complex(rprec), dimension(:, :), allocatable :: topw, bottomw
 
!### Temperature Parameters ###
 !## dynamic_sc.f90, subroutine interpolag_scalar and
 !## >>interpolag_scalar_Sdep
real(rprec), dimension(:, :, :), allocatable :: FF_KX, FF_XX, FF_KX2, FF_XX2

! ---
end module intermediate
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
