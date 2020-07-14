!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine init()
!--------------------------------------------------------------------!
!   initialize everything to make main program clean                                           
!--------------------------------------------------------------------!    
use param, only: theta_flag, pcon_flag, model, ubc, sgs, comm, ierr, flag_ibm
use sim_param, only: u, v, w
use io, only: screen_display
use intermediate
use bottombc, only: patches
use fft
use test_filtermodule
use topbc, only: sponge, setsponge
use ibm, only: init_ibm
implicit none
! ---
call init_namelist()            ! read all pararmeters from namelist 
! ---
call init_parameter()           ! domain variables initialization
! ---
call init_nondimensional()      ! normalize the variable
! ---
call allocate_output_variable() ! output variables initialization
call allocate_flow_variable()   ! flow variables initialization
! ---
call init_vel_field()
call init_temperature_field()
call init_concentration_field()
call data_exchange()            ! exchange data between neighboring zones
! ---
call patches()                  ! Initialize surface physical conditions
! --- formulate the fft plans
call init_fft_plan()            ! create fft plans and initialize the kx, ky arrays
! --- 
!call openfiles()                ! open output files
! --- initialize test filter
if ( sgs ) then
    call init_test_filter(2._rprec * filter_size, G_test)

    if (model == 3 .or. model == 5) then  !--scale dependent dynamic
        call init_test_filter(4._rprec * filter_size, G_test_test)
    end if
end if
! --- define the upper boundary condition (sponge damping layer)
if (ubc == 1) then
    call setsponge()
else
    sponge = 0._rprec
end if
! ---
if ( flag_ibm ) then
    call init_ibm()
end if
call init_lad()
call init_coef_press()
! --- 
call filter_data(u)
call filter_data(v)
call filter_data(w)
! ---
! call screen_display()           ! echo the simulation parameters
! ---
end subroutine init
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine init_namelist()
!--------------------------------------------------------------------!
!     Read all the pararmeters in namelist                                        
!--------------------------------------------------------------------!
use param
implicit none
! ---
open(fid_param, file=fn_param, form='formatted', status='old')
read(fid_param, nml=output_control)
read (fid_param, nml=domain_param)
read (fid_param, nml=time_param)
read (fid_param, nml=flow_param)
read (fid_param, nml=canopy_param)
read (fid_param, nml=ibm_param)

read (fid_param, nml=temperature_param)
read (fid_param, nml=ocean_param)  
read (fid_param, nml=con_param)

if (PCon_flag) then
    info_con%n_con => n_con
    if (settling) then
        allocate (vel_settling(n_con))
        allocate (ratio_dens(n_con))
        allocate (vol_spec(n_con))
        read(fid_param, nml=particle_param)
        info_con%vel_settling => vel_settling
        info_con%ratio_dens => ratio_dens
        info_con%vol_spec => vol_spec
    end if 
end if
close (fid_param)

if ( flag_canopy ) then
    open(fid_leaf, file=fn_leaf, form='formatted', status='old')
    read(fid_leaf, nml=leaf_param)
    close(fid_leaf)
end if
! ---
end subroutine init_namelist
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_parameter()
!-----------------------------------------------------------------------
!   domain variables initialization (Bicheng Chen 07/02/2015)  
!-----------------------------------------------------------------------
use param
use scalars_module, only:hetero_count_out
use io, only:gsize, lsize, start
implicit none
! ---
integer :: ind
! --- 
nz   = nzt/npz + 1

nxt2  = 3*nxt/2
nyt2  = 3*nyt/2
lhx  = nxt/2 + 1
lhy  = nyt/2 + 1
lhx2 = nxt2/2 + 1
lhy2 = nyt2/2 + 1

nxnpy  = nxt / npy
nynpy  = nyt / npy
nx2npy = nxt2/ npy
ny2npy = nyt2/ npy
nxhnpy = nxnpy/2 + 1
nxh2npy= nx2npy/2 + 1

dx   = lx_tot / nxt
dy   = ly_tot / nyt
dz   = lz_tot / (nzt - 1)

ly   = ly_tot / npy 
lz   = lz_tot / npz

! --- initialize other time variables
hetero_count_out = p_count
! ---
! if (flag_restart) then
!     info_time%flag_restart = flag_restart
!     info_time%nt_restart = nt_restart
! end if
! ---
gsize(1) = nxt
gsize(2) = nyt
gsize(3) = nzt-1
lsize(1) = nxt
lsize(2) = nynpy
lsize(3) = nz-1
start(1) = 0
start(2) = nynpy * coordy
start(3) = (nz-1) * coordz
! ---
if ( mod(nxt2, npy) /= 0 .and. mod(nyt2, npy) /= 0 ) then
    write(*, *) 'Grid numbers and number of processors in y are not compatible'
    call mpi_finalize(ierr)
    stop
end if
! ---
end subroutine init_parameter
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_nondimensional()
!-----------------------------------------------------------------------
!   normalize the variables
!-----------------------------------------------------------------------
use types, only:rprec
use param
use stokes_drift
use intermediate
implicit none
! --- 
ly = ly / z_i
lz = lz / z_i

lx_tot = lx_tot/z_i
ly_tot = ly_tot/z_i
lz_tot = lz_tot/z_i

dx = dx/z_i
dy = dy/z_i
dz = dz/z_i

call init_meshgrid()
call fringe_prep()
! --- 
dt = dt * u_scale / z_i
u_star = u_star / u_scale
! --- canopy normalization
if ( flag_canopy ) then
    hc  = hc  / z_i
    hgap = hgap / z_i
    lad = lad * z_i
end if
! --- Coriolis effect
if (coriolis_forcing) then
    coriol = freq_coriolis*z_i/u_scale
    ug = ug_dim/u_scale
    vg = vg_dim/u_scale
end if
! --- Mean pressure gradient force
if (use_mean_p_force) then
    if (mean_p_force == float_missing) then
        mean_p_force = 1._rprec / lz_tot
    else
        mean_p_force = mean_p_force / u_scale**2 * z_i
    end if
else 
    mean_p_force = 0._rprec
end if

! ---
if (ocean_flag .and. theta_flag) then
    alpha_w = alpha_w * T_scale     ! normalize thermal expansion rate
end if

! --- Ocean flag
call init_stokes_drift

if (ocean_flag) then
    coriol = freq_coriolis*z_i/u_scale
    ug = ug_dim/u_scale
    vg = vg_dim/u_scale
    !- if use dynamical stress, read friction velocity and its angle
    if (flag_dynStress) then
        inquire (iolength=len_ustar) ustar_dyn, agl_stress
        open (unit=fid_ustar, file=fn_ustar, access='direct', recl=len_ustar)
        read (unit=fid_ustar, rec=1) ustar_dyn, agl_stress
        close (fid_ustar)
        ustar_dyn = ustar_dyn / u_scale
    end if 
    rad_stress = agl_stress/180._rprec*pi
    
    ! >>Cross flow over all domain
    if (flag_crossVel) then
        u_cross = udim_cross / u_scale
        v_cross = vdim_cross / u_scale
    end if
end if
! ---
end subroutine init_nondimensional
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_meshgrid()
!-----------------------------------------------------------------------
!    
!-----------------------------------------------------------------------
use param
implicit none
! ---
integer :: i, jx, jy, jz
! ---
allocate(x(nxt), y(nynpy), zuvp(0:nz), zw(0:nz))
allocate(fz_uvp(0:nz), fz_w(0:nz))
! --- 
select case (meshtype)
case('uniform')
    do jz = 0, nz
        zw(jz) = (coordz*(nz-1) + jz - 1) * dz
        ! zw(jz) = (coordz*(nz-1) + jz - 1) * (lz_tot/(nzt - 1))
    end do
case('cosine')
    do jz = 0, nz
        zw(jz) = (1._rprec - cos((coordz*(nz-1) + jz - 1)/(nzt-1) * pi/2._rprec)) * lz_tot
    end do
case('pbl')     ! uniform in the inner sublayer, and stretched in the transition layer, and uniform in the outer layer
!   ratio_inner = ; ratio_transition = ; ratio_outer = 1._rprec - ratio_inner - ratio_transition
!   lz_inner = lz_tot*ratio_inner; lz_transition = lz_tot*ratio_transition; lz_outer = lz_tot*ratio_outer
!   nz_outer = nzt - nz_inner - nz_transition
!   dz = lz_inner/nz_inner; dz2 = lz_outer/nz_outer;
!   ratio_grid = (dz2/dz)**(1._rprec/dble(nz_transition))
!   do jz = 0, nz_inner
!       zt_w(jz) = (jz-1) * dz
!   end do
!   do jz = nz_inner, nz_inner+nz_transition
!       zt_w(jz) = zt_w(nz_inner) + ratio_grid**(jz-nz_inner) * dz
!   end do
!   do jz = nz_inner+nz_transition, nzt
!       zt_w(jz) = zt_w(nz_inner+nz_transition) + (jz-nz_inner-nz_transition)*dz2
!   end do
!   call scatter_data_z(zt_w, zw)
case default
    write(*,*) 'invalid mesh type'
    stop
end select
zuvp(0:nz) = 0.5_rprec * (zw(0:nz-1) + zw(1:nz))

do jz = 0, nz-1
    fz_uvp(jz) = (zuvp(jz+1) - zw(jz+1)) / (zuvp(jz+1) - zuvp(jz))
    fz_w(jz)   = (zw(jz+1) - zuvp(jz)) / (zw(jz+1) - zw(jz))
end do

do jx = 1, nxt 
    x(jx) = (jx - 1) * dx
end do

do jy = 1, nynpy
    y(jy) = (coordy*nynpy + jy - 1) * dy
end do
! ---
end subroutine init_meshgrid
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine fringe_prep()
!-----------------------------------------------------------------------
!    Fringe method preparation
!-----------------------------------------------------------------------
use param
implicit none
! ---
integer :: i
! --- for fringe method
if ( inflow ) then
    l_fringe = l_fringe / z_i
    jx_relax  = int(l_fringe/dx)
    jx_fringe = nxt - jx_relax      ! jx_relax is readin from the param.nml

    gsize_fringe(1) = jx_relax
    gsize_fringe(2) = nyt
    gsize_fringe(3) = nzt-1
    lsize_fringe(1) = jx_relax
    lsize_fringe(2) = nynpy
    lsize_fringe(3) = nz-1
    start_fringe(1) = 0
    start_fringe(2) = nynpy * coordy
    start_fringe(3) = (nz-1) * coordz

    if ( read_inflow_file ) then
        allocate(wgt(jx_relax))
        lx_s  = lx_tot - l_fringe       ! Stevens et al. (2014) Renewable Energy
        lx_pl = lx_tot - 0.25_rprec*l_fringe
    
        do i = 1, jx_relax
            if ( x(jx_fringe + i) .lt. lx_pl ) then
                wgt(i) = 0.5_rprec*( 1._rprec - cos(pi*(x(jx_fringe + i) - lx_s) / (lx_pl - lx_s)) )
            else
                wgt(i) = 1._rprec    
            end if
        end do
    end if

    if ( rank == 0 ) then
        write(*,*) "jx_relax and jx_fringe are:", jx_relax, jx_fringe
    end if
end if
! ---
end subroutine fringe_prep
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine allocate_output_variable()
!-----------------------------------------------------------------------
!   output variables initialization (Bicheng Chen 06/13/2016)
!-----------------------------------------------------------------------
use io
use param !, only:nxt, nyt, nz, nzt, use_avgslice, average_dim_num, USE_MPI
use intermediate !, only:average_dim_select_flag, dim1_size, dim2_size, &
                 !     dim1_global, dim2_global, avg_out
implicit none
! ---
average_dim_select_flag = 1 - average_dim_num/2
dim1_size   = average_dim_select_flag*(nxt - nz + 1) + nz - 1
dim2_size   = average_dim_select_flag*(nz - 2) + 1
dim1_global = average_dim_select_flag*(nxt - nzt + 1) + nzt - 1
dim2_global = average_dim_select_flag*(nzt - 2) + 1
  
if (use_avgslice) then
    if (average_dim_num == 1) then
        allocate (avg_out(1:nxt, 1:nzt - 1))
    else if (average_dim_num == 2) then
        allocate (ax_1d(nz - 1))
        allocate (avg_out_1d(nzt - 1))
    end if
    allocate (ap(nxt, nz - 1), au(nxt, nz - 1), av(nxt, nz - 1), aw(nxt, nz - 1),       &
              p2(nxt, nz - 1), u2(nxt, nz - 1), v2(nxt, nz - 1), w2(nxt, nz - 1),       &
              auw(nxt, nz - 1), avw(nxt, nz - 1), acs(nxt, nz - 1), adudz(nxt, nz - 1), &
              advdz(nxt, nz - 1), aCs_Ssim(nxt, nz - 1), abeta_sgs(nxt, nz - 1),       &
              abetaclip_sgs(nxt, nz - 1), atxx(nxt, nz - 1), atxz(nxt, nz - 1),        &
              atyy(nxt, nz - 1), atyz(nxt, nz - 1), atzz(nxt, nz - 1),                 &
              u3(nxt, nz - 1), v3(nxt, nz - 1), w3(nxt, nz - 1),                       &
              adudt(nxt, nz - 1), advdt(nxt, nz - 1), adwdt(nxt, nz - 1))
    allocate (awu2(nxt, nz - 1), awv2(nxt, nz - 1), awT2(nxt, nz - 1),     &
              auv(nxt, nz - 1), awp(nxt, nz - 1), autau13(nxt, nz - 1),    &
              avtau23(nxt, nz - 1), adissip(nxt, nz - 1))
    
    if (theta_flag) then 
        allocate (atheta(nxt, nz - 1), t2(nxt, nz - 1), q2(nxt, nz - 1),       &
                  asgs_t3(nxt, nz - 1), awt(nxt, nz - 1), adTdz(nxt, nz - 1),  &
                  anu_t(nxt, nz - 1), t3(nxt, nz - 1), var_t(nxt, nz - 1))
    end if
    
    if (PCon_FLAG) then 
        allocate (aPCon(nxt, nz - 1), PCon2(nxt, nz - 1), asgs_PCon3(nxt, nz - 1),     &
                  awPCon(nxt, nz - 1), adPCondz(nxt, nz - 1), var_PCon(nxt, nz - 1),   &
                  aKc_t(nxt, nz - 1), aCs2Sc(nxt, nz - 1), asgs_PCon1(nxt, nz - 1),    & 
                  auPCon(nxt, nz - 1))
    end if
end if
  
if ( rank == 0 ) then
    allocate (ubar_avg(1, nzt - 1), vbar_avg(1, nzt - 1),             &
              thetabar_avg(1, nzt - 1), Cs2bar_avg(1, nzt - 1),       &
              Nutbar_avg(1, nzt - 1), ubar_tot(1, nzt - 1),           &
              vbar_tot(1, nzt - 1), thetabar_tot(1, nzt - 1),         &
              Cs2bar_tot(1, nzt - 1), Nutbar_tot(1, nzt - 1))
    allocate (upr2bar_avg(3, nzt - 1), stressbar_avg(3, nzt - 1),     &
              upr2bar_tot(3, nzt - 1), stressbar_tot(3, nzt - 1))
    allocate (Eozbar_avg(1, nxt/2, nzt - 1), Eozbar_tot(1, nxt/2, nzt - 1))
end if

! ---
end subroutine allocate_output_variable
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine allocate_flow_variable()
!-----------------------------------------------------------------------
!   flow variables initialization (Bicheng Chen 06/07/2016)
!! >>Add dynamical wind stress and dynamical Stokes drift
!! >>(Bicheng Chen 10/20/2016)
!-----------------------------------------------------------------------
use param
use sim_param
use fft
use test_filtermodule, only:G_test, G_test_test
use sgsmodule
use bottombc, only:ustar_avg, zo, z_os, patch, d0,          &
                   T_s, q_s, phi_m, psi_m, phi_h, psi_h,    &
                   zo_PCon, PCon_s
use topbc, only:sponge
use intermediate
use scalars_module
implicit none
! ---
integer :: jx, jy, ind

integer, dimension(npz) :: sendcounts, recvcounts, displs

real(rprec), dimension(:, :, :, :), allocatable :: k_ttlCR_tot
character(20), dimension(:), allocatable :: unt_bg      ! initial unit of background concentration
character(20) :: unt_to  
real(rprec), dimension(:, :), allocatable :: con_bg
real(rprec), dimension(:, :, :, :), allocatable :: con_bg_3d
real(rprec), dimension(:, :), allocatable :: k_bg       ! chemical reaction rate of specified species
integer, parameter :: fid_bg = 822
namelist/conprof/  unt_bg, unt_to, con_bg
namelist/chemrate/ k_bg
! ---
allocate (u(nxt, nynpy, 0:nz), v(nxt, nynpy, 0:nz), w(nxt, nynpy, 0:nz))
allocate (ke(nxt, nynpy, 0:nz), ke_temp(nxt, nynpy, 0:nz), &
          dragx(nxt, nynpy, 0:nz), dragy(nxt, nynpy, 0:nz), dragz(nxt, nynpy, 0:nz))
allocate (dudx(nxt, nynpy, 0:nz), dudy(nxt, nynpy, 0:nz), dudz(nxt, nynpy, 0:nz),        &
          dvdx(nxt, nynpy, 0:nz), dvdy(nxt, nynpy, 0:nz), dvdz(nxt, nynpy, 0:nz),        &
          dwdx(nxt, nynpy, 0:nz), dwdy(nxt, nynpy, 0:nz), dwdz(nxt, nynpy, 0:nz),        &
          RHSx(nxt, nynpy, 0:nz), RHSy(nxt, nynpy, 0:nz), RHSz(nxt, nynpy, 0:nz),        &
          RHSx_f(nxt, nynpy, 0:nz), RHSy_f(nxt, nynpy, 0:nz), RHSz_f(nxt, nynpy, 0:nz))
allocate (dudt(nxt, nynpy, 0:nz), dvdt(nxt, nynpy, 0:nz), dwdt(nxt, nynpy, 0:nz))
allocate (dkedx(nxt, nynpy, 0:nz), dkedy(nxt, nynpy, 0:nz), dkedz(nxt, nynpy, 0:nz))
allocate (txx(nxt, nynpy, 0:nz), txy(nxt, nynpy, 0:nz), txz(nxt, nynpy, 0:nz),           &
          tyy(nxt, nynpy, 0:nz), tyz(nxt, nynpy, 0:nz), tzz(nxt, nynpy, 0:nz))
allocate (divtx(nxt, nynpy, 0:nz), divty(nxt, nynpy, 0:nz), divtz(nxt, nynpy, 0:nz))
allocate (p(nxt, nynpy, 0:nz))
allocate (dpdx(nxt, nynpy, nz), dpdy(nxt, nynpy, nz), dpdz(nxt, nynpy, nz))

allocate (kx(lhx), ky(lhx))
allocate (kx_2d(lhx, nyt), ky_2d(lhx, nyt), k2(lhx, nyt))
allocate (kx_2d_mpi(nxhnpy, nyt), ky_2d_mpi(nxhnpy, nyt), k2_mpi(nxhnpy, nyt))
allocate (dealias(nxhnpy, nyt))
allocate (G_test(nxhnpy, nyt), G_test_test(nxhnpy, nyt))

allocate (Nu_t(nxt, nynpy, 0:nz), dissip(nxt, nynpy, 0:nz), magS(nxt, nynpy, 0:nz),     &
          u_lag(nxt, nynpy, 0:nz),v_lag(nxt, nynpy, 0:nz), w_lag(nxt, nynpy, 0:nz))
allocate (F_LM(nxt, nynpy, nz), F_MM(nxt, nynpy, nz), F_QN(nxt, nynpy, nz),             &
          F_NN(nxt, nynpy, nz), Beta(nxt, nynpy, nz), Betaclip(nxt, nynpy, nz))
allocate (Beta_avg(nz), Betaclip_avg(nz))
allocate (Cs_opt2(nxt, nynpy, nz), Cs_opt2_avg(nxt, nynpy, nz), Cs_Ssim(nxt, nynpy, nz), &
          epsilon_lag(nxt, nynpy, nz), epsilon_lag2(nxt, nynpy, nz),                  &
          xlag(nxt, nynpy, nz), ylag(nxt, nynpy, nz), zlag(nxt, nynpy, nz))
allocate (F_KX(nxt, nynpy, nz), F_XX(nxt, nynpy, nz), F_KX2(nxt, nynpy, nz), F_XX2(nxt, nynpy, nz))

allocate (ustar_avg(nxt, nynpy), zo(nxt, nynpy), z_os(nxt, nynpy), patch(nxt, nynpy), d0(nxt, nynpy))
allocate (sponge(0:nz))

allocate (cross_x(nxt, nynpy, nz), cross_y(nxt, nynpy, nz), cross_z(nxt, nynpy, nz))
allocate (vort1(nxt, nynpy, 0:nz), vort2(nxt, nynpy, 0:nz), vort3(nxt, nynpy, 0:nz))
allocate (beta_scal(nxt, nynpy, 0:nz))

allocate (cross_x_big(nxt2, ny2npy, nz), cross_y_big(nxt2, ny2npy, nz), cross_z_big(nxt2, ny2npy, nz))
allocate (u_big(nxt2, ny2npy, 0:nz), v_big(nxt2, ny2npy, 0:nz), w_big(nxt2, ny2npy, 0:nz))
allocate (vort1_big(nxt2, ny2npy, 0:nz), vort2_big(nxt2, ny2npy, 0:nz), vort3_big(nxt2, ny2npy, 0:nz))

! --- press_stag_array
allocate (rH_x(nxt, nynpy, 0:nz), rH_y(nxt, nynpy, 0:nz), rH_z(nxt, nynpy, 0:nz))
allocate (H_x(nxhnpy, nyt, 0:nz), H_y(nxhnpy, nyt, 0:nz), H_z(nxhnpy, nyt, 0:nz))
allocate (rtopw(nxt, nynpy), rbottomw(nxt, nynpy))
allocate (topw(nxhnpy, nyt), bottomw(nxhnpy, nyt))

! --- SGS model
if ( sgs ) then
    if (model == 3) then
        allocate (L11(nxt, nynpy), L12(nxt, nynpy), L13(nxt, nynpy),                        &
                  L22(nxt, nynpy), L23(nxt, nynpy), L33(nxt, nynpy),                        &
                  Q11(nxt, nynpy), Q12(nxt, nynpy), Q13(nxt, nynpy),                        &
                  Q22(nxt, nynpy), Q23(nxt, nynpy), Q33(nxt, nynpy), S_bar(nxt, nynpy),     &
                  S11_bar(nxt, nynpy), S12_bar(nxt, nynpy), S13_bar(nxt, nynpy),            &
                  S22_bar(nxt, nynpy), S23_bar(nxt, nynpy), S33_bar(nxt, nynpy),            &
                  S_S11_bar(nxt, nynpy), S_S12_bar(nxt, nynpy), S_S13_bar(nxt, nynpy),      &
                  S_S22_bar(nxt, nynpy), S_S23_bar(nxt, nynpy), S_S33_bar(nxt, nynpy),      &
                  S_hat(nxt, nynpy),       &
                  S11_hat(nxt, nynpy), S12_hat(nxt, nynpy), S13_hat(nxt, nynpy),            &
                  S22_hat(nxt, nynpy), S23_hat(nxt, nynpy), S33_hat(nxt, nynpy),            &
                  S_S11_hat(nxt, nynpy), S_S12_hat(nxt, nynpy), S_S13_hat(nxt, nynpy),      &
                  S_S22_hat(nxt, nynpy), S_S23_hat(nxt, nynpy), S_S33_hat(nxt, nynpy),      &
                  u_bar(nxt, nynpy), v_bar(nxt, nynpy), w_bar(nxt, nynpy),                  &
                  u_hat(nxt, nynpy), v_hat(nxt, nynpy), w_hat(nxt, nynpy), S(nxt, nynpy))
        allocate (beta_sd(nz))
        M11 => Q11; M12 => Q12; M13 => Q13; M22 => Q22; M23 => Q23; M33 => Q33
    else if (model == 4) then
        allocate (u_pr(nxt, nynpy), v_pr(nxt, nynpy), w_pr(nxt, nynpy))
        allocate (w_nod(nxt, nynpy), S(nxt, nynpy), tempos(nxt, nynpy),                 &    
                  L11(nxt, nynpy), L12(nxt, nynpy), L13(nxt, nynpy),                    &
                  L22(nxt, nynpy), L23(nxt, nynpy), L33(nxt, nynpy),                    &
                  Q11(nxt, nynpy), Q12(nxt, nynpy), Q13(nxt, nynpy),                    &
                  Q22(nxt, nynpy), Q23(nxt, nynpy), Q33(nxt, nynpy),                    &
                  M11(nxt, nynpy), M12(nxt, nynpy), M13(nxt, nynpy),                    &
                  M22(nxt, nynpy), M23(nxt, nynpy), M33(nxt, nynpy),                    &
                  N11(nxt, nynpy), N12(nxt, nynpy), N13(nxt, nynpy),                    &
                  N22(nxt, nynpy), N23(nxt, nynpy), N33(nxt, nynpy),                    & 
                  LM(nxt, nynpy), MM(nxt, nynpy), QN(nxt, nynpy), NN(nxt, nynpy),       &
                  Tn(nxt, nynpy), epsi(nxt, nynpy), dumfac(nxt, nynpy), S_bar(nxt, nynpy),  &
                  S11_bar(nxt, nynpy), S12_bar(nxt, nynpy), S13_bar(nxt, nynpy),        &
                  S22_bar(nxt, nynpy), S23_bar(nxt, nynpy), S33_bar(nxt, nynpy),        &
                  S_S11_bar(nxt, nynpy), S_S12_bar(nxt, nynpy), S_S13_bar(nxt, nynpy),  &
                  S_S22_bar(nxt, nynpy), S_S23_bar(nxt, nynpy), S_S33_bar(nxt, nynpy),  &
                  S_hat(nxt, nynpy),    & 
                  S11_hat(nxt, nynpy), S12_hat(nxt, nynpy), S13_hat(nxt, nynpy),        &
                  S22_hat(nxt, nynpy), S23_hat(nxt, nynpy), S33_hat(nxt, nynpy),        &
                  S_S11_hat(nxt, nynpy), S_S12_hat(nxt, nynpy), S_S13_hat(nxt, nynpy),  &
                  S_S22_hat(nxt, nynpy), S_S23_hat(nxt, nynpy), S_S33_hat(nxt, nynpy),  &
                  u_bar(nxt, nynpy), v_bar(nxt, nynpy), w_bar(nxt, nynpy),              &
                  u_hat(nxt, nynpy), v_hat(nxt, nynpy), w_hat(nxt, nynpy), fourbeta(nxt, nynpy))
        allocate (xp(nxt, nynpy, 0:nz), yp(nxt, nynpy, 0:nz), zp(nxt, nynpy, 0:nz),     &
                  u_temp(nxt, nynpy, 0:nz), v_temp(nxt, nynpy, 0:nz))
        allocate (FF_LM(nxt + 2, nynpy + 2, nz + 2), FF_MM(nxt + 2, nynpy + 2, nz + 2))
  
    else if (model == 5) then
        allocate (visc(nxt, nynpy, nz), Cs_opt2_2d(nxt, nynpy, nz), Cs_opt2_4d(nxt, nynpy, nz))
        allocate (S(nxt, nynpy), tempos(nxt, nynpy),                                        & 
                  L11(nxt, nynpy), L12(nxt, nynpy), L13(nxt, nynpy),                        &
                  L22(nxt, nynpy), L23(nxt, nynpy), L33(nxt, nynpy),                        &
                  Q11(nxt, nynpy), Q12(nxt, nynpy), Q13(nxt, nynpy),                        &
                  Q22(nxt, nynpy), Q23(nxt, nynpy), Q33(nxt, nynpy),                        &
                  M11(nxt, nynpy), M12(nxt, nynpy), M13(nxt, nynpy),                        &
                  M22(nxt, nynpy), M23(nxt, nynpy), M33(nxt, nynpy),                        &
                  N11(nxt, nynpy), N12(nxt, nynpy), N13(nxt, nynpy),                        &
                  N22(nxt, nynpy), N23(nxt, nynpy), N33(nxt, nynpy),                        & 
                  LM(nxt, nynpy), MM(nxt, nynpy), QN(nxt, nynpy), NN(nxt, nynpy),           &
                  Tn(nxt, nynpy), epsi(nxt, nynpy), dumfac(nxt, nynpy), S_bar(nxt, nynpy),  &
                  S11_bar(nxt, nynpy), S12_bar(nxt, nynpy), S13_bar(nxt, nynpy),            &
                  S22_bar(nxt, nynpy), S23_bar(nxt, nynpy), S33_bar(nxt, nynpy),            &
                  S_S11_bar(nxt, nynpy), S_S12_bar(nxt, nynpy), S_S13_bar(nxt, nynpy),      & 
                  S_S22_bar(nxt, nynpy), S_S23_bar(nxt, nynpy), S_S33_bar(nxt, nynpy),      &
                  S_hat(nxt, nynpy),        & 
                  S11_hat(nxt, nynpy), S12_hat(nxt, nynpy), S13_hat(nxt, nynpy),            &
                  S22_hat(nxt, nynpy), S23_hat(nxt, nynpy), S33_hat(nxt, nynpy),            &
                  S_S11_hat(nxt, nynpy), S_S12_hat(nxt, nynpy), S_S13_hat(nxt, nynpy),      &
                  S_S22_hat(nxt, nynpy), S_S23_hat(nxt, nynpy), S_S33_hat(nxt, nynpy),      &  
                  u_bar(nxt, nynpy), v_bar(nxt, nynpy), w_bar(nxt, nynpy),                  &
                  u_hat(nxt, nynpy), v_hat(nxt, nynpy), w_hat(nxt, nynpy))
        allocate (LMvert(nz), MMvert(nz), QNvert(nz), NNvert(nz))
        allocate (xp(nxt, nynpy, 0:nz), yp(nxt, nynpy, 0:nz), zp(nxt, nynpy, 0:nz),         &
                  u_temp(nxt, nynpy, 0:nz), v_temp(nxt, nynpy, 0:nz))
        allocate (FF_LM(nxt + 2, nynpy + 2, nz + 2), FF_MM(nxt + 2, nynpy + 2, nz + 2),     &
                  FF_QN(nxt + 2, nynpy + 2, nz + 2), FF_NN(nxt + 2, nynpy + 2, nz + 2),     &
                  Beta_t(nxt + 2, nynpy + 2, nz + 2)) 
    end if
end if

! --- Temperature field
allocate (phi_m(nxt, nynpy), psi_m(nxt, nynpy), phi_h(nxt, nynpy), psi_h(nxt, nynpy))
allocate (L(nxt, nynpy), wstar(nxt, nynpy))

if (theta_flag) then
    allocate (T_s(nxt, nynpy), q_s(nxt, nynpy))
    allocate (theta(nxt, nynpy, 0:nz), q(nxt, nynpy, 0:nz))
    allocate (Pr_(nxt, nynpy, 0:nz), dTdz(nxt, nynpy, 0:nz), dqdz(nxt, nynpy, 0:nz),         &
              RHS_Tf(nxt, nynpy, 0:nz), RHS_T(nxt, nynpy, 0:nz), RHS_qf(nxt, nynpy, 0:nz),   &
              RHS_q(nxt, nynpy, 0:nz), sgs_t3(nxt, nynpy, 0:nz), sgs_q3(nxt, nynpy, 0:nz))    
    allocate (T_s_filtered(nxt, nynpy))
    
    ! - SGS model
    if (model_sc == 4) then
        allocate (FF_KX(nxt + 2, nynpy + 2, nz + 2), FF_XX(nxt + 2, nynpy + 2, nz + 2))
    else if (model_sc == 5) then
        allocate (FF_KX(nxt + 2, nynpy + 2, nz + 2), FF_XX(nxt + 2, nynpy + 2, nz + 2), &
                  FF_KX2(nxt + 2, nynpy + 2, nz + 2), FF_XX2(nxt + 2, nynpy + 2, nz + 2))
    end if
    !- use dynamical kinematic heat flux (Bicheng Chen)
    if (flag_dynWT) then
        inquire (iolength=len_wt) wt_s
        open (unit=fid_wt, file=fn_wt, access='direct', recl=len_wt)
        read (unit=fid_wt, rec=1) wt_s
        close (fid_wt)
    end if
end if

! --- concentration variables initialization
allocate (beta_pcon(nxt, 0:nynpy+1, 0:nz))  ! allocate the variable always needed
if (pcon_flag) then
    npcon = info_con%n_con
    allocate (zo_PCon(nxt, nynpy), PCon_s(nxt, nynpy))
    allocate (P_surf_flux(nxt, nynpy), deposition(nxt, nynpy),  &
              P_surf_flux_dep(nxt, nynpy), Real_dep(nxt, nynpy))
    allocate (matrix_x(nxt, nxt), matrix_y(nyt, nyt))
    allocate (dvector_x(nxt), dvector_y(nyt))

    if (settling) then
        allocate (settling_vel(npcon), densratio_pcon(npcon), v_pcon(npcon))
        settling_vel = info_con%vel_settling(1:npcon)/u_scale
        densratio_pcon = info_con%ratio_dens(1:npcon)
        v_pcon = info_con%vol_spec(1:npcon)*pcon_scale
    end if

    allocate (Kc_t(nxt, 0:nynpy+1, 0:nz), dPCondz(nxt, 0:nynpy+1, 0:nz),          &
              sgs_PCon3(nxt, 0:nynpy+1, 0:nz), res_PCon3(nxt, 0:nynpy+1, 0:nz),   &
              sgs_PCon1(nxt, 0:nynpy+1, 0:nz), Cs2Sc(nxt, 0:nynpy+1, 0:nz))
    
    !- allocate concentration
    allocate (PCon(nxt, 0:nynpy+1, 0:nz, npcon))
    allocate (RHS_PConf(nxt, 0:nynpy+1, 0:nz, npcon), RHS_PCon(nxt, 0:nynpy+1, 0:nz, npcon))

    call read_source_info()
end if

! ---
Cs_opt2_avg = 0._rprec
! ---
end subroutine allocate_flow_variable
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_vel_field()
!-----------------------------------------------------------------------
!   initialize velocity field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
implicit none
! ---
if ( restart_vel ) then
    if ( inflow .and. read_inflow_file .and. interp ) then
        call read_precursor_field()
    else 
        call read_flow_field()
    end if 
else
    if ( read_ic_profile ) then
        call ic_read()
    else
        if (theta_flag) then
            if ( intp_vel ) then
                call create_fft_plan_intp()
                call read_coarse_field()
                call intp_flow_field()
            else
                call ic_scal()
            end if
        else
            call ic()
        end if
    end if
end if

!- add cross velocity to the all velocity field (Bicheng Chen 09/07/2016)
if (flag_crossVel .and. ocean_flag) then
    u(1:nxt, 1:nynpy, 0:nz) = u(1:nxt, 1:nynpy, 0:nz) + u_cross
    v(1:nxt, 1:nynpy, 0:nz) = v(1:nxt, 1:nynpy, 0:nz) + v_cross
end if

! ---
end subroutine init_vel_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_temperature_field()
!-----------------------------------------------------------------------
!   initialize temperature field
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use scalars_module, only:RHS_T, sgs_t3
use bottombc, only:psi_m
implicit none
! ---
integer(MPI_OFFSET_KIND) :: disp0, disp2d, disp3d
real(rprec), dimension(nxt/2, nynpy, nz-1) :: theta_tmp, RHS_T_tmp
real(rprec), dimension(nxt/2, nynpy) :: sgs_t3_tmp, psi_m_tmp
! ---
character(80) :: fname
logical :: exst
! ---
if (theta_flag .and. restart_theta .and. restart_vel) then

    if ( inflow .and. read_inflow_file .and. interp ) then    
        write (fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        
        disp0 = 0
        disp2d = npy * sizeof(sgs_t3_tmp(1:nxt/2, 1:nynpy))
        disp3d = np  * sizeof(theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1))

        call read_file_mpi_3d_inflow(theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_3d_inflow(RHS_T_tmp(1:nxt/2, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_xy_inflow(sgs_t3_tmp(1:nxt/2, 1:nynpy), fname, disp0)
        disp0 = disp0 + disp2d
        call read_file_mpi_xy_inflow(psi_m_tmp, fname, disp0)

        theta(1:nxt/2, 1:nynpy, 1:nz-1) = theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        RHS_T(1:nxt/2, 1:nynpy, 1:nz-1) = RHS_T_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        sgs_t3(1:nxt/2, 1:nynpy, 1) = sgs_t3_tmp(1:nxt/2, 1:nynpy)
        psi_m(1:nxt/2, 1:nynpy) = psi_m_tmp(1:nxt/2, 1:nynpy)

        theta(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = theta_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        RHS_T(nxt/2+1:nxt, 1:nynpy, 1:nz-1) = RHS_T_tmp(1:nxt/2, 1:nynpy, 1:nz-1)
        sgs_t3(nxt/2+1:nxt, 1:nynpy, 1) = sgs_t3_tmp(1:nxt/2, 1:nynpy)
        psi_m(nxt/2+1:nxt, 1:nynpy) = psi_m_tmp(1:nxt/2, 1:nynpy)
    else
        write (fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        disp0 = 0
        disp2d = npy * sizeof(sgs_t3(1:nxt, 1:nynpy, 1))
        disp3d = np  * sizeof(theta(1:nxt, 1:nynpy, 1:nz-1))

        call read_file_mpi_3d(theta(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_3d(RHS_T(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0 + disp3d
        call read_file_mpi_xy(sgs_t3(1:nxt, 1:nynpy, 1), fname, disp0)
        disp0 = disp0 + disp2d
        call read_file_mpi_xy(psi_m, fname, disp0)
    end if

else if (theta_flag .and. restart_theta .and. (.not. restart_vel)) then !LV1
    print *, "Cannot initialize temperature field with initializing velocity."
    print *, "Stop"
    stop
else if ((.not. theta_flag) .and. restart_theta .and. restart_vel) then !LV1
    print *, "Cannot initialize temperature field with theta_flag."
    print *, "Stop"
    stop
else
    !- the initialization of temperature without file is in ic_scal
end if
! ---
end subroutine init_temperature_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_concentration_field()
!-----------------------------------------------------------------------
!   concentration variables initialization
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
use scalars_module, only:deposition, Real_dep, RHS_PCon, Kc_t
use sgsmodule, only:F_KX, F_XX, F_KX2, F_XX2
implicit none
! ---
integer(MPI_OFFSET_KIND) :: disp0, disp2d, disp3d
! ---
character(80) :: fname
logical :: exst
integer :: ipcon
! ---
if ( PCon_flag .and. inflow .and. write_inflow_file ) then
    call nutrient_inflow_read()
end if

if (PCon_flag .and. restart_pcon) then
    write (fname, '(a,i8.8,a)') path_restart//'con_tt', nums, '.out'
    disp0 = 0
    disp2d = npy * sizeof(deposition)
    disp3d = np  * sizeof(PCon(1:nxt, 1:nynpy, 1:nz-1, 1))
    do ipcon = 1, npcon
        call read_file_mpi_3d(PCon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)
        call read_file_mpi_3d(RHS_Pcon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0+npcon*disp3d)
        disp0 = disp0 + disp3d
    end do

    write (fname, '(a,i8.8,a)') path_restart//'con_diff_tt', nums, '.out'
    call read_file_mpi_xy(deposition, fname, 0)
    disp0 = disp2d
    call read_file_mpi_xy(Real_dep, fname, disp0)
    disp0 = disp0+disp2d
    call read_file_mpi_3d(Kc_t(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_KX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_XX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_KX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call read_file_mpi_3d(F_XX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
else if ((.not. pcon_flag) .and. restart_pcon) then
    print *, "Cannot initialize concentration field without pcon_flag. Stop."
    stop
else if (PCon_flag .and. .not. Restart_pcon) then
    !- Always initialize pollen concentration with zeros
    PCon(:, :, :, :) = 0._rprec
    if ( coordz == 0 ) then
        PCon(:, :, 0, :) = BOGUS
    end if
end if
! ---
end subroutine init_concentration_field
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine read_source_info
!-----------------------------------------------------------------------
!   initialize point sources and single sources (Bicheng Chen 05/22/2015)
!-----------------------------------------------------------------------
use param
use scalars_module
implicit none
! ---
character(80) :: fname
integer :: ipcon
! ---
if (pointsource) then
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write(fname, "(A,I0,A)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        !--- read in temp variable
        open(fid_src, file=fname, form='formatted', status='old')
        read(fid_src, nml=number_source)
        !print *, 'n_src ', n_src
        read(fid_src, nml=source)
        close(fid_src)
        !--- copy to real source variable
        con_src(ipcon)%n = n_src
        con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
        con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
        con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
        con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
    end do
  
    ! if (.not. single_point) then
    !     do ipcon = 1, npcon
    !         if (con_src(ipcon)%n /= 1) then
    !             print *, "the number of source is not 1 for non-single point sources."
    !             stop
    !         end if
    !     end do
    ! end if
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_y = int(con_src(ipcon)%iy/nynpy)
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iy_cpu = con_src(ipcon)%iy - con_src(ipcon)%icpu_y*nynpy
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do

else if (flag_1Dsrc) then
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write(fname, "(A,I0,A)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        !--- read in temp variable
        open(fid_src, file=fname, form='formatted', status='old')
        read(fid_src, nml=number_source)
        read(fid_src, nml=source_alignment)
        read(fid_src, nml=source)
        close(fid_src)
        !--- copy to real source variable
        select case (src_align)
        case ('lateral')
            con_src(ipcon)%n = n_src
            con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
            ! con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
            con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
            con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
            ! con_src(ipcon)%xs(1:n_src) = x_src(1:n_src)
            ! ! con_src(ipcon)%ys(1:n_src) = y_src(1:n_src)
            ! con_src(ipcon)%zs(1:n_src) = z_src(1:n_src)
            ! con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
        case ('streamwise')
            con_src(ipcon)%n = n_src
            ! con_src(ipcon)%ix(1:n_src) = ix_src(1:n_src)
            con_src(ipcon)%iy(1:n_src) = iy_src(1:n_src)
            con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
            con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
            ! ! con_src(ipcon)%xs(1:n_src) = x_src(1:n_src)
            ! con_src(ipcon)%ys(1:n_src) = y_src(1:n_src)
            ! con_src(ipcon)%zs(1:n_src) = z_src(1:n_src)
            ! con_src(ipcon)%rls(1:n_src) = rls_src(1:n_src)
        case default
            write (*, *) 'invalid source alignment type'
            stop
        end select
    end do

    !-- calculate the specific cpu coord for point sources
    !-- and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_y = int(con_src(ipcon)%iy/nynpy)
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iy_cpu = con_src(ipcon)%iy - con_src(ipcon)%icpu_y*nynpy
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
        con_src(ipcon)%rls = con_src(ipcon)%rls/(dx*z_i*dy*z_i*dz*z_i)/(PCon_scale/(z_i/u_scale))
    end do

else if (flag_2Dsrc) then
    !-- allocate variable
    allocate (con_src(npcon))
    !-- read source release data
    do ipcon = 1, info_con%n_con
        write (fname, "(A,I0,A)") trim(fnPre_src), ipcon, trim(fnSuf_src)
        !--- read in temp variable
        open(fid_src, file=fname, form='formatted', status='old')
        read(fid_src, nml=number_source)
        !print *, 'n_src ', n_src
        read(fid_src, nml=source)
        close(fid_src)
        !--- copy to real source variable
        con_src(ipcon)%n = n_src
        con_src(ipcon)%iz(1:n_src) = iz_src(1:n_src)
    end do

    !-- calculate the specific cpu coord for point sources
    !-- and normalize the source release
    do ipcon = 1, npcon
        con_src(ipcon)%icpu_z = int((con_src(ipcon)%iz - 1)/(nz - 1))
        con_src(ipcon)%iz_cpu = con_src(ipcon)%iz - con_src(ipcon)%icpu_z*(nz - 1)
    end do
end if
! ---
end subroutine read_source_info
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------