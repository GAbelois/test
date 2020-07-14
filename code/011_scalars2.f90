!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
module scalars_module2
!--------------------------------------------------------------------!
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter attribution (Bicheng Chen 06/12/2016)
!--------------------------------------------------------------------!    
use types, only:rprec
use param !, jt_global => jt  !--rename to avoid name clashes
 !--couldx also modify all routines to access jt
 !  from param module, not argument list
use bottombc !makes obukhov functions available
use sim_param, only:u, v, w, theta, q, pcon
use sgsmodule, only:Nu_t
use scalars_module, only:L, wstar, dTdz, dqdz, sgs_t3, sgs_q3,          &
                         DPCondz, sgs_PCon3, res_PCon3, Kc_t, Cs2Sc,    &
                         sgs_PCon1
use intermediate, only:average_dim_select_flag, dim1_size, dim2_size,   &
                       dim1_global, dim2_global
implicit none

contains

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------    
    subroutine scalar_in()
    !subroutine scalar_in(T_s,zo)
    !subroutine scalar_in(T_s,zo,suffix2,suffix3)
    !c This reads in the scalar input from a interpolated scalar
    !c file and assigns it to the x-y plane at z=0 i.e. at the ground
    !c This for the moment is used to read the momentum roughness, z0
    !c and temperature from the USDA remote-sensed data set. The data
    !c can be interpolated using bilinear/cubic/spline interpolation (MATLAB)
    !c for the grid size(nxt,nyt)
    !c Authored by Vijayant Kumar
    !c Last modified on April 11th, 2004

    !use param
    use bottombc, only:T_s, zo !Load the variables from bottombc and update in here
    implicit none
    ! ---
    integer :: ii, jj
    character(len=6) :: suffix2, suffix3
    character(*), parameter :: fmt_5169 = "(1400(E17.10))"
    write (suffix2, '(i6.6)') nxt ! create a string from nxt
!       call append_zeros_string(suffix2,nxt) ! Append leading zeros to string

    if (coarse_grain_flag) then
        write (suffix3, '(i6.6)') stencil_pts
!      call append_zeros_string(suffix3,stencil_pts)
        open (unit=77, file='../interp_data/coarse_grained/interp_temp_cg_' &
            //suffix2(4:6)//'pts_'//suffix3(4:6)//'.out', status='unknown')
        open (unit=78, file='../interp_data/coarse_grained/interp_z_m_cg_' &
            //suffix2(4:6)//'pts_'//suffix3(4:6)//'.out', status='unknown')

        print *, 'interp_temp_cg_'//suffix2(4:6)//'pts_' &
            //suffix3(4:6)//'.out loaded from scalar_in.f'

        do jj = 1, nyt
            read (77, fmt_5169) (T_s(ii, jj), ii=1, nxt)
            read (78, fmt_5169) (zo(ii, jj), ii=1, nxt)
        end do
        T_s = T_s + 273.15_rprec ! Convert T to Kelvin
        close (77)
        close (78)
    else

!open(unit=77,file='./interp_data/interp_temp_'//suffix2(4:6)&
!//'X'//suffix2(4:6)//'_cubic.out',status='unknown')
!open(unit=78,file='./interp_data/interp_z_m_'//suffix2(4:6)&
!//'X'//suffix2(4:6)//'_cubic.out',status='unknown')

        open (unit=77, file=path//'interp_data/interp_temp_'//suffix2(4:6) &
            //'X'//suffix2(4:6)//'_cubic.out', status='unknown')
        open (unit=78, file=path//'interp_data/interp_z_m_'//suffix2(4:6) &
            //'X'//suffix2(4:6)//'_cubic.out', status='unknown')
        print *, path//'interp_data/interp_temp_'//suffix2(4:6)//'X'//suffix2(4:6)//'_cubic.out'

        do jj = 1, nyt
!          do ii=1,nxt
            read (77, fmt_5169) (T_s(ii, jj), ii=1, nxt)
            read (78, fmt_5169) (zo(ii, jj), ii=1, nxt)
        end do
        print *, 'Mean T_s (K) = ', sum(T_s)/float(nxt*nyt)
        print *, 'Mean zo (m) = ', sum(zo)/float(nxt*nyt)
        close (77)
        close (78)
        T_s = T_s + 273.15_rprec ! Convert T to Kelvin
    end if
    ! ---
    end subroutine scalar_in
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine scalar_slice()
    !subroutine scalar_slice(w,t,q,sgs_t3,sgs_q3,dTdz,beta_scal,Nu_t,jt)
    !c This is exactly the same like the subroutine avgslice with the
    !c only difference being that it averages the scalar variables
    !c to find the y-averaged instantaneous x-z slices of variables
    !c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
    !c It also outputs the average covariance between wt and wq
    !use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3
    !use output_slice,only: collocate_MPI_averages
    use intermediate, only: atheta, t2, asgs_t3, awt, adTdz, anu_t, t3, var_t, &
                            avg_out, ax_1d, avg_out_1d
    implicit none
    integer:: i, j, k
    !real(kind=rprec),dimension(nxt,nz-1),save:: atheta,t2,q2,asgs_t3,awt
    !real(kind=rprec),dimension(nxt,nz-1),save:: adTdz,anu_t,t3,var_t
    !real(kind=rprec),dimension(nz-1),save::ax_1d
    !!real(kind=rprec),dimension(nxt,nz-1),save:: aq,adqdz,awq,q2,asgs_q3
    real(kind=rprec) :: ttheta1, tt2, tsgst, twt, tdTdz, fr
    real(kind=rprec) :: tnu_t, tt3
    real(kind=rprec), dimension(nynpy) :: arg1 
    !real(kind=rprec):: tq1,tq2,tsgsq,twq,tdqdz
    !real(kind=rprec),dimension(:,:),allocatable:: avg_out
    !real(kind=rprec),dimension(:),allocatable::avg_out_1d
    real(kind=rprec), dimension(0:nz-1) :: Tbar_profile
    character(*), parameter :: fmt_5168 = "(1400(E14.5))"

    fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)

    !if (jt .EQ. theta_init_time) then
    if (jt .eq. c_count) then
        atheta = 0._rprec; t2 = 0._rprec; asgs_t3 = 0._rprec
        awt = 0._rprec; adTdz = 0._rprec
        anu_t = 0._rprec; t3 = 0._rprec; var_t = 0._rprec
        !aq=0._rprec;q2=0._rprec;asgs_q3=0._rprec;awq=0._rprec;adqdz=0._rprec
    end if

    call calc_hor_avg2(theta(:,:,1:nz), Tbar_profile)


    do k = 1, nz - 1
    do i = 1, nxt
        ttheta1 = 0._rprec; tt2 = 0._rprec; tsgst = 0._rprec
        twt = 0._rprec; tdTdz = 0._rprec
        tnu_t = 0._rprec; tt3 = 0._rprec
        !tq1=0._rprec;tq2=0._rprec;tsgsq=0._rprec;twq=0._rprec;tdqdz=0._rprec

        do j = 1, nynpy
            if ( (k .eq. 1) .and. rank == 0 ) then
                arg1(j) = 0._rprec
                !         arg2=0._rprec
            else
                arg1(j) = (theta(i, j, k) + theta(i, j, k - 1))/2.
                !         arg2=(q(i,j,k)+q(i,j,k-1))/2.
            end if
        end do

        call avg_y(theta(i, :, k), ttheta1)
        call avg_y(sgs_t3(i, :, k), tsgst)
        call avg_y(dTdz(i, :, k), tdTdz)
        call avg_y(Nu_t(i, :, k), tnu_t)
        call avg_y(w(i, :, k)*arg1(:), twt)

        !var_t(i, k) = var_t(i, k) + fr*sum((theta(1:nxt, 1:nynpy, k) - sum(theta(1:nxt, 1:nynpy, k))/(nxt*nynpy))**2)/(nxt*nynpy)

        atheta(i, k) = atheta(i, k) + (fr)*ttheta1
        !t2(i,k)=t2(i,k)+(fr)*tt2
        asgs_t3(i, k) = asgs_t3(i, k) + (fr)*tsgst
        awt(i, k) = awt(i, k) + (fr)*twt
        adTdz(i, k) = adTdz(i, k) + (fr)*tdTdz
        anu_t(i, k) = anu_t(i, k) + (fr)*tnu_t
        !t3(i,k)=t3(i,k)+(fr)*tt3/nynpy
        
        ! Humidity variables
        !aq(i,k)=aq(i,k)+(fr)*tq1/nynpy
        !q2(i,k)=q2(i,k)+(fr)*tq2/nynpy
        !asgs_q3(i,k)=asgs_q3(i,k)+(fr)*tsgsq/nynpy
        !awq(i,k)=awq(i,k)+(fr)*twq/nynpy
        !adqdz(i,k)=adqdz(i,k)+(fr)*tdqdz/nynpy
        ! Humidity variables
    end do
    end do

    if (mod(jt, p_count) == 0) then
        if (average_dim_num == 1) then
            call collocate_MPI_averages_N(atheta, avg_out, 35, 'theta')
            call collocate_MPI_averages_N(t2, avg_out, 36, 't2')
            call collocate_MPI_averages_N(asgs_t3, avg_out, 37, 'sgs_t3')
            call collocate_MPI_averages_N(awt, avg_out, 38, 'wt')
            call collocate_MPI_averages_N(adTdz, avg_out, 39, 'dTdz')
            !    call collocate_MPI_averages_N(aq,avg_out,40,'q')
            !    call collocate_MPI_averages_N(q2,avg_out,41,'q2')
            !    call collocate_MPI_averages_N(adqdz,avg_out,42,'dqdz')
            !    call collocate_MPI_averages_N(asgs_q3,avg_out,43,'sgs_q3')
            !    call collocate_MPI_averages_N(awq,avg_out,44,'wq')
            call collocate_MPI_averages_N(anu_t, avg_out, 45, 'Nu_t')
            call collocate_MPI_averages_N(t3, avg_out, 46, 't3')
            !call collocate_MPI_averages_N(var_t, avg_out, 47, 'var_t')
        else if (average_dim_num == 2) then
            ax_1d = sum(atheta, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 35, 'theta')
            ax_1d = sum(t2, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 36, 't2')
            ax_1d = sum(asgs_t3, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 37, 'sgs_t3')
            ax_1d = sum(awt, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 38, 'wt')
            ax_1d = sum(adTdz, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 39, 'dTdz')
            ax_1d = sum(anu_t, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 45, 'Nu_t')
            ax_1d = sum(t3, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 46, 't3')
            !ax_1d = sum(var_t, dim=1)/nxt
            !call collocate_MPI_averages_N(ax_1d, avg_out_1d, 47, 'var_t')
        end if

        atheta = 0._rprec; t2 = 0._rprec; asgs_t3 = 0._rprec
        awt = 0._rprec; adTdz = 0._rprec;
        anu_t = 0._rprec; t3 = 0._rprec;
        !aq=0._rprec;q2=0._rprec;asgs_q3=0._rprec;awq=0._rprec;adqdz=0._rprec
        !var_t = 0._rprec
    end if
    ! ---
    end subroutine scalar_slice
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine budget_TKE_scalar
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! To compute the budget of TKE, we need the following terms:
    ! <u'w'>*d/dz(<U>), <v'w'>*d/dz(<V>), (g/<T>)*<w'T'>, [0.5*d/dz(<w'e>) where
    ! e=<u'^2>+<v'^2>+<w'^2>], d/dz(<w'p'>), dissipation (Nu_t*S^2),d/dz(<u'\tau{13}>)+
    ! d/dz(v'\tau{23}), <\tau{13}>d/dz(<U>)+<\tau{23}>d/dz(<V>)
    ! Of the eight terms above, we can directly compute terms 1,2,3,8 from the variables
    ! calculated/outputted in avgslice.f90 and scalar_slice.f90
    ! So, the rest 4 terms will be computed/outputted here
    ! Similarly for temperature variance budget, we need
    ! <w'T'>*d/dz(<T>), d/dz(<w*T^2>) and we already have term 1 from
    ! scalar_slice.f90. so we just need to compute term 2 here
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use sim_param, only:u, v, w, theta, p, txz, tyz
    use param, only:path, p_count, c_count, jt
    use sgsmodule, only:dissip
    use intermediate, only: awu2, awv2, awT2, auv, awp, autau13, avtau23,   &
                            adissip, ax_1d, avg_out_1d, avg_out
    use stokes_drift, only: ust, vst
    implicit none
    ! ---
    integer :: i, j, k, jz
    real(kind=rprec) :: twu2,twv2,twt2,tuv,twp,fr
    real(kind=rprec), dimension(nynpy) :: arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8
    real(kind=rprec) :: tutau13,tvtau23,tdissip
    real(kind=rprec), dimension(0:nz-1) :: ubar_profile,vbar_profile,Tbar_profile

    ! --- added by cyan
    real(kind=rprec), dimension(nxt, nynpy, 0:nz) :: p_org
    real(kind=rprec), dimension(0:nz-1) :: pbar_profile
    ! ---
    fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)

    do jz = 1, nz-1
        p_org(:, :,jz) = p(:, :,jz) - 0.5*(u(:, :, jz)*u(:, :, jz) +   &
                                           v(:, :, jz)*v(:, :, jz) +   &
                                           w(:, :, jz)*w(:, :, jz))    &
                                    -     (u(:, :, jz)*ust(jz)+v(:, :, jz)*vst(jz))   &
                                    - 0.5*(ust(jz)**2 + vst(jz)**2)
    end do
    
    call calc_hor_avg2(u(:,:,1:nz), ubar_profile)
    call calc_hor_avg2(v(:,:,1:nz), vbar_profile)
    call calc_hor_avg2(p_org(:,:,1:nz), pbar_profile)
    if (theta_flag) then
        call calc_hor_avg2(theta(:,:,1:nz), Tbar_profile)
    end if
    
    do k = 1, nz-1
    do i = 1, nxt
        twu2=0._rprec;twv2=0._rprec;twT2=0._rprec;tuv=0._rprec;twp=0._rprec;
        tutau13=0._rprec;tvtau23=0._rprec;tdissip=0._rprec
    
        do j = 1, nynpy
            if( k .eq. 1 .and. coordz == 0 ) then
                arg1 = 0._rprec; arg2 = 0._rprec; arg3 = 0._rprec
                arg4 = 0._rprec; arg5 = 0._rprec; arg6 = 0._rprec
                arg7 = 0._rprec
            else
    !           arg1(j)=(u(i,j,k)*u(i,j,k)+u(i,j,k-1)*u(i,j,k-1))/2.
    !           arg2(j)=(v(i,j,k)*v(i,j,k)+v(i,j,k-1)*v(i,j,k-1))/2.
    !           arg3(j)=(theta(i,j,k)*theta(i,j,k)+theta(i,j,k-1)*theta(i,j,k-1))/2.
                ! --- commented by cyan
                ! arg4(j)=(p(i,j,k)+p(i,j,k-1))/2.  
                arg4(j)=(p_org(i,j,k)-pbar_profile(k)+p_org(i,j,k-1)-pbar_profile(k-1))/2.
                ! ---/cyan
                arg5(j)=(u(i,j,k)-ubar_profile(k)+u(i,j,k-1)-ubar_profile(k-1))/2.
                arg6(j)=(v(i,j,k)-vbar_profile(k)+v(i,j,k-1)-vbar_profile(k-1))/2.
                if (theta_flag) then
                    arg7(j)=(theta(i,j,k)-Tbar_profile(k)+theta(i,j,k-1)-Tbar_profile(k-1))/2.
                end if
                arg8(j) = (u(i,j,k)-ubar_profile(k))*(v(i,j,k)-vbar_profile(k))
            end if
        end do
        
        call avg_y(w(i,:,k)*arg5(:)*arg5(:), twu2)  !directly computes <w(u')^2>
        call avg_y(w(i,:,k)*arg6(:)*arg6(:), twu2)  !directly computes <w(v')^2>
        ! Also note that since <w> = 0, there is no need to calculate it as
        ! <w^3> = <w'^3> and we are outputting <w^3> in avgslice
        if (theta_flag) call avg_y(w(i,:,k)*arg7(:), twT2)
        call avg_y(w(i,:,k)*arg4(:), twp) 
        ! <u'v'> is not as simple as <u'w'> since <u> .ne. 0 whereas <w>=0
        ! therefore, we directly calculate <u'v'> here
        call avg_y(arg8(:), tuv)
        call avg_y(arg5(:)*txz(i,:,k), tutau13) ! computes SGS transport of TKE i.e. <u'\tau_{13}>
        call avg_y(arg6(:)*tyz(i,:,k), tvtau23) ! computes SGS transport of TKE i.e. <v'\tau_{13}>
        call avg_y(dissip(i,:,k), tdissip)      ! outputs dissip calculated in sgs stag..
        
        awu2(i,k) = awu2(i,k)+(fr)*twu2
        awv2(i,k) = awv2(i,k)+(fr)*twv2
        awT2(i,k) = awT2(i,k)+(fr)*twT2
        awp(i,k)  = awp(i,k)+(fr)*twp
        auv(i,k)  = auv(i,k)+(fr)*tuv
        autau13(i,k) = autau13(i,k)+(fr)*tutau13
        avtau23(i,k) = avtau23(i,k)+(fr)*tvtau23
        adissip(i,k) = adissip(i,k)+(fr)*tdissip
    end do
    end do
    
    if (mod(jt,p_count)==0) then
        if (average_dim_num == 1) then
            call collocate_MPI_averages_N(awu2,avg_out,61,'wu2')
            call collocate_MPI_averages_N(awv2,avg_out,62,'wv2')
            call collocate_MPI_averages_N(awT2,avg_out,63,'wT2')
            call collocate_MPI_averages_N(awp,avg_out,64,'wp')
            call collocate_MPI_averages_N(auv,avg_out,65,'uv')
            call collocate_MPI_averages_N(autau13,avg_out,66,'utau13')
            call collocate_MPI_averages_N(avtau23,avg_out,67,'vtau23')
            call collocate_MPI_averages_N(adissip,avg_out,68,'dissip')
        else if(average_dim_num == 2) then
            ax_1d=sum(awu2,dim=1)/nxt            
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,61,'wu2')
            ax_1d=sum(awv2,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,62,'wv2')
            ax_1d=sum(awT2,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,63,'wT2')
            ax_1d=sum(awp,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,64,'wp')
            ax_1d=sum(auv,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,65,'uv')
            ax_1d=sum(autau13,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,66,'utau13')
            ax_1d=sum(avtau23,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,67,'vtau23')
            ax_1d=sum(adissip,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,68,'dissip')
        end if
        !VK Zero out the outputted averages !!
        awu2=0._rprec;awv2=0._rprec;awT2=0._rprec;awp=0._rprec;auv=0._rprec
        autau13=0._rprec;avtau23=0._rprec;adissip=0._rprec
    end if
    ! ---
    end subroutine budget_TKE_scalar
!--------------------------------------------------------------------!
!                                               
!--------------------------------------------------------------------!    
    subroutine pollen_slice()
    !c This is exactly the same like the subroutine avgslice with the
    !c only difference being that it averages the scalar variables
    !c to find the y-averaged instantaneous x-z slices of variables
    !c t,q,sgs_t3,sgs_q3 and their variances such as t2, q2.
    !c It also outputs the average covariance between wt and wq
    !use scalars_module,only: dTdz,dqdz,sgs_t3,sgs_q3
    !use output_slice,only: collocate_MPI_averages
    use intermediate, only:aPCon, PCon2, asgs_PCon3, awPCon, adPCondz, var_PCon, &
                aKc_t, aCs2Sc, asgs_PCon1, auPCon, avg_out, ax_1d, avg_out_1d
    implicit none
    ! ---
    integer:: i, j, k
    real(rprec):: tPCon, tPCon2, tsgs_PCon3, twPCon, tdPCondz, fr
    real(rprec):: tKc_t, tCs2Sc, tmp
    real(rprec):: tsgs_PCon1, tuPCon
    real(rprec), dimension(nynpy) :: arg1
    real(rprec), dimension(nz-1) :: arg2
    character(*), parameter :: fmt_5168 = "(1400(E14.5))"
    ! ---
    fr = (1._rprec/float(p_count))*float(c_count)

    if (jt .eq. c_count) then
        aPCon = 0._rprec; PCon2 = 0._rprec; asgs_PCon3 = 0._rprec; awPCon = 0._rprec;
        adPCondz = 0._rprec; var_PCon = 0._rprec; aKc_t = 0._rprec; aCs2Sc = 0._rprec
        asgs_PCon1 = 0._rprec; auPCon = 0._rprec
    end if

    call calc_hor_avg2(PCon(:,:,1:nz,npcon), arg2) 

    do k = 1, nz - 1
    do i = 1, nxt
        tPCon = 0._rprec; tPCon2 = 0._rprec; tsgs_PCon3 = 0._rprec; twPCon = 0._rprec
        tdPCondz = 0._rprec; tKc_t = 0._rprec; tCs2Sc = 0._rprec
        tsgs_PCon1 = 0._rprec; tuPCon = 0._rprec

        do j = 1, nynpy
            if ((k .eq. 1) .and. coordz == 0) then
                arg1 = 0._rprec
            else
                arg1 = (PCon(i, j, k, npcon) + PCon(i, j, k - 1, npcon))/2.
            end if
        end do

        call avg_y(PCon(i, :, k, npcon), tPCon)
        call avg_y(sgs_PCon3(i, :, k), tsgs_PCon3)
        call avg_y(dPCondz(i, :, k), tdPCondz)
        !call avg_y(Nu_t(i,:,k), tnu_t)
        call avg_y(Kc_t(i, :, k), tKc_t)
        call avg_y(Cs2Sc(i, :, k), tCs2Sc)
        call avg_y(PCon(i, :, k, npcon)*PCon(i, :, k, npcon), tPCon2)
        call avg_y(sgs_PCon1(i, :, k), tsgs_PCon1)
        call avg_y(u(i, :, k)*PCon(i, :, k, npcon), tuPCon)
        
        if (PCon_scheme == 1) then
            call avg_y(w(i, j, k)*arg1(:), twPCon)    
        else
            call avg_y(res_PCon3(i, :, k), twPCon)
        end if

        call calc_rms(PCon(1:nxt, 1:nynpy, k, npcon)-arg2(k), &
                      PCon(1:nxt, 1:nynpy, k, npcon)-arg2(k), tmp)
        var_PCon(i, k) = var_PCon(i, k) + (fr)*tmp

        aPCon(i, k) = aPCon(i, k) + (fr)*tPCon
        asgs_PCon3(i, k) = asgs_PCon3(i, k) + (fr)*tsgs_PCon3
        awPCon(i, k) = awPCon(i, k) + (fr)*twPCon
        adPCondz(i, k) = adPCondz(i, k) + (fr)*tdPcondz
        ! anu_t(i,k)=anu_t(i,k)+(fr)*tnu_t
        aKc_t(i, k) = aKc_t(i, k) + (fr)*tKc_t
        aCs2Sc(i, k) = aCs2Sc(i, k) + (fr)*tCs2Sc
        PCon2(i, k) = PCon2(i, k) + (fr)*tPCon2
        asgs_PCon1(i, k) = asgs_PCon1(i, k) + (fr)*tsgs_PCon1
        auPCon(i, k) = auPCon(i, k) + (fr)*tuPCon
    end do
    end do

    if (mod(jt, p_count) == 0) then
        if (average_dim_num == 1) then
            call collocate_MPI_averages_N(aPCon, avg_out, 35, 'PCon')
            call collocate_MPI_averages_N(PCon2, avg_out, 36, 'PCon2')
            call collocate_MPI_averages_N(asgs_PCon3, avg_out, 37, 'sgs_PCon3')
            call collocate_MPI_averages_N(awPCon, avg_out, 38, 'wPCon')
            call collocate_MPI_averages_N(adPCondz, avg_out, 39, 'dPCondz')
            call collocate_MPI_averages_N(asgs_PCon1, avg_out, 40, 'sgs_PCon1')
            call collocate_MPI_averages_N(auPCon, avg_out, 41, 'uPCon')
            !call collocate_MPI_averages_N(t3,avg_out,46,'PCon3')
            call collocate_MPI_averages_N(var_PCon, avg_out, 47, 'var_PCon');
            call collocate_MPI_averages_N(aKc_t, avg_out, 48, 'Kc_sgs');
            call collocate_MPI_averages_N(aCs2Sc, avg_out, 49, 'Cs2Sc_sgs');
        else if (average_dim_num == 2) then
            ax_1d = sum(aPCon, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 35, 'PCon')
            ax_1d = sum(PCon2, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 36, 'PCon2')
            ax_1d = sum(asgs_PCon3, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 37, 'sgs_PCon3')
            ax_1d = sum(awPCon, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 38, 'wPCon')
            ax_1d = sum(adPCondz, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 39, 'dPCondz')
            ax_1d = sum(asgs_PCon1, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 40, 'sgs_PCon1')
            ax_1d = sum(auPCon, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 41, 'uPCon')
            !ax_1d=sum(t3,dim=1)/nxt
            !call collocate_MPI_averages_N(ax_1d,avg_out_1d,46,'PCon3')
            ax_1d = sum(var_PCon, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 47, 'var_PCon');
            ax_1d = sum(aKc_t, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 48, 'Kc_sgs');
            ax_1d = sum(aCs2Sc, dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d, avg_out_1d, 49, 'Cs2Sc_sgs')
        end if

        aPCon = 0._rprec; PCon2 = 0._rprec; asgs_PCon3 = 0._rprec; awPCon = 0._rprec
        adPCondz = 0._rprec; var_PCon = 0._rprec; aKc_t = 0._rprec; aCs2Sc = 0._rprec
        asgs_PCon1 = 0._rprec; auPCon = 0._rprec;
    end if
    ! ---
    end subroutine pollen_slice
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
    subroutine collocate_MPI_averages_N(avg_var_proc, avg_var_tot_domain, file_ind, filename_str)
    !-----------------------------------------------------------------------
    !--The following subroutine does the collocation of the MPI arrays for
    ! averaging in avgslice and scalar_slice (in scalars_module2.f90)
    !-----------------------------------------------------------------------
    use param
    implicit none
    ! ---
    real(kind=rprec), dimension(dim1_size, dim2_size)     :: avg_var_proc
    real(kind=rprec), dimension(dim1_global, dim2_global) :: avg_var_tot_domain
    integer, intent(in) :: file_ind
    character(*), intent(in) :: filename_str
    
    integer :: ind1, ind2, jx
    character(len=256) :: local_filename
    character(*), parameter :: fmt_5168 = "(1400(E14.5))"

    local_filename = path//'result/aver_'//trim(filename_str)//'.out'
    avg_var_tot_domain = 0._rprec
   
    if ( average_dim_num .eq. 1 ) then
        do jx = 1, dim1_size
            call collect_data_z(avg_var_proc(jx, 1:nz-1), avg_var_tot_domain(jx, 1:nzt-1))    
        end do
    else if ( average_dim_num .eq. 2 ) then 
        call collect_data_z(avg_var_proc(1:nz-1,1), avg_var_tot_domain(1:nzt-1,1))
    end if

    if ( rank == 0 ) then
        open (file_ind, file=trim(local_filename), status="unknown", position="append")

        if (average_dim_num .eq. 1) then
            do ind2 = 1, nzt - 1
            do ind1 = 1, nxt
                if (abs(avg_var_tot_domain(ind1, ind2)) .lt. TINYS) then
                    avg_var_tot_domain(ind1, ind2) = 0._rprec
                end if
            end do
            end do
      
            do ind2 = 1, nzt - 1
                write (file_ind, fmt_5168) nums*dt, (avg_var_tot_domain(ind1, ind2), ind1=1, nxt)
            end do
        else if (average_dim_num .eq. 2) then
            do ind1 = 1, nzt - 1
                if (abs(avg_var_tot_domain(ind1, 1)) .lt. TINYS) then
                    avg_var_tot_domain(ind1, 1) = 0._rprec
                end if
            end do
            write (file_ind, fmt_5168) nums*dt, (avg_var_tot_domain(ind1, 1), ind1=1, nzt - 1)
        end if
    
        close (file_ind)
    end if
    ! ---
    end subroutine collocate_MPI_averages_N
! ---
end module scalars_module2



