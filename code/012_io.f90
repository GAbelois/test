!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------     
module io
!-----------------------------------------------------------------------
! Modified to fit namelist feature, use allocatable attribution and
! >> remove the parameter and save attribution (Bicheng Chen 06/15/2016)
!-----------------------------------------------------------------------     
use types, only:rprec
use param
use intermediate, only:average_dim_select_flag, dim1_size, dim2_size,&
                       dim1_global, dim2_global
implicit none
save
! ---
integer,parameter :: num_hour_out = 1
! TOMAS CHOR included base in param.nml
integer :: nwrite
!integer,parameter::base=200,nwrite=base

!!!!  io_spec=.true. output plan-averaged spectrum
logical,parameter :: io_spec=.true., output_fields_3d_flag=.false.
integer,parameter :: spec_write_freqz=3000 !fields_3d_write_freqz=p_count*6
integer,parameter :: spec_write_start=1,spec_write_end=1000000
!integer,parameter::spec_write_start=1,spec_write_end=24*base
!! --------------------------------------------------------------------
!! The following block defines parameters for instantaneous slice output
!! inst_slice_freqz controls the output frequency
!! The 5 slices outputted every inst_slice_freqz (i.e. u,v,w,T,Cs in this order) ...
!! ... are saved in the 3rd dimension of inst_slice_array for inst_array_lim times ...
!! ... following which this array is outputted as a binary file and the process
!! ... starts over
!! Therefore, each file will contain inst_array_lim x-z slices of 5 variables
!! This information is important for post-processing of these binary files
!! --------------------------------------------------------------------

logical, parameter :: inst_slice_flag = .false.
integer, parameter :: num_vars = 4                  ! works currently only for u,v,w,T due to the size difference in Cs
!integer, parameter :: slice_inst = (nzt-1)/2    ! sets the value of the y node for the chosen x-z inst slice
!integer, parameter :: inst_slice_freqz = 5         ! 
!integer, parameter :: inst_array_lim=200
!real(kind=rprec), dimension(nxt,nz+1,num_vars*inst_array_lim) :: inst_slice_array
!integer:: inst_slice_counter

logical, parameter :: cumulative_time = .true.
character(*), parameter :: fcumulative_time = '../readin/total_time.dat'

integer, parameter :: n_avg_stats = 100             ! interval for updates in avg_stats
character(*), parameter :: end_hdr_avg = '# end header'

!! --------------------------------------------------------------------
!! The following block defines parameters for use in avgslice and scalar_slice
!! --------------------------------------------------------------------
!integer,parameter :: average_dim_select_flag=1-(average_dim_num/2)
!! The variable average_dim_select_flag generates the following values based
!! on the value of average_dim_num in param.f90 :-
!! a) average_dim_num = 2 :
!! 	average_dim_select_flag = 0 ==> average the 3d array over x and y and output the z profile
!! b) average_dim_num = 1 :
!! 	average_dim_select_flag = 1 ==> average the 3d array over y and output the (x,z) profile
!integer, parameter :: dim1_size=average_dim_select_flag*(nxt-nz+1)+nz-1
!integer, parameter :: dim2_size=average_dim_select_flag*(nz-2)+1
!integer, parameter :: dim1_global=average_dim_select_flag*(nxt-nzt+1)+nzt-1
!integer, parameter :: dim2_global=average_dim_select_flag*(nzt-2)+1
!! --------------------------------------------------------------------
!! --------------------------------------------------------------------

! ---  time_spec>0 output time series spectrum (need additional calcu.)
integer, parameter :: time_spec = 0
integer :: n_obs, jt_total_init
integer, allocatable :: obs_pt(:,:)

! ---  io_mean=.true. output small domain time-averaged velocity
logical, parameter :: io_mean = .false.
integer :: jx_pls, jx_ple, width
integer :: jy_pls, jy_ple
real(rprec), dimension(:, :, :), allocatable :: mean_u, mean_v, mean_w,     &
                                                mean_u2, mean_v2, mean_w2

! --- mpi io
integer, dimension(3) :: gsize, lsize, start   ! dimensions of global and local variables 

!!!!  io_lambda2
!logical,parameter::io_lambda2=.false.
!real(rprec), dimension(:, :, :), allocatable :: lam2

! ---
contains
! --------------------------------------------------------------------
! file number 1-10 are used for temporary use
! 11-19 are basic output
! 20-40 are avgslice
! use >50 for debugging
! --------------------------------------------------------------------
!     subroutine openfiles()
! !-----------------------------------------------------------------------
! !   Open output files
! !-----------------------------------------------------------------------
!     use param, only:path
!     implicit none
!     ! --- to hold file names
!     !character(64) :: temp
!     !character(64) :: fCS1plan, fCS2plan, fCS4plan, fVISCplan,  &
!     !                 fDISSplan, fCS1Vplan, fCS2Vplan, fCS4Vplan

!     integer :: i
!     logical :: exst
!     ! ---
!     if (cumulative_time) then
!         inquire (file=fcumulative_time, exist=exst)
!         if (exst) then
!             open(1, file=fcumulative_time)
!             read(1, *) jt_total
!             jt_total_init = jt_total
!             close (1)
!         else  !--assume this is the first run on cumulative time
!             !write(*, *) 'file ', fcumulative_time, ' not found'
!             !write(*, *) 'assuming jt_total = 0'
!             jt_total = 0
!             jt_total_init = jt_total
!         end if
!     end if

!     if ( time_spec .gt. 0 ) then
!         open(15,file=path//'result/velspec.out',form='unformatted',position='append')
!         if(jt_total .eq. 0) rewind(15)
!     end if

!     if ( io_mean ) then
!         open(51,file=path//'result/mean_u.out',form='unformatted',position='append')
!         if ( jt_total .eq. 0 ) then
!             rewind(51)
!             write(51,*) jx_pls, jx_ple, jy_pls, jy_ple
!         end if
!     end if

! !    fCS1plan  = path // 'result/CS1plan.out'
! !    fCS2plan  = path // 'result/CS2plan.out'
! !    fCS4plan  = path // 'result/CS4plan.out'
! !    fVISCplan = path // 'result/VISCplan.out'
! !    fDISSplan = path // 'result/DISSplan.out'
! !    fCS1Vplan = path // 'result/CS1Vplan.out'
! !    fCS2Vplan = path // 'result/CS2Vplan.out'
! !    fCS4Vplan = path // 'result/CS4Vplan.out'
! !
! !    !$if ($MPI)
! !    ! --- append coordinate identifiers
! !    write(temp, '(".c",i0)') coordz
! !    fCS1plan  = trim (fCS1plan) // temp
! !    fCS2plan  = trim (fCS2plan) // temp
! !    fCS4plan  = trim (fCS4plan) // temp
! !    fVISCplan = trim (fVISCplan) // temp
! !    fDISSplan = trim (fDISSplan) // temp
! !    fCS1Vplan = trim (fCS1Vplan) // temp
! !    fCS2Vplan = trim (fCS2Vplan) // temp
! !    fCS4Vplan = trim (fCS4Vplan) // temp
! !    !$endif
! !
! !!open (90, file=fCS1plan, form='unformatted')
! !!open (91, file=fCS2plan, form='unformatted')
! !!open (92, file=fCS4plan, form='unformatted')
! !!open (93, file=fVISCplan, form='unformatted')
! !!open (94, file=fDISSplan, form='unformatted')
! !!open (95, file=fCS1Vplan, form='unformatted')
! !!open (96, file=fCS2Vplan, form='unformatted')
! !!open (97, file=fCS4Vplan, form='unformatted')

!     if ( time_spec .gt. 0 ) then
!         open(1, file=path//'obs.pt')
!         read(1,*) n_obs
!         allocate(obs_pt(1:2,n_obs))
!         do i = 1, n_obs
!             read(1,*) obs_pt(1:2,i)
!         end do
!         close(1)
!     end if
!     ! ---
!     end subroutine openfiles
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine screen_display()
!-----------------------------------------------------------------------
!   echo the simulation parameters to make sure they are 
!   inputed correctly
!-----------------------------------------------------------------------
    use param
    implicit none
    ! --- Print computational environment on screen
    if ( rank == 0 ) then
        open(1,file='../output/screen.out')
        write(1,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	    write(1,*) "++--------------------------------------------------++"
	    write(1,*) "++          OCEANIC CANOPY FLOW WITH LES            ++"
	    write(1,*) "++--------------------------------------------------++"
	    write(1,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	    write(1,*) ""
	    write(1,*) ""
        
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "++                 LES PARAMETERS                   ++"
        write(1,*) "++--------------------------------------------------++"
        if ( model .eq. 1 ) then  
            write(1,*) "    Smagrinsky model   " 
            write(1,*) "    Cs                =", cs
        else if ( model .eq. 2 ) then
            write(1,*) "    Dynamic smagrinsky model   "    
        else if ( model .eq. 3 ) then
            write(1,*) "    Scale-dependent dynamic smagrinsky model   " 
        else if ( model .eq. 4 ) then
            write(1,*) "    Lagrangian scale similar model   " 
        else if ( model .eq. 5 ) then
            write(1,*) "    Lagrangian scale-depandent model   "  
        end if
        write(1,*) "++--------------------------------------------------++"
        write(1,*) ""
        
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "++                BACIC PARAMETERS                  ++"
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "      nxt = ", nxt, "  nyt = ", nyt, "  NZ = ", nzt
        write(1,*) "++--------------------------------------------------++"
        write(1,*) "molecular viscosity (dimensional) = ", nu_molec
        write(1,*) "Number of timesteps               = ", nsteps
        write(1,*) "Time step size                    = ", dt
        write(1,*) "Number of processors              = ", npy*npz
        write(1,*) 'sampling stats every ', c_count, ' timesteps'
        write(1,*) 'writing stats every ', p_count, ' timesteps'
        write(1,*) "++--------------------------------------------------++"
        write(1,*) ""
        ! --- 
        close(1)
    end if
    ! ---
    end subroutine screen_display
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------     
    subroutine avg_stats ()
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
    use param
    use sim_param, only : u, v, w, txz, theta
    use sgsmodule
    use fft, only : kx_2d
    use intermediate, only: ubar_avg, vbar_avg, thetabar_avg, Cs2bar_avg,       &
                            Nutbar_avg, ubar_tot, vbar_tot, thetabar_tot,       &
                            Cs2bar_tot, Nutbar_tot, upr2bar_avg, stressbar_avg, &
                            upr2bar_tot, stressbar_tot, Eozbar_avg, Eozbar_tot
    implicit none
!--choose naming convention that does not conflict with qpost
    character(*), parameter :: fubar_avg      = path//'result/ubar-avg_stats.dat'
    character(*), parameter :: fvbar_avg      = path//'result/vbar-avg_stats.dat'
    character(*), parameter :: fupr2bar_avg   = path//'result/upr2bar-avg_stats.dat'
    character(*), parameter :: fstressbar_avg = path//'result/stressbar-avg_stats.dat'
    character(*), parameter :: fEozbar_avg    = path//'result/Eozbar-avg_stats.dat'
    character(*), parameter :: fthetabar_avg  = path//'result/thetabar-avg_stats.dat'
    character(*), parameter :: fCs2bar_avg    = path//'result/Cs2bar-avg_stats.dat'
    character(*), parameter :: fNutbar_avg    = path//'result/Nutbar-avg_stats.dat'

    integer, parameter :: hdr_len = 256
    character(hdr_len) :: Eozbar_hdr
  
    integer, save :: n_ubar_avg
    integer, save :: n_vbar_avg
    integer, save :: n_upr2bar_avg
    integer, save :: n_stressbar_avg
    integer, save :: n_Eozbar_avg
    integer, save :: n_thetabar_avg
    integer, save :: n_Cs2bar_avg
    integer, save :: n_Nutbar_avg
    integer :: jz

    logical, save :: init = .false.

    real(rprec), dimension(0:nz-1) :: ubar, vbar, wbar, thetabar, Nutbar
    real(rprec), dimension(nz-1) :: Cs2bar
    real(rprec), dimension(nxt, nynpy) :: upr, vpr, wpr
    real(rprec), dimension(3, nz-1) :: upr2bar, stressbar
    real(rprec) :: Eozbar(nxt/2, nz-1)
    real(rprec) :: z
    
    real(rprec) :: zu(1, nzt-1)
    real(rprec) :: kz_z(2, nxt/2)
    
    ! ---
    !--check whether or not to actually do anything
    !--motivation for doing this way is that it cleans up interface in main
    if (modulo(jt, n_avg_stats) == 0) then
        if ( rank == 0 ) then
            if (.not. init) then  !--initialization
                call init_avg (fubar_avg, 1, ubar_avg, n_ubar_avg)
                call init_avg (fvbar_avg, 1, vbar_avg, n_vbar_avg)
                call init_avg (fupr2bar_avg, 1, upr2bar_avg, n_upr2bar_avg)
                call init_avg (fstressbar_avg, 1, stressbar_avg, n_stressbar_avg)
      
                do jz = 1, nz-2
                    call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, jz), n_Eozbar_avg,  &
                                   leaveopn='yes')
                end do
                call init_avg (fEozbar_avg, 2, Eozbar_avg(:, :, nz-1), n_Eozbar_avg)
                call init_avg (fthetabar_avg, 1, thetabar_avg, n_thetabar_avg)
                call init_avg (fCs2bar_avg, 1, Cs2bar_avg, n_Cs2bar_avg)
                call init_avg (fNutbar_avg, 1, Nutbar_avg, n_Nutbar_avg)

                init = .true.
            end if
        end if
        
        call calc_hor_avg(u, ubar)
        call calc_hor_avg(v, vbar)
        call calc_hor_avg(w, wbar)
        if (theta_flag) call calc_hor_avg(theta, thetabar)
        call calc_hor_avg2(Cs_opt2, Cs2bar)
        call calc_hor_avg(Nu_t, Nutbar)
        call calc_hor_avg(txz, stressbar(2, :))

        do jz = 1, nz - 1
            if ( coordz == 0 .and. jz == 1 ) then
                upr = 0._rprec
                vpr = 0._rprec
                wpr = 0._rprec
            else
            !--see qpost for u/w-node interpolation
            !--convention will be to put up, vp, wp on w-nodes
                upr = 0.5_rprec * (u(1:nxt, 1:nynpy, jz)   - ubar(jz) +  &
                                   u(1:nxt, 1:nynpy, jz-1) - ubar(jz-1))
                vpr = 0.5_rprec * (v(1:nxt, 1:nynpy, jz)   - vbar(jz) +  &
                                   v(1:nxt, 1:nynpy, jz-1) - vbar(jz-1))
                wpr = w(1:nxt, 1:nynpy, jz) - wbar(jz)
            end if

            call calc_rms(upr, upr, upr2bar(1, jz))
            call calc_rms(vpr, vpr, upr2bar(2, jz))
            call calc_rms(wpr, wpr, upr2bar(3, jz))
            call calc_rms(upr, wpr, stressbar(1, jz))

            stressbar(3, jz) = sum (stressbar(1:2, jz))

            !--energy spectra
            call spectrum (u(:, :, jz), Eozbar(:, jz))  !--not /z yet
            z = (jz - 0.5_rprec) * dz       ! cyan should insert coordz or not
            Eozbar(:, jz) = Eozbar(:, jz) / z
        end do
        
        ! --- collect current stats into nzt sized arrays
        call collect_data_z(ubar(1:nz-1), ubar_tot(1, :))
        call collect_data_z(vbar(1:nz-1), vbar_tot(1, :))
        call collect_data_z(upr2bar(1, 1:nz-1), upr2bar_tot(1, :))
        call collect_data_z(upr2bar(2, 1:nz-1), upr2bar_tot(2, :))
        call collect_data_z(upr2bar(3, 1:nz-1), upr2bar_tot(3, :))
        call collect_data_z(stressbar(1,1:nz-1), stressbar_tot(1, :))
        call collect_data_z(stressbar(2,1:nz-1), stressbar_tot(2, :))
        call collect_data_z(stressbar(3,1:nz-1), stressbar_tot(3, :))
        call collect_data_z(Eozbar(:,1:nz-1), Eozbar_tot(1, :, :))
        if (theta_flag) then
            call collect_data_z(thetabar(1:nz-1), thetabar_tot(1, :))
        end if
        call collect_data_z(Cs2bar(1:nz-1), Cs2bar_tot(1, :))
        call collect_data_z(Nutbar(1:nz-1), Nutbar_tot(1, :))
        
        if ( rank == 0 ) then
            ! --- calculation of cumulative average stats
            ubar_avg = (n_ubar_avg * ubar_avg + ubar_tot) / (n_ubar_avg + 1)
            n_ubar_avg = n_ubar_avg + 1
        
            vbar_avg = (n_vbar_avg * vbar_avg + vbar_tot) / (n_vbar_avg + 1)
            n_vbar_avg = n_vbar_avg + 1
        
            upr2bar_avg = (n_upr2bar_avg * upr2bar_avg + upr2bar_tot) /  &
               (n_upr2bar_avg + 1)
            n_upr2bar_avg = n_upr2bar_avg + 1
        
            stressbar_avg = (n_stressbar_avg * stressbar_avg + stressbar_tot) /  &
               (n_stressbar_avg + 1)
            n_stressbar_avg = n_stressbar_avg + 1
        
            Eozbar_avg = (n_Eozbar_avg * Eozbar_avg + Eozbar_tot) / (n_Eozbar_avg + 1)
            n_Eozbar_avg = n_Eozbar_avg + 1
        
            thetabar_avg = (n_thetabar_avg * thetabar_avg + thetabar_tot) / (n_thetabar_avg + 1)
            n_thetabar_avg = n_thetabar_avg + 1
        
            Cs2bar_avg = (n_Cs2bar_avg * Cs2bar_avg + Cs2bar_tot) / (n_Cs2bar_avg + 1)
            n_Cs2bar_avg = n_Cs2bar_avg + 1
        
            Nutbar_avg = (n_Nutbar_avg * Nutbar_avg + Nutbar_tot) / (n_Nutbar_avg + 1)
            n_Nutbar_avg = n_Nutbar_avg + 1
        
            !--prepare list of z-coordinates
            forall (jz=1:nzt-1) zu(1, jz) = (jz - 0.5_rprec) * dz
        
            !--prepare  header, optional
        
            !--write out to file
            call write_avg (fubar_avg, n_ubar_avg, zu, ubar_avg)
            call write_avg (fvbar_avg, n_vbar_avg, zu, vbar_avg)
            call write_avg (fupr2bar_avg, n_upr2bar_avg, zu, upr2bar_avg)
            call write_avg (fstressbar_avg, n_stressbar_avg, zu, stressbar_avg)
            call write_avg (fthetabar_avg, n_thetabar_avg, zu, thetabar_avg)
            call write_avg (fCs2bar_avg, n_Cs2bar_avg, zu, Cs2bar_avg)
            call write_avg (fNutbar_avg, n_Nutbar_avg, zu, Nutbar_avg)
        
            !--this is a bit awkward: maybe just have another routine to do it right
            Eozbar_hdr = 'zone' !--this is for tecplot...
            kz_z(1, :) = kx_2d(1:nxt/2, 1) * zu(1, 1)
            kz_z(2, :) = zu(1, 1)
            call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z,   &
                            Eozbar_avg(:, :, 1), hdr=Eozbar_hdr)
        
            do jz = 2, nzt - 1
                kz_z(1, :) = kx_2d(1:nxt/2, 1) * zu(1, jz)
                kz_z(2, :) = zu(1, jz)
        
                call write_avg (fEozbar_avg, n_Eozbar_avg, kz_z,   &
                                Eozbar_avg(:, :, jz), hdr=Eozbar_hdr, position='append')
            end do
        
        end if
    end if
    ! ---
    end subroutine avg_stats
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------     
    subroutine init_avg (file_avg, n_ccol, a_avg, n_avg, leaveopn)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------   
    implicit none
    ! ---
    character(*), intent(in) :: file_avg
    integer, intent(in) :: n_ccol  !--num. coordz columns: x, y, etc.
    real(rprec), intent(out) :: a_avg(:, :)
    integer, intent(out) :: n_avg
    character (*), optional, intent (in) :: leaveopn
    ! ---
    character(128) :: buff
    logical :: exst, opn
    integer :: j
    real(rprec) :: z(n_ccol)
! ---
    inquire (file=file_avg, exist=exst, opened=opn)
    if (exst) then
        if (.not. opn) then
            open(1, file=file_avg)
            read(1, '(a)') buff

            if (buff(1:1) == '#') then
                read (buff(2:), *) n_avg
            else
                write (*, *) 'avg_stats: error'
                write (*, *) trim (file_avg), ' does not have expected format on line 1'
                stop  !--need to replace all stops with nice mpi exits
            end if
        end if

    !--skip data header lines here
        do
            read (1, '(a)') buff
            if (trim (buff) == trim (end_hdr_avg)) exit
        end do

        do j = 1, size (a_avg, 2)
            read (1, *) z, a_avg(:, j)  !--z is just placeholder here
        end do

        if (present (leaveopn)) then
            if (leaveopn /= 'yes') close (1)  !--case sensitive here
        else
            close (1)
        end if
    else
        n_avg = 0
        a_avg = 0._rprec
    end if
    ! ---
    end subroutine init_avg
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------     
    subroutine write_avg (file_avg, n_avg, x, a_avg, hdr, position)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    implicit none
    ! ---
    character(*), intent(in) :: file_avg
    integer, intent(in) :: n_avg
    real(rprec), intent(in) :: x(:, :)  !--coordz columns for x, y, etc
    real(rprec), intent(in) :: a_avg(:, :)

    character(*), optional, intent(in) :: hdr
    character(*), optional, intent(in) :: position

    character(64) :: r_fmt, fmt
    character(32) :: posn

    integer :: j

! --- check sizes compatible
    if (size (x, 2) /= size (a_avg, 2)) then
        write (*, *) 'write_avg: error with sizes of x, a_avg'
        stop
    end if

    if (present (position)) then
        posn = position
    else
        posn = 'rewind'
    end if

    open (1, file=file_avg, position=posn)
    if (trim (posn) /= 'append') then  !--case sensitive
        write (1, '(a,i0)') '# ', n_avg  !--not needed when appending
    end if

    if (present (hdr)) then
        !--write data header, if present
        write(1, '(a)') trim (hdr)
    end if

    !--write something to indicate end of header, always do this
    write(1, '(a)') end_hdr_avg

    !--define output format
    write(r_fmt, '(2(a,i0))') 'es', precision (1._rprec) + 7,  &
                               '.', precision (1._rprec)
    write(fmt, '(a,i0,3a)') '(', size (x, 1) + size (a_avg, 1),  &
                            '(1x,', trim (r_fmt), '))'
!--write to file
    do j = 1, size (a_avg, 2)
        write (1, fmt) x(:, j), a_avg(:, j)
    end do
    close (1)
    ! ---
    end subroutine write_avg
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine post_spec(jt_local)
    
    use sim_param,only:u,v,w,theta,pcon
    use param
  use fft

  implicit none
  real(kind=rprec),dimension(4,nxt/2,nz-1)::spectra_uvwT
  real(kind=rprec),dimension(4,nxt/2,nzt-1)::spectra_uvwT_tot
  integer,intent(in)::jt_local
  integer::k,jz!,z

  real(kind=rprec),dimension(nz-1)::z   ! Ying changed z to array (09/30/2010)
  ! Chamecki changed z from integer to real (08/04/2006)
  real(kind=rprec),dimension(nzt-1)::z_tot

  character(len=64)::fname1
  !$if ($MPI)
  !$define $lbz 0
  integer :: recvcounts(npz)
  integer :: displs(npz)
  !$else
  !$define $lbz 1
  !$endif

  !$if ($MPI)
  recvcounts = size (z)
  displs = coord_of_rank * recvcounts
  do jz=1,nz-1
    z(jz) = (coordz*(nz-1)+jz-0.5_rprec)*dz*z_i
  enddo
  call mpi_gatherv (z(1), size (z), MPI_RPREC,&
     z_tot(1), recvcounts, displs,       &
     MPI_RPREC, 0, comm, ierr)
  !$else
  !do jz=1,nz-1
  !  z(jz)=(jz-0.5_rprec)*dz*z_i
  !  z_tot(jz)=z(jz)
  !enddo
  !$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coordz == 0)) then
    write(fname1,'(a,a)') path//'output/spec_x','.dat'
    open(82,file=fname1,form='formatted')
    do jz = 1,nzt-1
      write(82,*) (real(kx_2d(k,1)/z_i*z_tot(jz)),k=1,nxt/2)
    enddo
    close(82)
  endif

!write(fname1,'(a,a)') path//'output/spec_x','.dat'
!open(82,file=fname1,form='formatted')
  do jz=1,nz-1
!do jz=$lbz,nz-1
!   z=(jz-0.5_rprec)*dz*z_i
!   write(82,*) (real(2*pi/lx_tot*k/z_i*z),k=1,nxt/2-1)
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! kx_2d=(2pi/lx_tot)*(0:nxt/2) => lx_tot is already non-dim by z_i
!! => k_x=[(2pi/lx_tot)*(0:nxt/2)/z_i] and
!! kz is given by = k_x*z => kz=[(2pi/lx_tot)*(0:nxt/2)*([1:nz]-0.5)*dz]
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   write(82,*) (real(kx_2d(k,1)/z_i*z),k=1,nxt/2)
! Calculate spectrum now !! for each slice in MPI case !
!   call spectrum(u(:, :, jz), spectra_u(:,jz))
!   call spectrum(v(:, :, jz), spectra_v(:,jz))
!   call spectrum(w(:, :, jz), spectra_w(:,jz))
!   call spectrum(theta(:,:,jz),spectra_theta(:,jz))
    call spectrum(u(:, :, jz), spectra_uvwT(1,:,jz))
    call spectrum(v(:, :, jz), spectra_uvwT(2,:,jz))
    call spectrum(w(:, :, jz), spectra_uvwT(3,:,jz))
    call spectrum(theta(:, :, jz), spectra_uvwT(4,:,jz))

    ! Replace temperature spectra by Pollen concentration
    ! Chamecki - 08/10/2006
    IF (PCon_FLAG) call spectrum(PCon(:, :, jz,npcon), spectra_uvwT(4,:,jz))

  enddo
!   close(82)
!   print *,'spectra_sample_U',spectra_u(:,4)
!   print *,'spectra_sample_U2',spectra_uvwT(1,:,4)
!   print *,'spectra_sample_V',spectra_v(:,4)
!   print *,'spectra_sample_V2',spectra_uvwT(2,:,4)
!   print *,'spectra_sample_W',spectra_w(:,4)
!   print *,'spectra_sample_W2',spectra_uvwT(3,:,4)
!   print *,'spectra_sample_T',spectra_theta(:,4)
!   print *,'spectra_sample_T2',spectra_uvwT(4,:,4)
  !$if ($MPI)
  recvcounts = size (spectra_uvwT)
  displs = coord_of_rank * recvcounts
  call mpi_gatherv (spectra_uvwT(1, 1,1), size (spectra_uvwT), MPI_RPREC,&
     spectra_uvwT_tot(1, 1, 1), recvcounts, displs,       &
     MPI_RPREC, 0, comm, ierr)
  !$else
  !spectra_uvwT_tot=spectra_uvwT
  !$endif

  if((.not. USE_MPI) .or. (USE_MPI .and. coordz == 0)) then
    write(fname1,'(A,i6.6,A)')path//'output/spec_uvwT_',jt_local,'.bin'
    open(83,file=fname1,form='unformatted')

!write(fname1,'(a,i6.6,a)')path//'output/spec_u',jt_local,'.dat'
!open(1,file=fname1,form='formatted')
!write(fname2,'(A,i6.6,A)')path//'output/spec_v',jt_local,'.dat'
!open(2,file=fname2,form='formatted')
!write(fname3,'(A,i6.6,A)')path//'output/spec_w',jt_local,'.dat'
!open(3,file=fname3,form='formatted')
!write(fname4,'(A,i6.6,A)')path//'output/spec_t',jt_local,'.dat'
!open(4,file=fname4,form='formatted')

    write(83) real(spectra_uvwT_tot(:,1:nxt/2,:))
    close(83)
!do jz=1,nz
!   write(1,*)real(spectra_u(2:nxt/2,jz))
!   write(2,*)real(spectra_v(2:nxt/2,jz))
!   write(3,*)real(spectra_w(2:nxt/2,jz))
!   write(4,*)real(spectra_theta(2:nxt/2,jz))
!enddo
!close(1);close(2);close(3);close(4)
  end if
    ! ---
    end subroutine post_spec
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
    subroutine spectrum(u, spec)
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    use fft
    implicit none
    ! ---
    real(kind=rprec), dimension(nxt, nynpy), intent(in) :: u
    real(kind=rprec), dimension(nxt/2), intent(out) :: spec  !--assumes Nyquist is 0
    ! ---
    integer :: jy, k 
    real(kind=rprec), dimension(nxt) :: vel_r, vel_c
    real(kind=rprec), dimension(nxt/2) :: tmp 

    integer*8, save :: plan
    logical, save :: init = .false.

    if (.not. init) then
        call dfftw_plan_r2r_1d(plan, nxt, vel_r, vel_c, FFTW_R2HC, FFTW_MEASURE)
        init = .true.
    end if

    ! initialize
    tmp(:) = 0._rprec
    do jy = 1, nynpy
        vel_r(:) = u(1:nxt,jy)/real(nxt,kind=rprec)
        ! check this normaliztion-part of forward
        ! call the fft
        call dfftw_execute_r2r(plan, vel_r, vel_c)
        ! compute magnitudes
        ! the 0.5 is the 1/2, all others are taken care of! (except maybe Nyquist)
        tmp(1) = tmp(1) + 0.5*vel_c(1)*vel_c(1)
    
        do k = 2, nxt/2
            tmp(k) = tmp(k) + vel_c(k)*vel_c(k) + vel_c(nxt+2-k)*vel_c(nxt+2-k)
            ! print *,'k,vel,spec',k,vel_c(k),tmp(k)
        end do
        !--assume Nyquist is 0
        !tmp(nxt/2+1)=tmp(nxt/2+1)+vel_c(nxt/2+1)*vel_c(nxt/2+1)
    end do
    tmp(:) = tmp(:)/real(nyt,kind=rprec) ! for average over nyt
    call mpi_allreduce(tmp, spec, nxt/2, mpi_rprec, mpi_sum, comm_ver, ierr)
    ! ---
    end subroutine spectrum  
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
!    subroutine timeseries_spec
!!-----------------------------------------------------------------------
!!   
!!-----------------------------------------------------------------------  
!    use sim_param,only:u,v,w,theta
!    implicit none
!    ! ---
!    integer :: jx, jy, jz, i
!    
!    if (mod(jt_total,time_spec)==0 .and. jt_total.gt.2000) then
!        jx = nxt/8
!        jy = nyt/2+1
!        jz = NZ/2
!!TSwrite(15)real(u(jx+nxt/24*2,jy,jz:jz+3)),real(u(jx+nxt/24*4,jy,jz:jz+3)),&
!!TS     real(u(jx+nxt/24*6,jy,jz:jz+3)),&
!!TS     real(v(jx+nxt/24*2,jy,jz:jz+3)),real(v(jx+nxt/24*4,jy,jz:jz+3)),&
!!TS     real(v(jx+nxt/24*6,jy,jz:jz+3)),&
!!TS     real(w(jx+nxt/24*2,jy,jz:jz+3)),real(w(jx+nxt/24*4,jy,jz:jz+3)),&
!!TS     real(w(jx+nxt/24*6,jy,jz:jz+3))
!    end if
!    ! ---
!    end subroutine timeseries_spec
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
!-----------------------------------------------------------------------
!   
!----------------------------------------------------------------------- 
    subroutine avgslice
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------  
    use sim_param,only:u,v,w,dudz,dvdz, txx, txz, tyy, tyz, tzz, p,dudt,dvdt,dwdt
    use param,only:path,p_count,c_count,jt
    use sgsmodule,only:Cs_opt2,Cs_Ssim,Beta_avg,Betaclip_avg
    use intermediate, only: ap, au, av, aw, p2, u2, v2, w2, auw, avw, acs,&
        adudz, advdz, aCs_Ssim, abeta_sgs, abetaclip_sgs, atxx, atxz, atyy,&
        atyz, atzz, u3, v3, w3, adudt, advdt, adwdt, ax_1d, avg_out_1d, avg_out
    !use scalars_module2,only:collocate_MPI_averages
    implicit none
    ! ---
    integer :: i, j, k
    real(rprec) :: tu1,tv1,tw1,ttxx,ttxz,ttyy,ttyz,ttzz,tdudz,tdvdz,    &
                   tu2,tv2,tw2,tp1,tp2,tuw,tvw,tCs,fr,                  &
                   tu3, tv3, tw3, tdudt, tdvdt, tdwdt,tCs_Ssim
    real(rprec), dimension(nynpy) :: arg1, arg2
    ! ---
    fr=(1._rprec/real(p_count,kind=rprec))*real(c_count,kind=rprec)
    do k = 1, nz-1
    do i = 1, nxt
        tu1=0._rprec;tv1=0._rprec;tw1=0._rprec;tp1=0._rprec
        ttxx=0._rprec;ttxz=0._rprec;ttyy=0._rprec;ttyz=0._rprec
        ttzz=0._rprec;tdudz=0._rprec;tdvdz=0._rprec;tu2=0._rprec
        tv2=0._rprec;tw2=0._rprec;tp2=0._rprec;tuw=0._rprec;tvw=0._rprec
        tCs=0._rprec;tCs_Ssim=0._rprec;tu3=0._rprec;tv3=0._rprec;tw3=0._rprec
        tdudt=0._rprec;tdvdt=0._rprec;tdwdt=0._rprec

        do j = 1, nynpy
            if ( (k .eq. 1) .and. coordz == 0 ) then
                arg1(j) = 0._rprec
                arg2(j) = 0._rprec
            else
                arg1(j) = (u(i,j,k)+u(i,j,k-1))/2.
                arg2(j) = (v(i,j,k)+v(i,j,k-1))/2.
            end if
        end do
            
        call avg_y(u(i,:,k), tu1)
        call avg_y(v(i,:,k), tv1)
        call avg_y(w(i,:,k), tw1)
        call avg_y(p(i,:,k), tp1)
        call avg_y(txx(i,:,k), ttxx)
        call avg_y(txz(i,:,k), ttxz)
        call avg_y(tyy(i,:,k), ttyy)
        call avg_y(tyz(i,:,k), ttyz)
        call avg_y(tzz(i,:,k), ttzz)
        call avg_y(dudz(i,:,k), tdudz)
        call avg_y(dvdz(i,:,k), tdvdz)
        call avg_y(u(i,:,k)*u(i,:,k), tu2)
        call avg_y(v(i,:,k)*v(i,:,k), tv2)
        call avg_y(w(i,:,k)*w(i,:,k), tw2)
        call avg_y(p(i,:,k)*p(i,:,k), tp2)
        call avg_y(sqrt(Cs_opt2(i,:,k)), tCs)
        call avg_y(sqrt(Cs_Ssim(i,:,k)), tCs_Ssim)
        call avg_y(w(i,:,k)*arg1(:), tuw)
        call avg_y(w(i,:,k)*arg2(:), tvw)
        call avg_y(u(i,:,k)*u(i,:,k)*u(i,:,k), tu3)
        call avg_y(v(i,:,k)*v(i,:,k)*v(i,:,k), tv3)
        call avg_y(w(i,:,k)*w(i,:,k)*w(i,:,k), tw3)
        call avg_y(dudt(i,:,k), tdudt)
        call avg_y(dvdt(i,:,k), tdvdt)
        call avg_y(dwdt(i,:,k), tdwdt)
      
        au(i,k)=au(i,k)+(fr)*tu1
        av(i,k)=av(i,k)+(fr)*tv1
        aw(i,k)=aw(i,k)+(fr)*tw1
        ap(i,k)=ap(i,k)+(fr)*tp1
        adudz(i,k)=adudz(i,k)+(fr)*tdudz
        advdz(i,k)=advdz(i,k)+(fr)*tdvdz
        u2(i,k)=u2(i,k)+(fr)*tu2
        v2(i,k)=v2(i,k)+(fr)*tv2
        w2(i,k)=w2(i,k)+(fr)*tw2
        atxx(i,k)=atxx(i,k)+(fr)*ttxx
        atxz(i,k)=atxz(i,k)+(fr)*ttxz
        atyy(i,k)=atyy(i,k)+(fr)*ttyy
        atyz(i,k)=atyz(i,k)+(fr)*ttyz
        atzz(i,k)=atzz(i,k)+(fr)*ttzz
        p2(i,k)=p2(i,k)+fr*tp2
        aCs(i,k)=aCs(i,k)+(fr)*tCs
        auw(i,k)=auw(i,k)+(fr)*tuw
        avw(i,k)=avw(i,k)+(fr)*tvw
        aCs_Ssim(i,k)=aCs_Ssim(i,k)+fr*tCs_Ssim
        abeta_sgs(i,k)=abeta_sgs(i,k)+fr*Beta_avg(k)
        abetaclip_sgs(i,k)=abetaclip_sgs(i,k)+fr*Betaclip_avg(k)
        u3(i,k)=u3(i,k)+(fr)*tu3
        v3(i,k)=v3(i,k)+(fr)*tv3
        w3(i,k)=w3(i,k)+(fr)*tw3
        adudt(i,k)=adudt(i,k)+(fr)*tdudt
        advdt(i,k)=advdt(i,k)+(fr)*tdvdt
        adwdt(i,k)=adwdt(i,k)+(fr)*tdwdt
    end do
    end do

    if (mod(jt,p_count)==0) then
        if (average_dim_num == 1) then
            call collocate_MPI_averages_N(au,avg_out,20,'u')
            call collocate_MPI_averages_N(av,avg_out,21,'v')
            call collocate_MPI_averages_N(aw,avg_out,22,'w')
            call collocate_MPI_averages_N(ap,avg_out,23,'p')
            call collocate_MPI_averages_N(u2,avg_out,24,'u2')
            call collocate_MPI_averages_N(v2,avg_out,25,'v2')
            call collocate_MPI_averages_N(w2,avg_out,26,'w2')
            call collocate_MPI_averages_N(p2,avg_out,32,'p2')
            call collocate_MPI_averages_N(atxx,avg_out,27,'txx')
            call collocate_MPI_averages_N(atxz,avg_out,28,'txz')
            call collocate_MPI_averages_N(atyy,avg_out,29,'tyy')
            call collocate_MPI_averages_N(atyz,avg_out,30,'tyz')
            call collocate_MPI_averages_N(atzz,avg_out,31,'tzz')
            call collocate_MPI_averages_N(auw,avg_out,33,'uw')
            call collocate_MPI_averages_N(avw,avg_out,34,'vw')
            call collocate_MPI_averages_N(aCs,avg_out,35,'Cs')
            call collocate_MPI_averages_N(adudz,avg_out,36,'dudz')
            call collocate_MPI_averages_N(advdz,avg_out,37,'dvdz')
            call collocate_MPI_averages_N(aCs_Ssim,avg_out,38,'Cs_Ssim')
            call collocate_MPI_averages_N(abeta_sgs,avg_out,39,'beta_sgs')
            call collocate_MPI_averages_N(abetaclip_sgs,avg_out,40,'betaclip_sgs');
            call collocate_MPI_averages_N(u3,avg_out,41,'u3')
            call collocate_MPI_averages_N(v3,avg_out,42,'v3')
            call collocate_MPI_averages_N(w3,avg_out,43,'w3')
            call collocate_MPI_averages_N(adudt,avg_out,44,'dudt')
            call collocate_MPI_averages_N(advdt,avg_out,45,'dvdt')
            call collocate_MPI_averages_N(adwdt,avg_out,46,'dwdt')
        else if (average_dim_num == 2) then
            ax_1d=sum(au,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,20,'u')
            ax_1d=sum(av,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,21,'v')
            ax_1d=sum(aw,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,22,'w')
            ax_1d=sum(ap,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,23,'p')
            ax_1d=sum(u2,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,24,'u2')
            ax_1d=sum(v2,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,25,'v2')
            ax_1d=sum(w2,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,26,'w2')
            ax_1d=sum(p2,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,32,'p2')
            ax_1d=sum(atxx,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,27,'txx')
            ax_1d=sum(atxz,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,28,'txz')
            ax_1d=sum(atyy,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,29,'tyy')
            ax_1d=sum(atyz,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,30,'tyz')
            ax_1d=sum(atzz,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,31,'tzz')
            ax_1d=sum(auw,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,33,'uw')
            ax_1d=sum(avw,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,34,'vw')
            ax_1d=sum(aCs,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,35,'Cs')
            ax_1d=sum(adudz,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,36,'dudz')
            ax_1d=sum(advdz,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,37,'dvdz')
            ax_1d=sum(aCs_Ssim,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,38,'Cs_Ssim')
            ax_1d=sum(abeta_sgs,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,39,'beta_sgs')
            ax_1d=sum(abetaclip_sgs,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,40,'betaclip_sgs');
            ax_1d=sum(u3,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,41,'u3')
            ax_1d=sum(v3,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,42,'v3')
            ax_1d=sum(w3,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,43,'w3')
            ax_1d=sum(adudt,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,44,'dudt')
            ax_1d=sum(advdt,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,45,'dvdt')
            ax_1d=sum(adwdt,dim=1)/nxt
            call collocate_MPI_averages_N(ax_1d,avg_out_1d,46,'dwdt')
        end if
    
        !VK Zero out the outputted averages !!
        au=0._rprec;av=0._rprec;aw=0._rprec;ap=0._rprec;u2=0._rprec;v2=0._rprec
        w2=0._rprec;atxx=0._rprec;atxz=0._rprec;atyy=0._rprec;atyz=0._rprec
        atzz=0._rprec;p2=0._rprec;auw=0._rprec;avw=0._rprec;aCs=0._rprec
        adudz=0._rprec;advdz=0._rprec;aCs_Ssim=0._rprec;abeta_sgs=0._rprec
        abetaclip_sgs=0._rprec;u3=0._rprec;v3=0._rprec;w3=0._rprec;
        adudt=0._rprec;advdt=0._rprec;adwdt=0._rprec;
    end if
    ! ---
    end subroutine avgslice
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
!     subroutine checkpoint ()
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
!     use param
!     use sim_param, only : u, v, w, theta, pcon
!     use scalars_module, only : deposition,Real_dep,P_surf_flux,P_surf_flux_dep
!     use bottombc, only: ustar_avg
!     use sgsmodule
!     implicit none
!     !
!     integer :: lun=1
!     !move from outputloop (Bicheng Chen 06/26/2015)
!     character (64) :: fname
!     integer :: jx, jy, jz, ipcon
!     integer, dimension(npz) :: sendcounts, recvcounts, displs
  
!     real (rprec), dimension(nxt, nyt, 1:nzt) :: u_tot
!     real (rprec), dimension(nxt, nyt, 1:nzt) :: v_tot
!     real (rprec), dimension(nxt, nyt, 1:nzt) :: w_tot
!     real (rprec), dimension(nxt, nyt, 1:nzt) :: theta_tot
!     real (rprec), dimension(nxt, nyt, 1:nzt, npcon) :: PCon_tot
!     real (rprec), dimension(nxt, nyt, 1:nzt) :: Cs_opt2_tot
!     real (rprec), dimension(nxt, nyt, 1:nzt) :: Nu_t_tot

!     character(*), parameter :: fmt_505 = "(12e12.4)"
!     ! ---
!     if (theta_flag .OR. PCon_FLAG) then
!         write (fname, '(a,i8.8,a)') path // 'result/vel_sc', jt_total, '.out'
!     else
!         write (fname, '(a,i8.8,a)') path // 'result/vel', jt_total, '.out'
!     end if

!     sendcounts = size (u(:,:,1:nz))
!     recvcounts = size (u(:,:,1:nz-1))
!     displs = coord_of_rank * recvcounts
!     call mpi_gatherv (u(1,1,1), sendcounts, MPI_RPREC,&
!                     u_tot(1,1,1), sendcounts, displs,       &
!                     MPI_RPREC, 0, comm, ierr)
!     call mpi_gatherv (v(1,1,1), sendcounts, MPI_RPREC,&
!                     v_tot(1,1,1), sendcounts, displs,       &
!                     MPI_RPREC, 0, comm, ierr)
!     call mpi_gatherv (w(1,1,1), sendcounts, MPI_RPREC,&
!                     w_tot(1,1,1), sendcounts, displs,       &
!                     MPI_RPREC, 0, comm, ierr)
!     call mpi_gatherv (Cs_opt2(1,1,1), sendcounts, MPI_RPREC,&
!                     Cs_opt2_tot(1,1,1), sendcounts, displs,       &
!                     MPI_RPREC, 0, comm, ierr)
!     call mpi_gatherv (Nu_t(1,1,1), sendcounts, MPI_RPREC,&
!                     Nu_t_tot(1,1,1), sendcounts, displs,       &
!                     MPI_RPREC, 0, comm, ierr)
!     if (theta_flag) then
!         call mpi_gatherv (theta(1,1,1), sendcounts, MPI_RPREC,&
!                         theta_tot(1,1,1), sendcounts, displs,       &
!                         MPI_RPREC, 0, comm, ierr)
!     end if
  
!     if (PCon_FLAG) then
!         do ipcon=1,npcon
!             call mpi_gatherv (PCon(1,1,1,ipcon), sendcounts, MPI_RPREC,&
!                             PCon_tot(1,1,1,ipcon), sendcounts, displs,       &
!                             MPI_RPREC, 0, comm, ierr)
!         end do
!     end if
  

!     if ( coordz == 0 ) then
!         open(lun,file=fname,form='unformatted')
!         if (theta_flag .AND. PCon_FLAG) then ! With scalars and pollen - Ying 11/01/2010
!             write (lun) u_tot(:, :, 1:nzt), v_tot(:, :, 1:nzt), w_tot(:, :, 1:nzt), &
!                         theta_tot(:,:,1:nzt), PCon_tot(:,:,1:nzt,1:npcon), &
!                         deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
!                         P_surf_flux(:,:), P_surf_flux_dep(:,:)
!         elseif (theta_flag) then !WITH SCALARS
!             write (lun) u_tot(:, :, 1:nzt), v_tot(:, :, 1:nzt), w_tot(:, :, 1:nzt), &
!                         theta_tot(:,:,1:nzt)
!         elseif (PCon_FLAG) then ! Pollen - Chamecki 08/21/2006
!             write (lun) u_tot(:, :, 1:nzt), v_tot(:, :, 1:nzt), w_tot(:, :, 1:nzt), &
!                         PCon_tot(:,:,1:nzt,1:npcon),   &
!                         deposition(:,:), Real_dep(:,:), ustar_avg(:,:), &
!                         P_surf_flux(:,:), P_surf_flux_dep(:,:)
! !!!!!            !DY Added by Di Yang for testing
! !!!!!            !  write(100,*) 'variables=x,y,z,PCon'
! !!!!!            open(1000000+jt_total)
! !!!!!            !  write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,Pcon2,Pcon3,Cs_opt2,Nu_t'
! !!!!!            write(1000000+jt_total,*) 'variables=x,y,z,u,v,w,Pcon1,PCon2,Cs_opt2,Nu_t'
! !!!!!            write(1000000+jt_total,*) 'zone t="',jt_total,'" i=',nxt,' j=',nyt,' k=',nzt,' f=point'
! !!!!!            do jz=1,nzt
! !!!!!            do jy=1,nyt
! !!!!!            do jx=1,nxt
! !!!!!                write(1000000+jt_total,fmt_505) (jx-1)*lx_tot*z_i/nxt,(jy-1)*ly_tot*z_i/nyt,(jz-1)*z_i/nzt, &
! !!!!!!                               PCon_tot(jx,jy,jz,1),PCon_tot(jx,jy,jz,2),PCon_tot(jx,jy,jz,3), &
! !!!!!                                PCon_tot(jx,jy,jz,1), PCon_tot(jx,jy,jz,2), &
! !!!!!                                Cs_opt2_tot(jx,jy,jz),Nu_t_tot(jx,jy,jz)
! !!!!!            enddo
! !!!!!            enddo
! !!!!!            enddo
! !!!!!            close(1000000+jt_total)
! !!!!!            !DY End here
!         else ! No SCALARS
!             write (lun) u_tot(:, :, 1:nzt), v_tot(:, :, 1:nzt), w_tot(:, :, 1:nzt)
!         end if
!     end if

!     ! move from outputloop (Bicheng Chen 06/26/2015)
!     if ( coordz == 0 ) then
!         close(lun)
!     end if
!     ! ---
!     end subroutine checkpoint
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
!     subroutine checkVSpl()
! !++++++++++++++++++++++++++++++++++++++++++++++++++
! !BC program to output velocity sample
!   use param
!   use sim_param, only : u, v, w, theta


!   implicit none

!   character(80) :: fname
!   integer :: iSpl, ixs, ixe, iys, iye, izs, ize

!   !$if ($MPI)
!   integer :: sendcounts(npz)
!   integer :: recvcounts(npz)
!   integer :: displs(npz)
!   !$endif

!   real (rprec), dimension(nxt, nyt, 1:nzt) :: u_tot
!   real (rprec), dimension(nxt, nyt, 1:nzt) :: v_tot
!   real (rprec), dimension(nxt, nyt, 1:nzt) :: w_tot
!   real (rprec), dimension(nxt, nyt, 1:nzt) :: theta_tot

! !---------------------------------------------------------------------
! !-BC gather the related variables
!   !$if ($MPI)
!   sendcounts = size (u(:,:,1:nz))
!   recvcounts = size (u(:,:,1:nz-1))
!   displs = coord_of_rank * recvcounts
!   call mpi_gatherv (u(1,1,1), sendcounts, MPI_RPREC,&
!      u_tot(1,1,1), sendcounts, displs,       &
!      MPI_RPREC, 0, comm, ierr)
!   call mpi_gatherv (v(1,1,1), sendcounts, MPI_RPREC,&
!      v_tot(1,1,1), sendcounts, displs,       &
!      MPI_RPREC, 0, comm, ierr)
!   call mpi_gatherv (w(1,1,1), sendcounts, MPI_RPREC,&
!      w_tot(1,1,1), sendcounts, displs,       &
!      MPI_RPREC, 0, comm, ierr)
!   if (theta_flag) then
!     call mpi_gatherv (theta(1,1,1), sendcounts, MPI_RPREC,&
!        theta_tot(1,1,1), sendcounts, displs,       &
!        MPI_RPREC, 0, comm, ierr)
!   endif
!   !$endif

!   if ((USE_MPI .and. coordz == 0)) then
!     do iSpl = 1, nVSpl
!       !-BC generate file name
!       if (theta_flag) then
!         write(fname, '(a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i8.8,a)')&
!            path//'output/sample/sampleVS_ix',&
!            xVSplS(iSpl), '-', xVSplE(iSpl),&
!            '_iy', yVSplS(iSpl), '-', yVSplE(iSpl),&
!            '_iz', zVSplS(iSpl), '-', zVSplE(iSpl),&
!            '_tt', jt_total, '.out'
!       else
!         write(fname, '(a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i8.8,a)')&
!            path//'output/sample/sampleV_ix',&
!            xVSplS(iSpl), '-', xVSplE(iSpl),&
!            '_iy', yVSplS(iSpl), '-', yVSplE(iSpl),&
!            '_iz', zVSplS(iSpl), '-', zVSplE(iSpl),&
!            '_tt', jt_total, '.out'
!       endif

!       ixs = xVSplS(iSpl)
!       ixe = xVSplE(iSpl)
!       iys = yVSplS(iSpl)
!       iye = yVSplE(iSpl)
!       izs = zVSplS(iSpl)
!       ize = zVSplE(iSpl)

!       !-BC write data
!       open(1, file = fname, form='unformatted')
!       if (theta_flag) then
!         write (1) u_tot(ixs:ixe, iys:iye, izs:ize),&
!            v_tot(ixs:ixe, iys:iye, izs:ize),&
!            w_tot(ixs:ixe, iys:iye, izs:ize),&
!            theta_tot(ixs:ixe, iys:iye, izs:ize)
!       else
!         write (1) u_tot(ixs:ixe, iys:iye, izs:ize),&
!            v_tot(ixs:ixe, iys:iye, izs:ize),&
!            w_tot(ixs:ixe, iys:iye, izs:ize)
!       endif
!       close(1)
!     enddo

!   elseif (.not. USE_MPI) then
!     do iSpl = 1, nCSpl
!       !-BC generate file name
!       if (theta_flag) then
!         write(fname, '(a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i8.8,a)')&
!            path//'output/sample/sampleVS_ix',&
!            xVSplS(iSpl), '-', xVSplE(iSpl),&
!            '_iy', yVSplS(iSpl), '-', yVSplE(iSpl),&
!            '_iz', zVSplS(iSpl), '-', zVSplE(iSpl),&
!            '_tt', jt_total, '.out'
!       else
!         write(fname, '(a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i8.8,a)')&
!            path//'output/sample/sampleV_ix',&
!            xVSplS(iSpl), '-', xVSplE(iSpl),&
!            '_iy', yVSplS(iSpl), '-', yVSplE(iSpl),&
!            '_iz', zVSplS(iSpl), '-', zVSplE(iSpl),&
!            '_tt', jt_total, '.out'
!       endif

!       !-BC write data
!       open(1, file = fname, form='unformatted')
!       ixs = xVSplS(iSpl)
!       ixe = xVSplE(iSpl)
!       iys = yVSplS(iSpl)
!       iye = yVSplE(iSpl)
!       izs = zVSplS(iSpl)
!       ize = zVSplE(iSpl)
!       if (theta_flag) then
!         write (1) u(ixs:ixe, iys:iye, izs:ize),&
!            v(ixs:ixe, iys:iye, izs:ize),&
!            w(ixs:ixe, iys:iye, izs:ize),&
!            theta(ixs:ixe, iys:iye, izs:ize)
!       else
!         write (1) u(ixs:ixe, iys:iye, izs:ize),&
!            v(ixs:ixe, iys:iye, izs:ize),&
!            w(ixs:ixe, iys:iye, izs:ize)
!       endif
!       close(1)
!     enddo

!   endif

!     end subroutine checkVSpl
! !-----------------------------------------------------------------------
! !   
! !-----------------------------------------------------------------------
!     subroutine checkCSpl()
!   use param
!   use sim_param, only : u, pcon

!   implicit none

!   character(80) :: fname
!   integer :: iSpl, ixs, ixe, iys, iye, izs, ize

!   real(8), dimension(:,:,:,:), allocatable :: PConSpl
!   integer :: ipcon

!   !$if ($MPI)
!   integer :: sendcounts(npz)
!   integer :: recvcounts(npz)
!   integer :: displs(npz)
!   !$endif

!   real (rprec), dimension(nxt, nyt, 1:nzt, npcon) :: PCon_tot

! !---------------------------------------------------------------------
! !-BC gather the data
!   !$if ($MPI)
!   sendcounts = size (u(:,:,1:nz))
!   recvcounts = size (u(:,:,1:nz-1))
!   displs = coord_of_rank * recvcounts
!   if (PCon_FLAG) then
!     do ipcon=1,npcon
!       call mpi_gatherv (PCon(1,1,1,ipcon), sendcounts, MPI_RPREC,&
!          PCon_tot(1,1,1,ipcon), sendcounts, displs,       &
!          MPI_RPREC, 0, comm, ierr)
!     enddo
!   endif
!   !$endif

!   if ((USE_MPI .and. coordz == 0)) then
!     do iSpl = 1, nCSpl
!       !-BC generate the file name
!       write(fname, '(a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i8.8,a)')&
!          path//'output/sample/sampleC_ix',&
!          xCSplS(iSpl), '-', xCSplE(iSpl),&
!          '_iy', yCSplS(iSpl), '-', yCSplE(iSpl),&
!          '_iz', zCSplS(iSpl), '-', zCSplE(iSpl),&
!          '_tt', jt_total, '.out'

!       !-BC allocate the variables
!       ixs = xCSplS(iSpl)
!       ixe = xCSplE(iSpl)
!       iys = yCSplS(iSpl)
!       iye = yCSplE(iSpl)
!       izs = zCSplS(iSpl)
!       ize = zCSplE(iSpl)
!       allocate(PConSpl(ixs:ixe, iys:iye, izs:ize, npcon))

!       !-gather the output data
!       do ipcon = 1, npcon
!         PConSpl(ixs:ixe, iys:iye ,izs:ize,ipcon) = &
!            PCon_tot(ixs:ixe, iys:iye,izs:ize,ipcon)
!       enddo

!       !-BC output the data
!       open(1, file = fname, form='unformatted')
!       write(1) PConSpl
!       close(1)
!       deallocate(PConSpl)
!     enddo

!   elseif (.not. USE_MPI) then
!     do iSpl = 1, nCSpl
!       !-BC generate the file name
!       write(fname, '(a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i3.3,a,i8.8,a)')&
!          path//'output/sample/sampleC_ix',&
!          xCSplS(iSpl), '-', xCSplE(iSpl),&
!          '_iy', yCSplS(iSpl), '-', yCSplE(iSpl),&
!          '_iz', zCSplS(iSpl), '-', zCSplE(iSpl),&
!          '_tt', jt_total, '.out'

!       !-BC allocate the variables
!       ixs = xCSplS(iSpl)
!       ixe = xCSplE(iSpl)
!       iys = yCSplS(iSpl)
!       iye = yCSplE(iSpl)
!       izs = zCSplS(iSpl)
!       ize = zCSplE(iSpl)
!       allocate(PConSpl(ixs:ixe, iys:iye, izs:ize, npcon))

!       !-gather the output data
!       do ipcon = 1, npcon
!         PConSpl(ixs:ixe, iys:iye ,izs:ize,ipcon) = &
!            PCon(ixs:ixe, iys:iye,izs:ize,ipcon)
!       enddo

!       !-BC output the data
!       open(1, file = fname, form='unformatted')
!       write(1) PConSpl
!       close(1)
!       deallocate(PConSpl)
!     enddo
!   endif

!     end subroutine checkCSpl
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------    
    subroutine output_final(lun_opt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--lun_opt gives option for writing to different unit, and is used by
!  inflow_write
!--assumes the unit to write to (lun_default or lun_opt is already
!  open for sequential unformatted writing
!--this routine also closes the unit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    ! ---
    integer, intent(in), optional :: lun_opt  !--if present, write to unit lun
    ! ---
    integer, parameter :: lun_default = 17
    integer :: lun
    ! ---
    if (present (lun_opt)) then
        lun = lun_opt
    else
        lun = lun_default
    end if

    call checkpoint_final (lun)

    if (cumulative_time) then
    !--only do this for true final output, not intermediate recording
        open (1, file=fcumulative_time)
        !  if ((DYN_init .ne. 1) .AND. (GABLS_diurnal_test)) jt_total=jt_total_init
        write (1, *) nums
        close (1)
    end if
    ! ---
    end subroutine output_final
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    subroutine checkpoint_final (fid_out)
    use param
    use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta, pcon
    use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX, F_XX,F_KX2,F_XX2
    use scalars_module, only : RHS_T,RHS_PCon,sgs_t3,deposition,Real_dep,Kc_t
    use bottombc, only : psi_m

    implicit none
    ! ---
    integer, intent (in) :: fid_out
    ! ---
    integer:: ipcon
    logical :: flag_dir
    integer(MPI_OFFSET_KIND) :: disp0, disp3d, disp2d
    character(100) :: cmd
    character(80) :: fname
    ! --- create directory for restart file
    inquire(directory=path_restart, exist=flag_dir)
    if (.not. flag_dir) then
        write(cmd,'(a)') 'mkdir '//path_restart
        call system(trim(cmd))
    end if
    ! ---
    disp3d = np  * sizeof(u(1:nxt, 1:nynpy, 1:nz-1))
    disp2d = npy * sizeof(u(1:nxt, 1:nynpy, 1))

    ! --- output velocity field
    write(fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
    disp0 = 0
    call write_file_mpi_3d(u(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(v(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(w(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(RHSx(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(RHSy(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(RHSz(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(Cs_opt2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_LM(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_MM(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_QN(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    disp0 = disp0+disp3d
    call write_file_mpi_3d(F_NN(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)

    ! --- output the temperature field
    if (theta_flag) then
        write(fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
        disp0 = 0
        call write_file_mpi_3d(theta(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(RHS_T(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_xy(sgs_t3(1:nxt, 1:nynpy, 1), fname, disp0)
        disp0 = disp0+disp2d
        call write_file_mpi_xy(psi_m, fname, disp0)
    end if

    if (PCon_flag) then
        ! --- output the concentration
        write(fname, '(a,i8.8,a)') path_restart//'con_tt', nums, '.out'
        disp0 = 0
        do ipcon = 1, npcon
            call write_file_mpi_3d(Pcon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0)
            call write_file_mpi_3d(RHS_Pcon(1:nxt, 1:nynpy, 1:nz-1, ipcon), fname, disp0+npcon*disp3d)
            disp0 = disp0 + disp3d
        end do

        ! --- output related variables for concentration
        write(fname, '(a,i8.8,a)') path_restart//'con_diff_tt', nums, '.out'
        call write_file_mpi_xy(deposition, fname, 0)
        disp0 = disp2d
        call write_file_mpi_xy(Real_dep, fname, disp0)
        disp0 = disp0+disp2d
        call write_file_mpi_xy(Kc_t(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_KX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_XX(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_KX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
        disp0 = disp0+disp3d
        call write_file_mpi_3d(F_XX2(1:nxt, 1:nynpy, 1:nz-1), fname, disp0)
    end if
    ! ---
    end subroutine checkpoint_final
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
    ! subroutine checkpoint_final_old (fid_out)
    ! !use param, only : nxt,nyt,nz,nzt,theta_flag,PCon_FLAG
    ! use param
    ! use sim_param, only : u, v, w, RHSx, RHSy, RHSz, theta, pcon
    ! use sgsmodule, only : Cs_opt2, F_LM, F_MM, F_QN, F_NN, F_KX, F_XX,F_KX2,F_XX2
    ! use scalars_module, only : RHS_T,RHS_PCon,sgs_t3,deposition,Real_dep,Kc_t
    ! use bottombc, only : psi_m

    ! implicit none
    ! ! ---
    ! integer, intent (in) :: fid_out
    ! integer:: ipcon, jx
    ! character (80) :: fname
    ! logical :: flag_dir
    ! character (100) :: cmd

    ! !integer, dimension(npz) :: sendcounts, recvcounts, displs

    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: u_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: v_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: w_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: RHSx_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: RHSy_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: RHSz_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt) :: theta_tot
    ! real(kind=rprec), dimension(nxt, nyt, 1:nzt, npcon) :: PCon_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: Cs_opt2_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_LM_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_MM_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_QN_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_NN_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: RHS_T_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt,npcon) :: RHS_PCon_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: Kc_t_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_KX_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_XX_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_KX2_tot
    ! real(kind=rprec), dimension(nxt,nyt,nzt) :: F_XX2_tot

    ! real(kind=rprec), dimension(nxt, nyt) :: sgs_t3_s_hor
    ! real(kind=rprec), dimension(nxt, nyt) :: psi_m_hor
    ! ! --- exchange ghost node data to prevent BOGUS value when
    ! ! >>call output_final (Bicheng Chen)
    ! if ( coordz == npz - 1 ) then
    !     RHSx(:, :, nz) = RHSx(:, :, nz - 1)
    !     RHSy(:, :, nz) = RHSy(:, :, nz - 1)
    !     RHSz(:, :, nz) = RHSz(:, :, nz - 1)
    ! end if
    ! call mpi_sendrecv(RHSx(1, 1, 1),  nxt*nynpy, MPI_RPREC, down, 1, &
    !                   RHSx(1, 1, nz), nxt*nynpy, MPI_RPREC, up,   1, &
    !                   comm, status, ierr)
    ! call mpi_sendrecv(RHSy(1, 1, 1),  nxt*nynpy, MPI_RPREC, down, 2, &
    !                   RHSy(1, 1, nz), nxt*nynpy, MPI_RPREC, up,   2, &
    !                   comm, status, ierr)
    ! call mpi_sendrecv(RHSz(1, 1, 1),  nxt*nynpy, MPI_RPREC, down, 3, &
    !                   RHSz(1, 1, nz), nxt*nynpy, MPI_RPREC, up,   3, &
    !                   comm, status, ierr)
    ! ! ---
    ! call collect_data_3d(u, u_tot)
    ! call collect_data_3d(v, v_tot)
    ! call collect_data_3d(w, w_tot)
    ! call collect_data_3d(RHSx, RHSx_tot)
    ! call collect_data_3d(RHSy, RHSy_tot)
    ! call collect_data_3d(RHSz, RHSz_tot)
    ! call collect_data2_3d(Cs_opt2, Cs_opt2_tot)
    ! call collect_data2_3d(F_LM, F_LM_tot)
    ! call collect_data2_3d(F_MM, F_MM_tot)
    ! call collect_data2_3d(F_QN, F_QN_tot)
    ! call collect_data2_3d(F_NN, F_NN_tot)
    
    ! if (theta_flag) then
    !     call collect_data_3d(theta, theta_tot)
    !     call collect_data_3d(RHS_T, RHS_T_tot)

    !     do jx = 1, nxt
    !         call collect_data_y(sgs_t3(jx,:,1), sgs_t3_s_hor(jx,:))
    !     end do
    !     do jx = 1, nxt
    !         call collect_data_y(psi_m(jx,:), psi_m_hor(jx,:))
    !     end do
    ! end if

    ! if (PCon_flag) then
    !     do ipcon = 1, npcon
    !         call collect_data4_3d(PCon(:,:,:,ipcon), PCon_tot(:,:,:,ipcon))
    !         call collect_data4_3d(RHS_PCon(:,:,:,ipcon), RHS_PCon_tot(:,:,:,ipcon))
    !     end do
    !     call collect_data4_3d(Kc_t, Kc_t_tot)
    !     call collect_data3_3d(F_KX, F_KX_tot)
    !     call collect_data3_3d(F_XX, F_XX_tot)
    !     call collect_data3_3d(F_KX2, F_KX2_tot)
    !     call collect_data3_3d(F_XX2, F_XX2_tot)
    ! end if

    ! if ( rank == 0 ) then
    ! !- create directory for restart file
    !     inquire(directory=path_restart, exist=flag_dir)
    ! !inquire(file=path_restart, exist=flag_dir)
    !     if (.not. flag_dir) then
    !         write(cmd,'(a)') 'mkdir '//path_restart
    !         call system(trim(cmd))
    !     end if

    ! !- output velocity field
    !     ! write(fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
    !     ! open(fid_out, file=fname, form='unformatted')
    !     ! write(fid_out) u_tot(:,:,1:nzt), v_tot(:,:,1:nzt), w_tot(:,:,1:nzt)
    !     ! close(fid_out)
    !     write(fname, '(a,i8.8,a)') path_restart//'vel_tt', nums, '.out'
    !     open(fid_out, file=fname, form='unformatted')
    !     write(fid_out) u_tot(:,:,1:nzt), v_tot(:,:,1:nzt),&
    !         w_tot(:,:,1:nzt), RHSx_tot(:, :, 1:nzt), RHSy_tot(:, :, 1:nzt),&
    !         RHSz_tot(:, :, 1:nzt), Cs_opt2_tot,&
    !         F_LM_tot, F_MM_tot, F_QN_tot, F_NN_tot
    !     close(fid_out)

    ! !- output the temperature field
    !     if (theta_flag) then
    !         write(fname, '(a,i8.8,a)') path_restart//'temp_tt', nums, '.out'
    !         open(fid_out, file=fname, form='unformatted')
    !         write(fid_out) theta_tot(:,:,1:nzt), RHS_T_tot(:,:,1:nzt),  &
    !                        sgs_t3_s_hor(:,:), psi_m_hor
    !         close(fid_out)
    !     end if

    !     !- output the concentration field
    !     if (PCon_flag) then
    !         !--- output the concentration
    !         write(fname, '(a,i8.8,a)') path_restart//'con_tt', nums, '.out'
    !         open(fid_out, file=fname, form='unformatted')
    !         write(fid_out) PCon_tot(:,:,1:nzt,1:npcon),&
    !                        RHS_PCon_tot(:,:,1:nzt,1:npcon)
    !         close(fid_out)
    !         !--- output related variables for concentration
    !         write(fname, '(a,i8.8,a)') path_restart//'con_diff_tt', nums, '.out'
    !         open(fid_out, file=fname, form='unformatted')
    !         write(fid_out) deposition(:,:), Real_dep(:,:),&
    !                        Kc_t_tot, F_KX_tot, F_XX_tot, F_KX2_tot, F_XX2_tot
    !         close(fid_out)
    !     end if
    ! end if
    ! ! ---
    ! end subroutine checkpoint_final_old
! ---
end module io