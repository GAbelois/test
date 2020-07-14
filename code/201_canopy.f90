!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------
subroutine init_lad()
!--------------------------------------------------------------------!
!    initialize the leaf area density profile                                          
!--------------------------------------------------------------------! 
use types, only:rprec
use param
use canopy, only:lwest, least, lsouth, lnorth, LBB, LGL
use ibm, only: z_tpg
implicit none
! ---
integer :: jx, jy, jz, jh
real(rprec) :: tmp
! ---
if (flag_canopy) then
    allocate(a_leafz_uv(nxt, nynpy, nz))
    allocate(a_leafz_w(nxt, nynpy, nz))
    a_leafz_uv = 0._rprec
    a_leafz_w  = 0._rprec

    if (ocean_flag) then
        if (use_field_scale) then
            select case (canopy_conf)
            case ('stripe')
                do jz = 1, nz-1
                do jy = 1, nynpy
                do jx = 1, nxt
                    if (x(jx) .ge. lwest/z_i  .and. x(jx) .le. lx_tot - least/z_i  .and.    &
                        y(jy) .ge. lsouth/z_i .and. y(jy) .le. ly_tot - lnorth/z_i .and.    &
                        mod(x(jx)*z_i-lwest, LBB) .ge. 0._rprec .and. mod(x(jx)*z_i-lwest, LBB) .lt. LGL-dx*z_i/2._rprec ) then
                        
                        do jh = 1, nc-1
                            if ( zuvp(jz) .ge. hgt_leaf(jh)/z_i .and. zuvp(jz) .le. hgt_leaf(jh+1)/z_i ) then
                                tmp = (hgt_leaf(jh+1) - zuvp(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                                a_leafz_uv(jx, jy, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                            end if

                            if ( zw(jz) .ge. hgt_leaf(jh)/z_i .and. zw(jz) .le. hgt_leaf(jh+1)/z_i ) then
                                tmp = (hgt_leaf(jh+1) - zw(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                                a_leafz_w(jx, jy, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                            end if
                        end do
                    end if
                end do
                end do
                end do
            case ('block')
                do jz = 1, nz-1
                    do jy = 1, nynpy
                    do jx = 1, nxt
                        if ( x(jx) .ge. lwest/z_i  .and. x(jx) .le. lx_tot - least/z_i  .and.    &
                             y(jy) .ge. lsouth/z_i .and. y(jy) .le. ly_tot - lnorth/z_i ) then       
                            do jh = 1, nc-1
                                if ( zuvp(jz) .ge. hgt_leaf(jh)/z_i .and. zuvp(jz) .le. hgt_leaf(jh+1)/z_i ) then
                                    tmp = (hgt_leaf(jh+1) - zuvp(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                                    a_leafz_uv(jx, jy, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                                end if
    
                                if ( zw(jz) .ge. hgt_leaf(jh)/z_i .and. zw(jz) .le. hgt_leaf(jh+1)/z_i ) then
                                    tmp = (hgt_leaf(jh+1) - zw(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                                    a_leafz_w(jx, jy, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                                end if
                            end do
                        end if
                    end do
                    end do
                    end do
            case default
                write(*,*) "Canopy configutation NOT supported!"
                stop
            end select
        end if
    else
        if (flag_ibm) then
        ! If use immersed boundary method, the altitude is with respect 
        ! to the phi=0 (immersed boundary)
            do jz = 1, nz-1
                do jy = 1, nynpy
                do jx = 1, nxt
                    if (zuvp(jz)-z_tpg(jx,jy) .ge. hgt_leaf(jh)/z_i .and.   &
                        zuvp(jz)-z_tpg(jx,jy) .le. hgt_leaf(jh+1)/z_i) then
                        tmp = (hgt_leaf(jh+1) - zuvp(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                        a_leafz_uv(jx, jy, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i
                    end if

                    if (zw(jz)-z_tpg(jx,jy) .ge. hgt_leaf(jh)/z_i .and.     &
                        zw(jz)-z_tpg(jx,jy) .lt. hgt_leaf(jh+1)/z_i) then
                        tmp = (hgt_leaf(jh+1) - zw(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                        a_leafz_w(:, :, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                    end if
                end do
                end do
            end do
        else
            if (use_field_scale) then
                do jz = 1, nz-1
                    do jh = 1, nc-1
                        if ( zuvp(jz) .ge. hgt_leaf(jh)/z_i .and. zuvp(jz) .le. hgt_leaf(jh+1)/z_i ) then
                            tmp = (hgt_leaf(jh+1) - zuvp(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                            a_leafz_uv(:, :, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                        end if

                        if ( zw(jz) .ge. hgt_leaf(jh)/z_i .and. zw(jz) .le. hgt_leaf(jh+1)/z_i ) then
                            tmp = (hgt_leaf(jh+1) - zw(jz)*z_i) / (hgt_leaf(jh+1) - hgt_leaf(jh))
                            a_leafz_w(:, :, jz) = tmp*lad_leaf(jh)*z_i + (1._rprec-tmp)*lad_leaf(jh+1)*z_i 
                        end if
                    end do
                end do
            end if
        end if
    end if
end if
! ---
end subroutine init_lad
!--------------------------------------------------------------------!
!                                              
!--------------------------------------------------------------------!
subroutine canopy_model 
!-----------------------------------------------------------------------
!   pressure gradient and canopy drag
!-----------------------------------------------------------------------
use types, only:rprec
use param
use sim_param
implicit none
! ---
real(kind=rprec) :: V_uv, V_w
integer :: jx, jy, jz
! ---

dragx = 0._rprec
dragy = 0._rprec
dragz = 0._rprec

if (flag_canopy) then
    if (use_field_scale) then
        do jz = 1, nz - 1
        do jy = 1, nynpy
        do jx = 1, nxt
            V_uv = (  u(jx,jy,jz)**2._rprec &
                 +    v(jx,jy,jz)**2._rprec &
                 +  ((w(jx,jy,jz)+w(jx,jy,jz+1))/2._rprec)**2._rprec)**0.5_rprec
            V_w  = (((u(jx,jy,jz)+u(jx,jy,jz-1))/2._rprec)**2._rprec &
                 +  ((v(jx,jy,jz)+v(jx,jy,jz-1))/2._rprec)**2._rprec &
                 +    w(jx,jy,jz)**2._rprec                         )**0.5_rprec

            dragx(jx, jy, jz) = -0.5_rprec*Cd_const*a_leafz_uv(jx, jy, jz)*V_uv**(vel_exp-1.)*u(jx, jy, jz) 
            dragy(jx, jy, jz) = -0.5_rprec*Cd_const*a_leafz_uv(jx, jy, jz)*V_uv**(vel_exp-1.)*v(jx, jy, jz)
            dragz(jx, jy, jz) = -0.5_rprec*Cd_const*a_leafz_w(jx, jy, jz)*V_w**(vel_exp-1.)*w(jx, jy, jz)
        end do
        end do
        end do
    end if 
end if
! ---
end subroutine canopy_model
!-----------------------------------------------------------------------
!   
!-----------------------------------------------------------------------

    



