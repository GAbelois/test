!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module canopy
!-----------------------------------------------------------------------
! canopy-related parameters
!-----------------------------------------------------------------------
use types, only:rprec
implicit none
save
public
! --- the upstream, downstream, south, and north domain size of the canopy
real(rprec), parameter :: lwest = 150._rprec
real(rprec), parameter :: least = 250._rprec
real(rprec), parameter :: lsouth = 0._rprec
real(rprec), parameter :: lnorth = 0._rprec
real(rprec), parameter :: LBB    = 26.0_rprec            ! the spacing between adjacent backbones
real(rprec), parameter :: LGL    = 8.0_rprec            ! the length of the growth line
! --- nutrient uptake related parameters
real(rprec), parameter :: Vmax = 752._rprec/3600._rprec  ! maximum uptake rate, units: umoles m-2 s-1
real(rprec), parameter :: Km   = 10.2_rprec*1000._rprec  ! half-saturation constant, units: umoles/m-3, unit conversion from liter to m^3
! ---
real(kind=rprec), dimension(:), allocatable :: nutrient_inflow
! ---
end module 
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 