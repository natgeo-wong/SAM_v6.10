! This subroutine solves for perturbation omega based on WTG approximation as
! described in the appendix of Blossey, Bretherton & Wyant (2009):
!
!   http://dx.doi.org/10.3894/JAMES.2009.1.8
!
!  - original version: Peter Blossey (pblossey at gmail dot com), 
!                      September 2009 based on earlier implementaion in SAM.
!
! Note 1: this version handles the situation where the model top is
!   below the tropopause gracefully.  The boundary condition at the
!   tropopause (omega'=0) is enforced by adding additional levels
!   between the model top and the tropopause (whose height is assumed
!   to be 100 hPa if it is above the model top).
!
! Note 2: Below the tropopause, the temperature will be modified
!   by the WTG method through large-scale vertical motion.
!   Above the tropopause, the temperature should be nudged 
!   to the observed profile on a longish (~2 day?) timescale.
!
! Note 3: The wrapper routine allows the code to handle models
!   that index their pressure, temperature and moisture soundings
!   from either the surface-upwards or the top-downwards.  The driver
!   routine assumes that these soundings are indexed from the surface
!   upwards as in SAM, the model for which this routine was first
!   written.

subroutine wtg_qjrms2005(nzm, nz, z, theta_ref, theta_model, tabs_model, t_wtg, w_wtg)

implicit none

! ======= inputs =======
integer, intent(in) :: nzm ! number of model levels
integer, intent(in) :: nz  ! number of vertical levels
real, intent(in) :: z(nz) ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: theta_ref(nzm) ! reference potential temperature sounding in K
real, intent(in) :: theta_model(nzm) ! model potential temperature sounding in K (domain-mean for LES)

! WTG potential temperature relaxation timescale (t_wtg)
!   default is 1 day (t_wtg = 86400 s)
real, intent(in) :: t_wtg     ! potential temperature relaxation timescale (seconds)

! ======= output =======
real, intent(out) :: w_wtg(nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.

real :: ztrop ! Height of tropopause level (m)
real, parameter :: pi = 3.141592653589793 ! from MATLAB, format long.

if (z(nz) < 1.e4) then

  if(masterproc) then
    write(*,*) '******* ERROR in WTG Scheme of Raymond and Zeng (2005) *******'
    write(*,*) 'The model top is less than 10 km high. The WTG module of'
    write(*,*) 'Raymond and Zeng (2005) has not yet been configured to run on'
    write(*,*) 'a short domain.'
  end if
  call task_abort()

end if

! ===== find index of cold point tropopause in vertical. =====
! reverse pressure coordinate, and find index
!   of cold point tropopause in the vertical.
ktrop = nzm+1 ! default is top of model/atmosphere (counting from surface)
min_temp = tabs_model(nzm) 
do k = 1,nzm
  if(tabs_model(k).lt.min_temp) then
    ktrop = k
    min_temp = tabs_model(k)
  end if
end do

ztrop = z(ktrop)

! ===== find index of boundary layer top =====
! the boundary layer is defined to be the bottom 1km layer of the atmosphere
kbl = 1 ! set to be the model bottom
do k = nzm,1,-1
  if (z(k)>1000) then
    kbl = k-1
  end if
end do

do k = kbl,ktrop

  w_wtg(k) = sin(pi*z(k)/ztrop) * (theta_model(k)-theta_ref(k)) / t_wtg * &
              (z(k+1)-z(k-1)) / (theta_model(k+1)-theta_model(k-1))

end do

end subroutine wtg_james2009