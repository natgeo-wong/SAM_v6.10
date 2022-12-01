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

subroutine wtg_decomp2022(masterproc, nzm, nz, z, &
                          theta_ref, theta_model, tabs_model, tau_wtg, &
                          tau_h, tau_f, &
                          dowtgLBL, boundstatic, dthetadz_min, w_wtg, &
                          wwtgr, wwtgh, wwtgf, wwtgc)

implicit none

! ======= inputs =======
logical, intent(in) :: masterproc ! For printing out to logs
integer, intent(in) :: nzm ! number of model levels
integer, intent(in) :: nz  ! number of vertical levels
real, intent(in) :: z(nz) ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: theta_ref(nzm) ! reference potential temperature sounding in K
real, intent(in) :: theta_model(nzm) ! model potential temperature profile in K (domain-mean for LES)
real, intent(in) :: tabs_model(nzm) ! model temperature profile in K (domain-mean for LES)

! WTG potential temperature relaxation timescale (tau_wtg)
!   default is 1 day^-1 (tau_wtg = 1/86400 s^-1)
real, intent(in) :: tau_wtg ! potential temperature relaxation timescale (s^-1)
real, intent(in) :: tau_h ! potential temperature relaxation timescale (s^-1)
real, intent(in) :: tau_f ! potential temperature relaxation timescale (s^-1)

logical, intent(in) :: dowtgLBL    ! Calculate w_wtg at boundary layer instead of linear interpolation
logical, intent(in) :: boundstatic ! Restrict lower bound for static stability
real, intent(in) :: dthetadz_min   ! if boundstatic = .true., what is the minimum bound?

! ======= output =======
real, intent(out) :: w_wtg(nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.
real, intent(out) :: wwtgr(nzm) ! Raw w_wtg, assuming constant = 1
real, intent(out) :: wwtgh(nzm) ! Half-sine component of w_wtg
real, intent(out) :: wwtgf(nzm) ! Full-sine component of w_wtg
real, intent(out) :: wwtgc(nzm) ! Coefficient for the decompositions


! ======= local variables =======
integer :: k
integer :: ktrop
integer :: kbl
real :: min_temp ! temporary variable used to find cold point of model sounding.
real :: ztrop ! Height of tropopause level (m)
real :: dthetadz ! Static Stability
real, parameter :: pi = 3.141592653589793 ! from MATLAB, format long.
real :: theta_diff ! Static Stability
real :: ttheta_h ! decomposition constant into half-sine curve
real :: ttheta_f ! decomposition constant into full-sine curve

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
if(.NOT.dowtgLBL) then
  do k = nzm,1,-1
    if (z(k)>1000) then
      kbl = k
    end if
  end do
end if

ttheta_h = (theta_model(kbl)-theta_ref(kbl)) * sin(pi*z(kbl)  /ztrop) * z(kbl)
ttheta_f = (theta_model(kbl)-theta_ref(kbl)) * sin(pi*z(kbl)*2/ztrop) * z(kbl)

do k = (kbl+1),ktrop

  ttheta_h = ttheta_h + ((theta_model(k)   - theta_ref(k))   * sin(pi*z(k)  /ztrop)  &
                      +  (theta_model(k-1) - theta_ref(k-1)) * sin(pi*z(k-1)/ztrop)) &
                      *  (z(k)-z(k-1))

  ttheta_f = ttheta_f + ((theta_model(k)   - theta_ref(k))   * sin(pi*z(k)*2  /ztrop)  &
                      +  (theta_model(k-1) - theta_ref(k-1)) * sin(pi*z(k-1)*2/ztrop)) &
                      *  (z(k)-z(k-1))

end do

ttheta_h = ttheta_h / ztrop
ttheta_f = ttheta_f / ztrop

if(.NOT.dowtgLBL) then
  do k = kbl,ktrop

    dthetadz = (theta_model(k+1)-theta_model(k-1)) / (z(k+1)-z(k-1))
    if (boundstatic.AND.(dthetadz.lt.dthetadz_min).AND.(z(k)>5000)) dthetadz = dthetadz_min
    ! According to Raymond and Zeng (2005) model feedbacks in the upper troposphere can
    ! otherwise result in very weak static stabilities and unrealistically
    ! large values of w_wtg
    wwtgr(k) = theta_model(k) - theta_ref(k)           * tau_wtg / dthetadz
    wwtgh(k) = tau_h * ttheta_h * sin(pi*z(k)  /ztrop) * tau_wtg / dthetadz
    wwtgf(k) = tau_f * ttheta_f * sin(pi*z(k)*2/ztrop) * tau_wtg / dthetadz
    w_wtg(k) = wwtgh(k) + wwtgf(k)

  end do
  do k = 1,(kbl-1)
    wwtgr(k) = wwtgr(kbl) * z(k) / z(kbl)
    wwtgh(k) = wwtgh(kbl) * z(k) / z(kbl)
    wwtgf(k) = wwtgf(kbl) * z(k) / z(kbl)
    w_wtg(k) = w_wtg(kbl) * z(k) / z(kbl)
  end do
else
  do k = 2,ktrop

    dthetadz = (theta_model(k+1)-theta_model(k-1)) / (z(k+1)-z(k-1))
    if (boundstatic.AND.(dthetadz.lt.dthetadz_min)) dthetadz = dthetadz_min
    ! According to Raymond and Zeng (2005) model feedbacks in the upper troposphere can
    ! otherwise result in very weak static stabilities and unrealistically
    ! large values of w_wtg
    wwtgr(k) = theta_model(k) - theta_ref(k)           * tau_wtg / dthetadz
    wwtgh(k) = tau_h * ttheta_h * sin(pi*z(k)  /ztrop) * tau_wtg / dthetadz
    wwtgf(k) = tau_f * ttheta_f * sin(pi*z(k)*2/ztrop) * tau_wtg / dthetadz
    w_wtg(k) = wwtgh(k) + wwtgf(k)

  end do
  dthetadz = (theta_model(2)-theta_model(1)) / (z(2) - z(1))
  if (boundstatic.AND.(dthetadz.lt.dthetadz_min)) dthetadz = dthetadz_min
  wwtgr(1) = theta_model(1) - theta_ref(1)           * tau_wtg / dthetadz
  wwtgh(1) = tau_h * ttheta_h * sin(pi*z(1)  /ztrop) * tau_wtg / dthetadz
  wwtgf(1) = tau_f * ttheta_f * sin(pi*z(1)*2/ztrop) * tau_wtg / dthetadz
  w_wtg(1) = wwtgh(1) + wwtgf(1)
end if

wwtgc(1) = ttheta_h
wwtgc(2) = ttheta_f

end subroutine wtg_decomp2022
