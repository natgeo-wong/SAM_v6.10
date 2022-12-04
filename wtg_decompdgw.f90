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

subroutine wtg_decompdgw(masterproc, nzm, nz, z, &
                          pres_ref, tabs_ref, qv_ref, tabs_model, qv_model, qcond_model, &
                          lambda_wtg, am_wtg, &
                          wtgscale_vertmodenum, wtgscale_vertmodescl, &
                          w_wtg, wwtgc, ktrop)

implicit none

! ======= inputs =======
logical, intent(in) :: masterproc ! For printing out to logs
integer, intent(in) :: nzm ! number of model levels
integer, intent(in) :: nz  ! number of vertical levels
real, intent(in) :: z(nz) ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: pres_ref(nzm) ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: tabs_ref(nzm) ! reference temperature sounding in K
real, intent(in) :: qv_ref(nzm)   ! reference water vapor sounding in kg/kg
real, intent(in) :: tabs_model(nzm) ! model temperature sounding in K (domain-mean for LES)
real, intent(in) :: qv_model(nzm)   ! model water vapor sounding in kg/kg (domain-mean for LES)
real, intent(in) :: qcond_model(nzm) ! model condensate (liq+ice) sounding in kg/kg (domain-mean for LES)

real, intent(in) :: lambda_wtg ! WTG length scale in m (JAMES2009 value = 650.e6 m)
real, intent(in) :: am_wtg ! WTG momentum damping rate in 1/s at p/pref (default = 1./86400. /s)
integer, intent(in) :: wtgscale_vertmodenum ! response scaling for 1st vertical mode
real, intent(in) :: wtgscale_vertmodescl(wtgscale_vertmodenum) ! response scaling for 2nd vertical mode

! ======= output =======
real, intent(out) :: w_wtg(nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.
real, intent(out) :: wwtgc(nz) ! Coefficient for the decompositions
integer, intent(out) :: ktrop ! index of interface just above the cold point.


! ======= local variables =======
integer :: k
integer :: kbl
integer :: inum
real :: min_temp ! temporary variable used to find cold point of model sounding.
real :: ztrop ! Height of tropopause level (m)
real :: dthetadz ! Static Stability
real, parameter :: grav = 9.81 ! from MATLAB, format long.
real, parameter :: pi   = 3.141592653589793 ! from MATLAB, format long.
real, parameter :: rgas = 287. ! Gas constant for dry air, J/kg/K
real :: thetad1 ! Static Stability
real :: thetad2 ! Static Stability
real :: tv_model(nzm) ! virtual temperature of model sounding in K
real :: tv_ref(nzm) !  virtual temperature of reference sounding in K

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

! ===== find virtual temperature =====
do k = 1,ktrop
  ! virtual temperature of model sounding
  tv_model(k) = tabs_model(k)*( 1. + 0.61*qv_model(k) - qcond_model(k) )

  ! virtual temperature of reference sounding (assume unsaturated)
  tv_ref(k) = tabs_ref(k)*( 1. + 0.61*qv_ref(k) )
end do

! ===== calculate vertical mode coefficients =====

thetad1 = (tv_model(1) - tv_ref(1)) * pres_ref(1) / tabs_ref(1)**2
do inum = 1,wtgscale_vertmodenum
  wwtgc(inum) = thetad1 * sin(pi*z(1)*inum/ztrop) * z(1)
end do

do k = 2,ktrop

  thetad1 = (tv_model(k)   - tv_ref(k))   * pres_ref(k)   / tabs_ref(k)**2
  thetad2 = (tv_model(k-1) - tv_ref(k-1)) * pres_ref(k-1) / tabs_ref(k-1)**2
  do inum = 1,wtgscale_vertmodenum
    wwtgc(inum) = wwtgc(inum) + &
                  (thetad1*sin(pi*z(k)*inum/ztrop) + thetad2*sin(pi*z(k)*inum/ztrop)) * &
                  (z(k)-z(k-1))
  end do

end do

wwtgc = wwtgc * ztrop * grav**2 / (lambda_wtg*2)**2 / am_wtg / rgas * -1

do k = 1,ktrop

  w_wtg(k) = 0
  do inum = 1,wtgscale_vertmodenum
    w_wtg(k) = w_wtg(k) + &
               wtgscale_vertmodescl(inum) * wwtgc(inum) * sin(pi*z(k)*inum/ztrop)
  end do

end do

wwtgc(wtgscale_vertmodenum+1) = ztrop

end subroutine wtg_decompdgw
