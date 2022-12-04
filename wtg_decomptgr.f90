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

subroutine wtg_decomptgr(masterproc, nzm, nz, z, &
                          theta_ref, theta_model, tabs_model, tau_wtg, &
                          wtgscale_vertmodenum, wtgscale_vertmodescl, &
                          dowtgLBL, boundstatic, dthetadz_min, &
                          w_wtg, wwtgr, wwtgc)

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
integer, intent(in) :: wtgscale_vertmodenum ! response scaling for 1st vertical mode
real, intent(in) :: wtgscale_vertmodescl(wtgscale_vertmodenum) ! response scaling for 2nd vertical mode

logical, intent(in) :: dowtgLBL    ! Calculate w_wtg at boundary layer instead of linear interpolation
logical, intent(in) :: boundstatic ! Restrict lower bound for static stability
real, intent(in) :: dthetadz_min   ! if boundstatic = .true., what is the minimum bound?

! ======= output =======
real, intent(out) :: w_wtg(nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.
real, intent(out) :: wwtgr(nzm) ! Raw w_wtg, assuming constant = 1
real, intent(out) :: wwtgc(nz) ! Coefficient for the decompositions


! ======= local variables =======
integer :: k
integer :: ktrop
integer :: kbl
integer :: inum
real :: min_temp ! temporary variable used to find cold point of model sounding.
real :: ztrop ! Height of tropopause level (m)
real, parameter :: pi = 3.141592653589793 ! from MATLAB, format long.
real :: thetad1 ! Static Stability
real :: thetad2 ! Static Stability
real :: dthetadz(nzm) ! virtual temperature of model sounding in K

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

! ===== calculate static stability up to ztrop =====

dthetadz(1) = (theta_model(2)-theta_model(1)) / (z(2)-z(1))
if (boundstatic.AND.(dthetadz(1).lt.dthetadz_min).AND.(z(k)>5000)) dthetadz(1) = dthetadz_min
do k = 2,ktrop
  dthetadz(k) = (theta_model(k+1)-theta_model(k-1)) / (z(k+1)-z(k-1))
  if (boundstatic.AND.(dthetadz(k).lt.dthetadz_min).AND.(z(k)>5000)) dthetadz(k) = dthetadz_min
end do

! ===== calculate vertical mode coefficients =====

thetad1 = (theta_model(1) - theta_ref(1)) / dthetadz(1)
do inum = 1,wtgscale_vertmodenum
  wwtgc(inum) = thetad1 * sin(pi*z(1)*inum/ztrop) * z(1)
end do

do k = 2,ktrop

  thetad1 = (theta_model(k)   - theta_ref(k))   / dthetadz(k)
  thetad2 = (theta_model(k-1) - theta_ref(k-1)) / dthetadz(k-1)
  do inum = 1,wtgscale_vertmodenum
    wwtgc(inum) = wwtgc(inum) + &
                  (thetad1*sin(pi*z(k)*inum/ztrop) + thetad2*sin(pi*z(k)*inum/ztrop)) * &
                  (z(k)-z(k-1))
  end do

end do

wwtgc = wwtgc / ztrop * tau_wtg

do k = 1,ktrop

  w_wtg(k) = 0
  do inum = 1,wtgscale_vertmodenum
    w_wtg(k) = w_wtg(k) + &
               wtgscale_vertmodescl(inum) * wwtgc(inum) * sin(pi*z(k)*inum/ztrop)
  end do

end do

if(.NOT.dowtgLBL) then
  do k = kbl,ktrop
    wwtgr(k) = (theta_model(k) - theta_ref(k)) * tau_wtg / dthetadz(k)
  end do
  do k = 1,(kbl-1)
    wwtgr(k) = wwtgr(kbl) * z(k) / z(kbl)
  end do
else
  do k = 1,ktrop
    wwtgr(k) = (theta_model(k) - theta_ref(k)) * tau_wtg / dthetadz(k)
  end do
end if

wwtgc(wtgscale_vertmodenum+1) = ztrop

end subroutine wtg_decomptgr
