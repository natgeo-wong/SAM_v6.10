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
                          theta_ref, theta_model, tabs_model, ttheta_wtg, &
                          ttheta_a, ttheta_ab, &
                          boundstatic, dthetadz_min, w_wtg)

implicit none

! ======= inputs =======
logical, intent(in) :: masterproc ! For printing out to logs
integer, intent(in) :: nzm ! number of model levels
integer, intent(in) :: nz  ! number of vertical levels
real, intent(in) :: z(nz) ! pressure of model levels in Pa (domain-mean for LES)
real, intent(in) :: theta_ref(nzm) ! reference potential temperature sounding in K
real, intent(in) :: theta_model(nzm) ! model potential temperature sounding in K (domain-mean for LES)
real, intent(in) :: tabs_model(nzm) ! model temperature sounding in K (domain-mean for LES)

! WTG potential temperature relaxation timescale (ttheta_wtg)
!   default is 1 day^-1 (ttheta_wtg = 1/86400 s^-1)
real, intent(in) :: ttheta_wtg ! potential temperature relaxation timescale (s^-1)
real, intent(in) :: ttheta_a ! potential temperature relaxation timescale (s^-1)
real, intent(in) :: ttheta_b ! potential temperature relaxation timescale (s^-1)

logical, intent(in) :: boundstatic ! Restrict lower bound for static stability
real, intent(in) :: dthetadz_min   ! if boundstatic = .true., what is the minimum bound?

! ======= output =======
real, intent(out) :: w_wtg(nzm) ! WTG large-scale pressure velocity in Pa/s on model levels.


! ======= local variables =======
integer :: k
integer :: ktrop
integer :: kbl
real :: min_temp ! temporary variable used to find cold point of model sounding.
real :: ztrop ! Height of tropopause level (m)
real :: dthetadz ! Static Stability
real, parameter :: pi = 3.141592653589793 ! from MATLAB, format long.
real :: theta_diff ! Static Stability

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
    kbl = k
  end if
end do

a = (theta_model(kbl)-theta_ref(kbl)) * sin(pi*z(k)  /ztrop) * z(kbl)
b = (theta_model(kbl)-theta_ref(kbl)) * sin(pi*z(k)*2/ztrop) * z(kbl)

do k = (kbl+1),ktrop

  a = a + ((theta_model(k)-theta_ref(k)) * sin(pi*z(k)  /ztrop) - &
           (theta_model(k-1)-theta_ref(k-1)) * sin(pi*z(k-1)/ztrop)) * (z(k)-z(k-1))

  b = b + ((theta_model(k)-theta_ref(k)) * sin(pi*z(k)*2/ztrop) - &
           (theta_model(k-1)-theta_ref(k-1)) * sin(pi*z(k-1)*2/ztrop)) * (z(k)-z(k-1))

end do

a = a / ztrop
b = b / ztrop

do k = kbl,ktrop

  dthetadz = (theta_model(k+1)-theta_model(k-1)) / (z(k+1)-z(k-1))
  if (boundstatic.AND.(dthetadz.lt.dthetadz_min).AND.(z(k)>5000)) dthetadz = dthetadz_min
  ! According to Raymond and Zeng (2005) model feedbacks in the upper troposphere can
  ! otherwise result in very weak static stabilities and unrealistically
  ! large values of w_wtg
  theta_diff = theta_model(k)-theta_ref(k)
  w_wtg(k) = (ttheta_a * a * theta_diff * sin(pi*z(k)  /ztrop) + &
              ttheta_b * b * theta_diff * sin(pi*z(k)*2/ztrop)) * ttheta_wtg / dthetadz

end do

do k = 1,(kbl-1)
  w_wtg(k) = w_wtg(kbl) * z(k) / z(kbl)
end do

end subroutine wtg_decomp2022
