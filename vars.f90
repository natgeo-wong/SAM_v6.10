module vars

use grid

implicit none
!--------------------------------------------------------------------
! prognostic variables:

real u   (dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm) ! x-wind
real v   (dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm) ! y-wind
real w   (dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz ) ! z-wind
real t   (dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! liquid/ice water static energy

!--------------------------------------------------------------------
! diagnostic variables:

real p      (0:nx, (1-YES3D):ny, nzm)     ! perturbation pressure (from Poison eq)
real tabs   (nx, ny, nzm)                 ! temperature
real qv      (nx, ny, nzm)                ! water vapor
real qcl     (nx, ny, nzm)                ! liquid water  (condensate)
real qpl     (nx, ny, nzm)                ! liquid water  (precipitation)
real qci     (nx, ny, nzm)                ! ice water  (condensate)
real qpi     (nx, ny, nzm)                ! ice water  (precipitation)

!--------------------------------------------------------------------
! time-tendencies for prognostic variables

real dudt   (nxp1, ny, nzm, 3)
real dvdt   (nx, nyp1, nzm, 3)
real dwdt   (nx, ny, nz,  3)

!----------------------------------------------------------------
! Temporary storage array:

	real misc(nx, ny, nz)
!------------------------------------------------------------------
! fluxes at the top and bottom of the domain:

real fluxbu (nx, ny), fluxbv (nx, ny), fluxbt (nx, ny)
real fluxbq (nx, ny), fluxtu (nx, ny), fluxtv (nx, ny)
real fluxtt (nx, ny), fluxtq (nx, ny), fzero  (nx, ny)
real precsfc(nx,ny) ! surface precip. rate

! parameters from surface transfer scheme: friction velocity,
!   surface transfer coef for water and surface qsat.  These are
!   useful for computing surface fluxes of isoptopic water.
real ustar (nx, ny), fluxbq_coef (nx, ny), qsat_surf (nx, ny), u10arr(nx, ny)

!-----------------------------------------------------------------
! profiles

real   t0(nzm), q0(nzm), qv0(nzm), tabs0(nzm), tl0(nzm), &
       tv0(nzm), u0(nzm), v0(nzm), &
       tg0(nzm), tp0(nzm), qg0(nzm), ug0(nzm), vg0(nzm), p0(nzm), &
       t01(nzm), q01(nzm), qp0(nzm), qn0(nzm), &
       usounding0(nzm), vsounding0(nzm), tg0_ref(nzm), qg0_ref(nzm), &
       qcldliq0(nzm), qpcpliq0(nzm), qcldice0(nzm), qpcpice0(nzm)
!----------------------------------------------------------------
! "observed" (read from snd file) surface characteristics

real ::  sstobs = 0., lhobs = 0., shobs = 0.
!----------------------------------------------------------------
!  Domain top stuff:

real   gamt0    ! gradient of t() at the top,K/m
real   gamq0    ! gradient of q() at the top,g/g/m

!-----------------------------------------------------------------
! reference vertical profiles:

real   prespot(nzm)  ! (1000./pres)**R/cp
real   prespoti(nzm) ! initial prespot, to convert pottemp soundings to temp for LS forcing
real   rho(nzm)	  ! air density at pressure levels,kg/m3
real   rhow(nz)   ! air density at vertical velocity levels,kg/m3
real   bet(nzm)	  ! = ggr/tv0
real   gamaz(nzm) ! ggr/cp*z
real   wsub(nz)   ! Large-scale subsidence velocity,m/s
real   qtend(nzm) ! Large-scale tendency for total water
real   ttend(nzm) ! Large-scale tendency for temp.

!---------------------------------------------------------------------
! Large-scale and surface forcing:

integer nlsf	! number of large-scale forcing profiles
integer nrfc	! number of radiative forcing profiles
integer nsfc	! number of surface forcing profiles
integer nsnd	! number of observed soundings
integer nzlsf	! number of large-scale forcing profiles
integer nzrfc	! number of radiative forcing profiles
integer nzsnd	! number of observed soundings

real, allocatable :: dqls(:,:) ! Large-scale tendency for total water
real, allocatable :: dtls(:,:) ! Large-scale tendency for temp.
real, allocatable :: ugls(:,:) ! Large-scale wind in X-direction
real, allocatable :: vgls(:,:) ! Large-scale wind in Y-direction
real, allocatable :: wgls(:,:) ! Large-scale subsidence velocity,m/s
real, allocatable :: pres0ls(:)! Surface pressure, mb
real, allocatable :: zls(:,:)  ! Height
real, allocatable :: pls(:,:)  ! Pressure
real, allocatable :: dayls(:)  ! Large-scale forcing arrays time (days)
real, allocatable :: utraj_ls(:)  ! Zonal velocity of forcings (for Lagrangian-type simulations)
real, allocatable :: vtraj_ls(:)  ! Meridional velocity of forcings (for Lagrangian-type simulations)
real, allocatable :: dtrfc(:,:)! Radiative tendency for pot. temp.
real, allocatable :: dayrfc(:) ! Radiative forcing arrays time (days)
real, allocatable :: prfc(:,:) ! Pressure/Height
real, allocatable :: sstsfc(:) ! SSTs
real, allocatable :: shsfc(:)   ! Sensible heat flux,W/m2
real, allocatable :: lhsfc(:)  ! Latent heat flux,W/m2
real, allocatable :: tausfc(:) ! Surface drag,m2/s2
real, allocatable :: daysfc(:) ! Surface forcing arrays time (days)
real, allocatable :: usnd(:,:) ! Observed zonal wind
real, allocatable :: vsnd(:,:) ! Observed meriod wind
real, allocatable :: tsnd(:,:) ! Observed Abs. temperature
real, allocatable :: qsnd(:,:) ! Observed Moisture
real, allocatable :: tsnd_ref(:,:) ! Abs. temperature from reference sounding (from SCAM IOP file)
real, allocatable :: qsnd_ref(:,:) ! moisture from reference sounding (from SCAM IOP file)
real, allocatable :: tsnd_init(:,:) ! Abs. temperature sounding for initialization (from SCAM IOP file)
real, allocatable :: qsnd_init(:,:) ! moisture sounding for initializaton (from SCAM IOP file)
real, allocatable :: tsnd_upwind(:,:) ! Upwind Abs. temperature sounding for relaxational forcing (from SCAM IOP file)
real, allocatable :: qsnd_upwind(:,:) ! Upwind moisture sounding for relaxational forcing (from SCAM IOP file)
real, allocatable :: o3snd_mmr(:,:) ! Ozone sounding (mass mixing ratio) input from SCAM IOP file (optional)
real, allocatable :: zsnd(:,:) ! Height
real, allocatable :: psnd(:,:) ! Pressure
real, allocatable :: daysnd(:) ! number of sounding samples

!bloss(2019-10): Aerosol sounding input through IOP netcdf file
real, allocatable :: AccumAerosolMass_snd(:,:)   ! mass mixing ratio (kg/kg) of accumulation mode aerosol
real, allocatable :: AccumAerosolNumber_snd(:,:) ! number mixing ratio (#/kg) of accumulation mode aerosol
real, allocatable :: AitkenAerosolMass_snd(:,:)   ! mass mixing ratio (kg/kg) of aitken mode aerosol
real, allocatable :: AitkenAerosolNumber_snd(:,:) ! number mixing ratio (#/kg) of aitken mode aerosol

!bloss(2019-10): Isotope sounding input through IOP netcdf file
real, allocatable :: deltaDVapor_snd(:,:)   ! deltaD of water vapor with respect to SMOW (per mil)
real, allocatable :: deltaO18Vapor_snd(:,:)   ! deltaO18 of water vapor with respect to SMOW (per mil)
real, allocatable :: deltaO17Vapor_snd(:,:)   ! deltaO17 of water vapor with respect to SMOW (per mil)

!bloss(2019-10): Reference sounding for water isotopes, input through IOP netcdf file
!   The reference sounding could be used for things like the isotopic composition
!   of moisture converged into the column from nearby columns.
real, allocatable :: deltaDVapor_ref(:,:)   ! deltaD of water vapor with respect to SMOW (per mil)
real, allocatable :: deltaO18Vapor_ref(:,:)   ! deltaO18 of water vapor with respect to SMOW (per mil)
real, allocatable :: deltaO17Vapor_ref(:,:)   ! deltaO17 of water vapor with respect to SMOW (per mil)

!---------------------------------------------------------------------
!  Horizontally varying stuff (as a function of xy)
!
real sstxy(0:nx,(1-YES3D):ny)	!  surface temperature xy-distribution (perturbation from t00)
real fcory(0:ny)      !  Coriolis parameter xy-distribution
real fcorzy(0:ny)      !  z-Coriolis parameter xy-distribution
real latitude(nx,ny)	     ! latitude (degrees)
real longitude(nx,ny)	     ! longitude(degrees)
real prec_xy(nx,ny) ! mean precip. rate for outout
real shf_xy(nx,ny) ! mean precip. rate for outout
real lhf_xy(nx,ny) ! mean precip. rate for outout
real lwns_xy(nx,ny) ! mean net lw at SFC
real swns_xy(nx,ny) ! mean net sw at SFC
real lwnsc_xy(nx,ny) ! clear-sky mean net lw at SFC
real swnsc_xy(nx,ny) ! clear-sky mean net sw at SFC
real lwnt_xy(nx,ny) ! mean net lw at TOA
real swnt_xy(nx,ny) ! mean net sw at TOA
real lwntc_xy(nx,ny) ! clear-sky mean net lw at TOA
real swntc_xy(nx,ny) ! clear-sky mean net sw at TOA
real solin_xy(nx,ny) ! solar TOA insolation
real pw_xy(nx,ny)   ! precipitable water
real cw_xy(nx,ny)   ! cloud water path
real iw_xy(nx,ny)   ! ice water path
real cld_xy(nx,ny)   ! cloud frequency
real u200_xy(nx,ny) ! u-wind at 200 mb
real usfc_xy(nx,ny) ! u-wind at at the surface
real v200_xy(nx,ny) ! v-wind at 200 mb
real vsfc_xy(nx,ny) ! v-wind at the surface
real w500_xy(nx,ny) ! w at 500 mb
real qocean_xy(nx,ny) ! ocean cooling in W/m2

real, parameter :: t00 = 300.   ! constant offset for sstxy

!----------------------------------------------------------------------
!	Vertical profiles of quantities sampled for statitistics purposes:

real &
    twle(nz), twsb(nz), precflux(nz), &
    uwle(nz), uwsb(nz), vwle(nz), vwsb(nz), &
    radlwup(nz), radlwdn(nz), radswup(nz), radswdn(nz), &
    radqrlw(nz), radqrsw(nz), &
    tkeleadv(nz), tkelepress(nz), tkelediss(nz), tkelediff(nz),tkelebuoy(nz), &
    t2leadv(nz),t2legrad(nz),t2lediff(nz),t2leprec(nz),t2lediss(nz), &
    q2leadv(nz),q2legrad(nz),q2lediff(nz),q2leprec(nz),q2lediss(nz), &
    twleadv(nz),twlediff(nz),twlepres(nz),twlebuoy(nz),twleprec(nz), &
    qwleadv(nz),qwlediff(nz),qwlepres(nz),qwlebuoy(nz),qwleprec(nz), &
    momleadv(nz,3),momlepress(nz,3),momlebuoy(nz,3), &
    momlediff(nz,3),tadv(nz),tdiff(nz),tlat(nz), tlatqi(nz),qifall(nz),qpfall(nz)

! scalar statistics output in *.stat file.
!    ISCCP, MODIS, MISR statistics only produced if simulator is enabled in PARAMETERS namelist
real :: &
     w_max = 0., u_max = 0., s_acld = 0., s_acldcold = 0., s_ar = 0., s_arthr = 0., s_sst = 0., &
     ncmn = 0., nrmn = 0., z_inv = 0., z_cb = 0., z_ct = 0., z_cbmn = 0., z_ctmn = 0., &
     z2_inv = 0., z2_cb = 0., z2_ct = 0., cwpmean = 0., cwp2 = 0., &
     precmean = 0., prec2 = 0., precmax = 0., nrainy = 0., ncloudy = 0., &
     s_acldl = 0., s_acldm = 0., s_acldh = 0.,  &
     s_acldisccp = 0., s_acldlisccp = 0., s_acldmisccp = 0., s_acldhisccp = 0., s_ptopisccp = 0., &
     s_acldmodis = 0., s_acldlmodis = 0., s_acldmmodis = 0., s_acldhmodis = 0., s_ptopmodis = 0., &
     s_acldmisr = 0., s_ztopmisr = 0., s_relmodis = 0., s_reimodis = 0., s_lwpmodis = 0., s_iwpmodis = 0., &
     s_tbisccp = 0., s_tbclrisccp = 0., s_acldliqmodis = 0., s_acldicemodis = 0., &
     s_cldtauisccp = 0.,s_cldtaumodis = 0.,s_cldtaulmodis = 0.,s_cldtauimodis = 0.,s_cldalbisccp = 0., &
     s_flns = 0.,s_flnt = 0.,s_flntoa = 0.,s_flnsc = 0.,s_flntoac = 0.,s_flds = 0.,s_fsns = 0., &
     s_fsnt = 0.,s_fsntoa = 0.,s_fsnsc = 0.,s_fsntoac = 0.,s_fsds = 0.,s_solin = 0.

! register functions:


        real, external :: esatw,esati,dtesatw,dtesati,qsatw,qsati,dtqsatw,dtqsati
!        integer, external :: lenstr, bytes_in_rec

! energy conservation diagnostics:

  real(8) total_water_before, total_water_after
  real(8) total_water_evap, total_water_prec, total_water_ls

!===========================================================================
! UW ADDITIONS

! conditional average statistics, subsumes cloud_factor, core_factor, coredn_factor
integer :: ncondavg, icondavg_cld, icondavg_cor, icondavg_cordn, &
     icondavg_satdn, icondavg_satup, icondavg_env
real, allocatable :: condavg_factor(:,:) ! replaces cloud_factor, core_factor
real, allocatable :: condavg_mask(:,:,:,:) ! indicator array for various conditional averages
character(LEN=8), allocatable :: condavgname(:) ! array of short names
character(LEN=25), allocatable :: condavglongname(:) ! array of long names

real   qlsvadv(nzm) ! Large-scale vertical advection tendency for total water
real   tlsvadv(nzm) ! Large-scale vertical advection tendency for temperature
real   ulsvadv(nzm) ! Large-scale vertical advection tendency for zonal velocity
real   vlsvadv(nzm) ! Large-scale vertical advection tendency for meridional velocity

real   qnudge(nzm) ! Nudging of horiz.-averaged total water profile
real   tnudge(nzm) ! Nudging of horiz.-averaged temperature profile
real   unudge(nzm) ! Nudging of horiz.-averaged zonal velocity
real   vnudge(nzm) ! Nudging of horiz.-averaged meridional velocity

real   qstor(nzm) ! Storage of horiz.-averaged total water profile
real   tstor(nzm) ! Storage of horiz.-averaged temperature profile
real   ustor(nzm) ! Storage of horiz.-averaged zonal velocity
real   vstor(nzm) ! Storage of horiz.-averaged meridional velocity

real   utendcor(nzm) ! coriolis acceleration of zonal velocity
real   vtendcor(nzm) ! coriolis acceleration of meridional velocity

! 850 mbar horizontal winds
real u850_xy(nx,ny) ! zonal velocity at 850 mb
real v850_xy(nx,ny) ! meridional velocity at 850 mb

! Surface pressure
real psfc_xy(nx,ny) ! pressure (in millibar) at lowest grid point

! Saturated water vapor path, useful for computing column relative humidity
real swvp_xy(nx,ny)  ! saturated water vapor path (wrt water)

! Cloud and echo top heights, and cloud top temperature (instantaneous)
real cloudtopheight(nx,ny), echotopheight(nx,ny), cloudtoptemp(nx,ny)

! WTG am coefficient, for gradual WTG implementation from RCE to full damping state
! (added by Nathanael Wong on 2021/01/17)
real twtg
real twtgmax
real :: am_wtg_time

logical :: IsInitializedRestartFilename = .false.
CHARACTER(LEN=256) :: RestartFilename, RestartFilenameSave
CHARACTER(LEN=256) :: MiscRestartFilename, MiscRestartFilenameSave

! END UW ADDITIONS
!===========================================================================

end module vars
