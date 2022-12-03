module params

use grid, only: nzm

implicit none

!   Constants:

real, parameter :: cp = 1004.             ! Specific heat of air, J/kg/K
real, parameter :: ggr = 9.81             ! Gravity acceleration, m/s2
real, parameter :: lcond = 2.5104e+06     ! Latent heat of condensation, J/kg
real, parameter :: lfus = 0.3336e+06      ! Latent heat of fusion, J/kg
real, parameter :: lsub = 2.8440e+06      ! Latent heat of sublimation, J/kg
real, parameter :: rv = 461.              ! Gas constant for water vapor, J/kg/K
real, parameter :: rgas = 287.            ! Gas constant for dry air, J/kg/K
real, parameter :: diffelq = 2.21e-05     ! Diffusivity of water vapor, m2/s
real, parameter :: therco = 2.40e-02      ! Thermal conductivity of air, J/m/s/K
real, parameter :: muelq = 1.717e-05      ! Dynamic viscosity of air

real, parameter :: fac_cond = lcond/cp
real, parameter :: fac_fus = lfus/cp
real, parameter :: fac_sub = lsub/cp

real, parameter ::  pi = 3.141592653589793

!
! internally set parameters:

real   epsv     ! = (1-eps)/eps, where eps= Rv/Ra, or =0. if dosmoke=.true.
logical:: dosubsidence = .false.
real fcorz      ! Vertical Coriolis parameter
real coszrs
real salt_factor ! correction factor for water vapor saturation over sea-water

!bloss: options to allow a larger number of substeps and to keep a margin of safety away from the CFL limit
integer:: ncycle_max = 50 ! maximum number of subcycling within dt
integer:: ncycle_min = 1  ! minimum number of subcycling within dt
integer:: ncycle0 = 1     ! number of subcycles for first timestep of simulation.
real :: cfl_safety_factor = 1.5  ! ensure that CFL < CFL_stability_limit / cfl_safety_factor

!----------------------------------------------
! Parameters set by PARAMETERS namelist:
! Initialized to default values.
!----------------------------------------------

real:: ug = 0.        ! Velocity of the Domain's drift in x direction
real:: vg	= 0.        ! Velocity of the Domain's drift in y direction
real:: fcor = -999.   ! Coriolis parameter
real:: longitude0 = 0.    ! latitude of the domain's center
real:: latitude0  = 0.    ! longitude of the domain's center
real:: nxco2 = 1.         ! factor to modify co2 concentration
logical:: doradlat = .false.
logical:: doradlon = .false.

real(8):: tabs_s =0.	! surface temperature,K
real:: delta_sst = 0.   ! amplitude of sin-pattern of sst about tabs_s (ocean_type=1)
real:: depth_slab_ocean = 2. ! thickness of the slab-ocean (m)
real:: Szero = 0.  ! mean ocean transport (W/m2)
real:: deltaS = 0. ! amplitude of linear variation of ocean transport (W/m2)
real:: timesimpleocean = 0. ! time to start simple ocean

real::   fluxt0 =0.  ! surface sensible flux, Km/s
real::   fluxq0 =0.  ! surface latent flux, m/s
real::   tau0   =0.  ! surface stress, m2/s2
real::   z0     =0.035	! roughness length
real::   soil_wetness =1.! wetness coeff for soil (from 0 to 1.)
integer:: ocean_type =0 ! type of SST forcing
logical:: cem =.false.    ! flag for Cloud Ensemble Model
logical:: les =.false.    ! flag for Large-Eddy Simulation
logical:: ocean =.false.  ! flag indicating that surface is water
logical:: land =.false.   ! flag indicating that surface is land
logical:: sfc_flx_fxd =.false. ! surface sensible flux is fixed
logical:: sfc_tau_fxd =.false.! surface drag is fixed

real:: timelargescale =0. ! time to start large-scale forcing

! nudging boundaries (between z1 and z2, where z2 > z1):
real:: nudging_uv_z1 =-1., nudging_uv_z2 = 1000000.
real:: nudging_t_z1 =-1., nudging_t_z2 = 1000000.
real:: nudging_q_z1 =-1., nudging_q_z2 = 1000000.
real:: tauls = 99999999.    ! nudging-to-large-scaler-profile time-scale
real:: tautqls = 99999999.! nudging-to-large-scaler-profile time-scale for scalars

logical:: dodamping = .false.
logical:: doupperbound = .false.
logical:: docloud = .false.
logical:: doprecip = .false.
logical:: dolongwave = .false.
logical:: doshortwave = .false.
logical:: dosgs = .false.
logical:: docoriolis = .false.
logical:: docoriolisz = .false.
logical:: dofplane = .true.
logical:: dosurface = .false.
logical:: dolargescale = .false.
logical:: doradforcing = .false.
logical:: dosfcforcing = .false.
logical:: doradsimple = .false.
logical:: donudging_uv = .false.
logical:: donudging_tq = .false.
logical:: donudging_t = .false.
logical:: donudging_q = .false.
logical:: doensemble = .false.
logical:: dowallx = .false.
logical:: dowally = .false.
logical:: docolumn = .false.
logical:: docup = .false.
logical:: doperpetual = .false.
logical:: doseasons = .false.
logical:: doradhomo = .false.
logical:: dosfchomo = .false.
logical:: dossthomo = .false.
logical:: dodynamicocean = .false.
logical:: dosolarconstant = .false.
logical:: dotracers = .false.
logical:: dosmoke = .false.
logical:: notracegases = .false.
logical:: doseawater = .true. !bloss(Set UW default to true!!) .false.

! Specify solar constant and zenith angle for perpetual insolation.
! Based onn Tompkins and Graig (1998)
! Note that if doperpetual=.true. and dosolarconstant=.false.
! the insolation will be set to the daily-averaged value on day0.
real:: solar_constant = 685. ! solar constant (in W/m2)
real:: zenith_angle = 51.7   ! zenith angle (in degrees)

integer:: nensemble =0   ! the number of subensemble set of perturbations
integer:: perturb_type  = 0 ! type of initial noise in setperturb()

! Initial bubble parameters. Activated when perturb_type = 2
  real:: bubble_x0 = 0.
  real:: bubble_y0 = 0.
  real:: bubble_z0 = 0.
  real:: bubble_radius_hor = 0.
  real:: bubble_radius_ver = 0.
  real:: bubble_dtemp = 0.
  real:: bubble_dq = 0.

! Option for simple treatment of aerosol radiative effect in RRTMG by Tak Yamaguchi
logical:: doradaerosimple = .false.

! Option for writing restart file only at end of simulation
logical:: dorestart_last_only = .false.

!!! KUANG_LAB OPTIONS
logical:: dompiensemble = .false. ! Subdomains defined in domains.f90 are run separately

! Damped Gravity Wave and Temperature Gradient Relaxation Implementations
logical :: dodgw = .false.
logical :: dotgr = .false.
logical :: dowtg_decomp = .false.
real :: wtgscale_time = 0. ! period over which theta relaxation timescale scales from infinity to ttheta_wtg.  Express as fraction of time over which WTG large-scale forcing is implemented.  So if WTG/Large-scale is turned on for 100 days, twtg_scale = 1/4 means that the scaling up to WTG occurs over 25 days.

logical :: dowtg_raymondzeng_QJRMS2005 = .false.
logical :: dowtg_daleuetal_JAMES2015 = .false.
logical :: dowtg_decompdgw = .false.
logical :: dowtg_decomptgr = .false.

real :: tau_wtg = 1. ! Relaxation timescale (in hours) for WTG Approximation of Raymond and Zeng [2005]
logical :: dowtgLBL = .false.
logical :: boundstatic = .true. ! Restrict the static stability lower bound to prevent unrealistically large values of w_wtg
real :: dthetadz_min = 1.e-3 ! if boundstatic = .true., what is the minimum bound? Default from Raymond & Zeng [2005] is 1.e-3 K/km

integer :: wtgscale_vertmodenum = 2 ! adjustment coefficient for first vertical mode
real, dimension(wtgscale_vertmodenum) :: wtgscale_vertmodescl = (/1., 1./) ! adjustment coefficient for second vertical mode

logical :: doradtendency = .false. ! Radiative tendencies as per Pauluis & Garner [2006]
real :: troptend = 1.5 ! Convective tendency in Pauluis & Garner [2006]

logical :: dooceantimeperturb = .false.
real, dimension(5) :: tabs_ptscale = (/0., 0., 0., 0., 0./) ! Vector of sinusoidal periods, units in days
real, dimension(5) :: tabs_pamp    = (/0., 0., 0., 0., 0./) ! Vector of sinusoidal amplitudes, units in K
real, dimension(5) :: tabs_pphase  = (/0., 0., 0., 0., 0./) ! Vector of sinusoidal amplitudes, units in K

end module params
