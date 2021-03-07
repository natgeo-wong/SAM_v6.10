module microphysics

! main interface to Morrison microphysics.
! original implementation by Peter Blossey, UW

use params, only: rgas, rv, cp, lcond, lsub, fac_cond, fac_sub, ggr, &
     ug, vg, & ! domain advection velocities for ship plume.
     pi, doprecip, docloud

use src_scavenging, only: memory, init_scavenging, m2011_scavenging, scav_cloud_2m, qxeps
use aerosol_utils, only: DryAerosolMassFraction, PartitionAerosolMass
use mirage_wateruptake, only: modal_aero_kohler

use grid, only: nx,ny,nzm,nz, &  !grid dimensions; nzm = nz-1 # of scalar lvls
     dimx1_s,dimx2_s,dimy1_s,dimy2_s, & ! actual scalar-array dimensions in x,y
     z, zi, pres, adz, dz, dx, nx_gl, dy, ny_gl, &
     time, dt, dtfactor, ncycle, nstat, nstatfrq, nrestart, day, &
     rank, dompi, masterproc, nsubdomains, nsave3D, save3Dbin, &
     dostatis, doSAMconditionals, dosatupdnconditionals, &
     case, caseid, &
     save2Dbin, save2Davg, nsave2D, &
     compute_reffc, compute_reffi, compute_reffl, &
     do_chunked_energy_budgets, nsaveMSE

use vars, only: pres, rho, dtn, w, t, tabs, qv, qcl, qpl, qci, qpi, &
     gamaz, rhow, tabs0, t0, q0, qv0, &
     fluxbq, fluxtq, u10arr, precsfc, prec_xy, &
     precflux, qpfall, tlat, tlatqi, &
     nrainy, nrmn, ncmn, &
     condavg_mask, ncondavg, condavgname, condavglongname, &
     nstep, nstatis, nprint, icycle, total_water_prec, &
     AccumAerosolMass_snd, AccumAerosolNumber_snd, &
     nsnd,nzsnd,daysnd,zsnd,psnd
     
use module_mp_GRAUPEL, only: GRAUPEL_INIT, M2005MICRO_GRAUPEL, polysvp
use micro_params
use radar_simulator_types, only: class_param

implicit none

logical :: isallocatedMICRO = .false., isallocatedSCAV3D = .false., isallocatedMKBUDGET = .false.

integer :: nmicro_fields ! total number of prognostic water vars

real, allocatable, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
!bloss: Move all other indices to micro_params.f90 so that they can be used inside module_mp_graupel.f90
integer :: index_water_vapor ! separate water vapor index used by SAM

integer :: loc_nan(1)=0
integer :: glb_nan(1)=0

real, allocatable, dimension(:) :: lfac
integer, allocatable, dimension(:) :: flag_wmass, flag_precip, flag_number, flag_advect
integer, allocatable, dimension(:) :: flag_micro3Dout

! number of fields output from micro_write_fields2D, micro_write_fields3D
integer :: nfields2D_micro=0 
integer :: nfields3D_micro=0 

integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

real, allocatable, dimension(:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
real, allocatable, dimension(:,:,:) :: reffc, reffi, reffs, reffr
real, allocatable, dimension(:,:,:) :: &
     CloudLiquidMassMixingRatio, CloudLiquidGammaExponent, CloudLiquidLambda, &
     CloudIceMassMixingRatio, SnowMassMixingRatio
     
!bloss(2021-01): dummy arrays that allow multi-category P3 ice to work with RAD_RRTM/
integer :: nCat_ice_P3
real, allocatable, dimension(:,:,:,:) :: IceMassMixingRatio_P3, ReffIce_P3

real, allocatable, dimension(:,:,:) :: cloudliq

real, allocatable, dimension(:,:,:,:) :: mtend3d !brnr for 3d microtendency output
! variables accumulating precipitation and sedimentation tendencies for use in mse.f90
real, allocatable, dimension(:,:), public :: prec_accum, prec_ice_accum
real, allocatable, dimension(:,:,:), public :: qtot_sed, qice_sed

!bloss(2020-10): Make it possible to output budgets for specific
!  microphysical variables.  These are averaged over subdomains using the
!  framework in SRC/mse.f90 (enabled with do_chunked_energy_budgets=.true.
!  in the PARAMETERS namelist and do_vapor_chunk_budgets=.true. in the 
!  MICRO_THOMPSON namelist) and provide a view of how microphysical and
!  energy budgets vary across the domain and in time.
!   -  n_mkbudget is the number of budgets to be computed based on
!       individual microphysical variables in micro_field or the sum of
!       different micro_field variables.
!   -  n_mkbudget_extra allows individual microphysical or other tendencies
!       such as autoconversion rates, rain evaporation or aerosol scavenging.
!
logical :: do_chunk_mkbudget = .false., do_mkbudget_extra = .false.
integer :: n_mkbudget = 0, n_mkbudget_extra = 1
integer, allocatable, dimension(:,:) :: flag_mkbudget ! which fields should be included in each chunk budget output in mse.f90
real, allocatable, dimension(:,:,:,:) :: mkbudget_sed, mkbudget_extra
character*8, allocatable, dimension(:) :: mkbudget_name, mkbudget_extra_name
character*80, allocatable, dimension(:) :: mkbudget_longname, mkbudget_extra_longname
character*10, allocatable, dimension(:) :: mkbudget_units, mkbudget_extra_units

!bloss(2020-10): Add 3d outputs of scavenging tendencies
logical, save :: isinitializedSCAV = .false.
integer, save :: iscv = -1
real, allocatable, dimension(:,:,:,:) :: scav3d ! last index: qadr, qadcl, nadr, nadcl
character*8, allocatable, dimension(:) :: scavname
character*80, allocatable, dimension(:) :: scavlongname
character*10, allocatable, dimension(:) :: scavunits


real, allocatable, dimension(:,:) :: & ! statistical arrays
     mkwle, & ! resolved vertical flux
     mkwsb, & ! SGS vertical flux
     mksed, & ! sedimentation vertical flux
     mkadv, & ! tendency due to vertical advection
     mkdiff, &! tendency due to vertical diffusion
     mklsadv, & ! tendency due to large-scale vertical advection
     mkstor, & ! storage term BRNR 7/12/2012
     mfrac, & ! fraction of domain with microphysical quantity > 1.e-6
     stend, & ! tendency due to sedimentation
     mtend, & ! tendency due to microphysical processes (other than sedimentation)
     trtau, & ! optical depths of various species
     micro_proc_rates, & ! process rates of individual microphysical interactions
     mk0, &
     mk_ref

logical, allocatable, dimension(:) :: is_water_vapor

real, allocatable, dimension(:) :: tmtend

real :: sfcpcp, sfcicepcp

! arrays with names/units for microphysical outputs in statistics.
character*3, allocatable, dimension(:) :: mkname
character*80, allocatable, dimension(:) :: mklongname
character*10, allocatable, dimension(:) :: mkunits
real, allocatable, dimension(:) :: mkoutputscale

!options for coupling cloud radiative properties to information
!  from the microphysics
logical :: douse_reffc = .true., douse_reffi = .true., douse_reffl = .false.
logical :: dosnow_radiatively_active = .true.
logical :: dorrtm_cloud_optics_from_effrad_LegacyOption = .false.

! You can also have some additional, diagnostic, arrays, for example, total
! nonprecipitating cloud water, etc:

!bloss: array which holds temperature tendency due to microphysics
real, allocatable, dimension(:,:,:), SAVE :: tmtend3d

!brnr (Jun 11): add option for prognostic aerosol
!bloss (Feb 2018): Most aerosol options now sit in micro_params.f90.
!logical :: doprogaertest deprecate this option and use aerinitmode instead
real :: dum, dum2, dum3, dum4, dum5, dumM1, dumM2

!bloss(2019-07): Statistical outputs for scavenging tendencies
real, allocatable, dimension(:) :: scvtndqadclstat,scvtndqadrstat, &
     scvtndnadclstat,scvtndnadrstat

! array variable for aerosol sigma_g
real, dimension(1) :: aer_sig_arr

!brnr (Feb 13): option for shiptrack in 2D simulations
logical :: shiptrack_timeflag !flag for if shiptrack has been created

real, dimension(:,:), allocatable :: aer_nflux_avg_xy, aer_qflux_avg_xy

!bloss (Apr 09): Add option for output of cloud radar reflectivity.
!                Computed using quickbeam cloud radar simulator.
!                Will be output as histogram in statistics file 
!                (with nradar_bins bins between -40 and 20 dBZ) and
!                in 3D files as a full, instantaneous 3D field.
integer :: nbins_cloudradar
real, allocatable, dimension(:) :: binedges_cloudradar
real, allocatable, dimension(:,:) :: hist_cloudradar
real, allocatable, dimension(:,:,:) :: dBZ_cloudradar
character(LEN=8), allocatable, dimension(:) :: binname_cloudradar
type(class_param) :: hp_cloudradar ! hydrometeor class settings
integer :: nhclass

integer :: nsizes_cloudradar
logical :: dostatis_quickbeam
real :: factor_quickbeam

logical :: ldebug = .false.

CONTAINS
!----------------------------------------------------------------------
function micro_scheme_name()
  character(len=32) :: micro_scheme_name
  ! Return the scheme name, normally the same as the directory name with leading "MICRO_" removed  
  micro_scheme_name = "m2005_pa" 
end function   micro_scheme_name
!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
subroutine micro_setparm()
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder
  integer :: count
  
   NAMELIST /MICRO_M2005_PA/ &
      dototalwater, &       ! use total water variable (vapor + cloud liquid)
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! graupel species has qualities of hail
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization in place of KK(2000)
      dopredictNc, &        ! prediction of cloud droplet number
      doprogaerosol, &      ! use Berner et al. (2011) prognostic aerosol !brnr
      !doprogaertest, &     ! initialize test aerosol distribution for prognostic aerosol (deprecated)
      aerinitmode, &        ! integer to set type of initial aerosol distribution
      aerfluxmode, &        ! integer to set type of surface flux for aerosol
      fluxQAd0, &           ! value for fixed surface aerosol mass flux
      fluxNAd0, &           ! value for fixed surface aerosol number flux
      sflux_nacc_coef, &    ! coefficient for surface aerosol number flux
      sflux_racc, &         ! characteristic radius of aerosol produced by surface flux (microns)
      whitecap_coef, &      ! white cap coefficient in Clarke et al (2006) surface aerosol flux formula
      do_m2011_scavenge, &  ! use Muhlbauer M2011 interstitial scavenging scheme
      doscavcloud2m, &      ! use 2 moment version of cloud scavenging
      do_scav_3d_output, &
      dodissip3d,    &      ! use 3d turbulent dissipation rate for cloud scavenging
      nscv, &               ! number of time steps between scavenging calculations
      doprecoff, &          ! turn off autoconversion in KK warm rain !brnr
      dosedoff, &           ! turn off sedimentation !brnr
      doevapnc, &           ! allow evaporation to remove nc !brnr
      dospecifyaerosol, &   ! specify two modes of (sulfate) aerosol
      doshiptrack2D, &      ! pulse of aerosol for shiptrack simulation
      aer_ship_n, &         ! number concentration in shiptrack pulse
      shiptrack_time0, &    ! initial time in seconds for pulse release
      doshipv2, &             ! model ship emissions as surface flux 
      shipv2_radius, &         ! Width of shiptrack plume at surface 
      shipv2_aer_n, &         ! number concentration in shiptrack plume
      shipv2_aer_r, &         ! number concentration in shiptrack plume
      shipv2_time0, &         ! initial time in seconds for plume release
      shipv2_timef, &         ! final time in seconds for plume release
      shipv2_mode, &      ! type of ship track (0=stationary ship once through domain)
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for dospecifyaerosol=.false. (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for dospecifyaer...=.true.
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      doPOCinit_Na, &       ! perturb aerosol concentration to try and force POC formation
      doPOCinit_qt, &       ! perturb qt to try and force POC formation
      POCperturb_Na, &     ! amount by which to perturb POC region Na (#/mg)
      POCperturb_qt, &     ! amount by which to perturb POC region Qt (g/kg)
      POC_trans_length, &   ! length in meters of transition between POC and OVC
      OVC_length, &         ! desired portion of domain set to be initially OVC
      dofix_pgam, pgam_fixed, & ! option to specify pgam (exponent of cloud water's gamma distn)
      douse_reffc, &        ! use computed effective radius in radiation computation
      douse_reffl, &        ! use computed effective radius for cloud and drizzle in radiation (Brnr)
      douse_reffi, &        ! use computed effective ice size in radiation computation
      massfactor, &         ! factor by which to rescale the aerosol mass field
      dorrtm_cloud_optics_from_effrad_LegacyOption, & 
      dosnow_radiatively_active, &
      do_scale_dependence_of_autoconv, &  ! allow heuristic scaling based on dx
      do_scale_dependence_of_activation, &! both default to true.
      do_output_micro_process_rates, &
      do_aerosol_chunk_budgets, &
      doreflectivity_cloudradar, & ! use quickbeam cloud radar simulator
      binwidth_cloudradar, & ! histogram bin width in dBZ (from -40 to 20 dBZ, default=5dBZ)
      min_dBZbin_cloudradar, & ! min end of reflectivity histogram in dBZ (default=0 dBZ)
      max_dBZbin_cloudradar, & ! max end of reflectivity histogram in dBZ (default=70 dBZ)
      freq_cloudradar, & ! frequency of cloud radar (default = 94 GHz for cloudsat)
      surface_cloudradar, & ! location of cloud radar (1=surface, 0=space)
      usegasabs_cloudradar, & ! include gas absorption in reflectivity computations (1=yes, 0=no)
      doray_cloudradar, & ! do Rayleigh computations of reflectivity for comparison (1=yes, 0=no)
      InCloudSuperSatPercentForScavenging, & ! Fix in-cloud supersaturation for scavenging code
      hygro, & ! hygroscopicity of aerosol
      DoScavengeAerosolDryDiameter ! Testing option to evaluate how use of wet diameter affects scavenging.  Do not use in science runs.
!bloss: add later??
!bloss      do_micro3dout         ! flag to output 3d snapshot for a bunch of micro tendencies


   !bloss: Create dummy namelist, so that we can figure out error code
   !       for a mising namelist.  This lets us differentiate between
   !       missing namelists and those with an error within the namelist.
   NAMELIST /BNCUIODSBJCB/ place_holder

   !bloss(2015-02): Default values for namelist variables moved to micro_params.f90

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
  
  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_M2005_PA namelist
  read (55,MICRO_M2005_PA,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_M2005_PA namelist'
        rewind(55)
        read (55,MICRO_M2005_PA) ! this should give a useful error message
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '****** No MICRO_M2005_PA namelist in prm file *********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

   if(.not.doicemicro) dograupel=.false.

  if(doprogaerosol) then
     ! make some consistency checks for prognostic aerosol setup

     if(.NOT.dopredictNC) then
        if(masterproc) write(*,*) 'dopredictnc must be .true. for doprogaerosol to be used.'
        call task_abort()
     end if

     if((aerfluxmode.eq.2).AND.(sflux_nacc_coef.eq.0.)) then
        if(masterproc) then 
           write(*,*) 'When using aerfluxmode==2, use sflux_nacc_coef to set'
           write(*,*) '  the coefficient of the surface aerosol number flux and'
           write(*,*) '  sflux_racc to set the characteristic radius of the flux.'
           write(*,*) 'If you really want to run with sflux_nacc_coef==0.,'
           write(*,*) '  disable this check in SRC/MICRO_M2005_PA/microphysics.f90'
        end if
        call task_abort()
     end if

     nfields2d_micro = 12

  end if

   ! write namelist values out to file for documentation
   if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.nml', form='formatted', position='append')    
      write (unit=55,nml=MICRO_M2005_PA,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if

   ! scale values of parameters for m2005micro
   aer_rm1 = 1.e-6*aer_rm1 ! convert from um to m
   aer_rm2 = 1.e-6*aer_rm2 
   aer_n1 = 1.e6*aer_n1 ! convert from #/cm3 to #/m3
   aer_n2 = 1.e6*aer_n2
   aer_ship_n = 1.e6*aer_ship_n 
  
  nmicro_fields = 1 ! start with water vapor and cloud water mass mixing ratio
  if(docloud) then
    if(.NOT.dototalwater) nmicro_fields = nmicro_fields + 1 ! add cloud water mixing ratio
     if(dopredictNc) nmicro_fields = nmicro_fields + 1 ! add cloud water number concentration (if desired)
     if(doprogaerosol) then 
!bloss(2020-11):        nmicro_fields = nmicro_fields + 3 !brnr add single wet and dry prognostic aerosol modes with 2 moments
       !bloss: prognostic wet+dry aerosol number and mass,
       !       diagnostic dry aerosol number and mass, wet (in-cloud-droplet) aerosol mass
        nmicro_fields = nmicro_fields + 5
        dospecifyaerosol = .true.
        if(.NOT.doprecoff.AND.doprecip) then 
           nmicro_fields = nmicro_fields + 1 !brnr add field for aerosol mass in rain
        end if
     end if
  end if
  if(doprecip)    nmicro_fields = nmicro_fields + 2 ! add rain mass and number (if desired)
  if(doicemicro)  nmicro_fields = nmicro_fields + 4 ! add snow and cloud ice number and mass (if desired)
  if(dograupel)   nmicro_fields = nmicro_fields + 2 ! add graupel mass and number (if desired)

  ! specify index of various quantities in micro_field array
  !  *** note that not all of these may be used if(.not.doicemicro) ***
  if(dototalwater) then
    iqv = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
    !bloss/qt  iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
    !bloss/qt: cloud liquid water not prognosed
    count = 1
  else
    iqv = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
    iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
    count = 2
  end if
  
  if(dopredictNc) then
    incl = count + 1  ! cloud water number mixing ratio [#/kg dry air]
    count = count + 1
  end if

  if(doprecip) then
    iqr = count + 1 ! rain mass mixing ratio [kg H2O / kg dry air]
    inr = count + 2 ! rain number mixing ratio [#/kg dry air]
    count = count + 2

    if(doprogaerosol) then
      iqar = count + 1 ! ! wet aerosol mass mixing ratio in rain [kg aerosol in rain/kg dry air]
      count = count + 1
    end if
  end if

  if(doicemicro) then
     iqci = count + 1  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = count + 2  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = count + 3   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = count + 4   ! snow number mixing ratio [#/kg dry air]
     count = count + 4
   end if

   if(dograupel) then
     iqg = count + 1   ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = count + 2  ! graupel number mixing ratio [#/kg dry air]
     count = count + 2
  end if

  !bloss(2020-11): Move prognostic aerosol to end of micro_field()
  if(doprogaerosol) then
    ! NOTE: prognostic variables are dry+wet aerosol mass and number
    !  The mass mixing ratios of dry and wet aerosol and dry aerosol number
    !    are diagnostic variables.
    !  Here wet aerosol means aerosol associated within cloud droplets
    !    while dry aerosol is in cloud-free air or interstitial aerosols
    !    within the cloud layer.
    inacc = count + 1 ! dry + wet aerosol mass mixing ratio [kg aerosol/kg dry air]
    iqacc = count + 2 ! dry aerosol mass mixing ratio [kg aerosol/kg dry air]
    ! diagnostic variables
    inad = count + 3 ! dry aerosol number mixing ratio [#/kg dry air]
    iqad = count + 4 ! dry aerosol mass mixing ratio [kg aerosol/kg dry air]
    iqaw = count + 5 ! wet aerosol mass mixing ratio [kg activated aerosol/kg dry air]
    count = count + 5
  end if

  if(count.ne.nmicro_fields) then
    if(masterproc) write(*,*) 'Error in MICRO_M2005_PA, micro_setparm():'
    if(masterproc) write(*,*) 'Mismatch between number of microphysical species and indices to specific species.'
    call task_abort()
  end if

  ! stop if icemicro is specified without precip -- we do not support this right now.
  if((doicemicro).and.(.not.doprecip)) then
     if(masterproc) write(*,*) 'Morrison 2005 Microphysics does not support both doice and .not.doprecip'
     call task_abort()
  end if
  index_water_vapor = iqv ! set SAM water vapor flag

  if(do_output_micro_process_rates) then
    ! set up number of process rate outputs from module_mp_graupel
    nmicro_proc = nmicro_process_rates_warm_mass &
         + nmicro_process_rates_warm_number
    if(doprogaerosol) then
      nmicro_proc = nmicro_proc + nmicro_process_rates_progaer_mass &
           + nmicro_process_rates_progaer_number
    end if
    if(doicemicro) then
      nmicro_proc = nmicro_proc + nmicro_process_rates_cold_mass &
           + nmicro_process_rates_cold_number
    end if
  else
    nmicro_proc = 1 ! default to one if no outputs
  end if

  if(do_chunked_energy_budgets.AND.do_aerosol_chunk_budgets.AND.doprogaerosol.AND.(.NOT.isallocatedMKBUDGET)) then
    ! using the mse.f90 routine, output chunk averaged budgets for total accumulation mode aerosol number and mass
    do_chunk_mkbudget = .true.
    n_mkbudget = 2 + 2 ! total aerosol number, total aerosol mass, water vapor, rain
    if(.NOT.dototalwater) n_mkbudget = n_mkbudget + 1 ! add cloud liquid
    do_mkbudget_extra = .true.
    n_mkbudget_extra =  7 ! autoconversion (3), accretion (2), rain evaporation (2)
    if(do_m2011_scavenge) n_mkbudget_extra =  n_mkbudget_extra + 4
    allocate(flag_mkbudget(nmicro_fields,n_mkbudget), &
         mkbudget_sed(nx,ny,nzm,n_mkbudget), mkbudget_extra(nx,ny,nzm,n_mkbudget_extra), &
         mkbudget_name(n_mkbudget), mkbudget_extra_name(n_mkbudget_extra), &
         mkbudget_longname(n_mkbudget), mkbudget_extra_longname(n_mkbudget_extra), &
         mkbudget_units(n_mkbudget), mkbudget_extra_units(n_mkbudget_extra), &
         STAT=ierr)
    if(ierr.ne.0) STOP 'when allocating mkbudget arrays in MICRO_M2005_PA/microphysics.f90'
    flag_mkbudget(:,:) = 0
    ! first budget is for total aerosol number = NAc+NR = NAd + NC + NR
    flag_mkbudget(inacc,1) = 1
    flag_mkbudget(inr,1) = 1
    ! second budget is for total aerosol mass = QAc+QAr = QAd+QAw+QAr
    flag_mkbudget(iqacc,2) = 1
    flag_mkbudget(iqar,2) = 1
    ! third budget is for water vapor (or qv+qcl if dototalwater=.true.)
    flag_mkbudget(index_water_vapor,3) = 1    
    ! fourth budget is for rain
    flag_mkbudget(iqr,4) = 1    
    if(.NOT.dototalwater) then
      ! fifth budget is for cloud liquid
      flag_mkbudget(iqcl,5) = 1    
    end if

    mkbudget_sed(:,:,:,:) = 0.
    isallocatedMKBUDGET = .true.
  end if

  if(.not.isallocatedMICRO) then
     ! allocate microphysical variables
     allocate(micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,nmicro_fields), &
          fluxbmk(nx,ny,nmicro_fields), fluxtmk(nx,ny,nmicro_fields), &
          reffc(nx,ny,nzm), reffr(nx,ny,nzm), reffi(nx,ny,nzm), reffs(nx,ny,nzm), &
          CloudLiquidMassMixingRatio(nx,ny,nzm), CloudLiquidGammaExponent(nx,ny,nzm), &
          CloudLiquidLambda(nx,ny,nzm), &
          CloudIceMassMixingRatio(nx,ny,nzm), SnowMassMixingRatio(nx,ny,nzm), &
          mkwle(nz,nmicro_fields), mkwsb(nz,nmicro_fields), &
          mkadv(nz,nmicro_fields), mkdiff(nz,nmicro_fields), &
          mklsadv(nz,nmicro_fields), mkstor(nz,nmicro_fields), &
          stend(nzm,nmicro_fields), mtend(nzm,nmicro_fields), &
          mfrac(nzm,nmicro_fields), trtau(nzm,nmicro_fields), &
          micro_proc_rates(nzm,nmicro_proc), &
          mksed(nzm,nmicro_fields), tmtend(nzm), &
          cloudliq(nx,ny,nzm), &
          qtot_sed(nx,ny,nzm), qice_sed(nx,ny,nzm), & ! for budgets in mse.f90
          prec_accum(nx,ny), prec_ice_accum(nx,ny), & ! for budgets in mse.f90
          aer_nflux_avg_xy(nx,ny), aer_qflux_avg_xy(nx,ny), & ! for 2d progaer outputs
          tmtend3d(nx,ny,nzm), flag_micro3Dout(nmicro_fields), &
          flag_wmass(nmicro_fields), flag_precip(nmicro_fields), &
          flag_advect(nmicro_fields), &
          flag_number(nmicro_fields), lfac(nmicro_fields), &
          mkname(nmicro_fields), mklongname(nmicro_fields), &
          mkunits(nmicro_fields), mkoutputscale(nmicro_fields), &
          mk0(nzm,nmicro_fields), is_water_vapor(nmicro_fields), &
          scvtndqadclstat(nzm),scvtndqadrstat(nzm), &
          scvtndnadclstat(nzm),scvtndnadrstat(nzm), &
          STAT=ierr)
     if(ierr.ne.0) then
        write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
        call task_abort()
     else
        isallocatedMICRO = .true.
     end if

     ! zero out statistics variables associated with cloud ice sedimentation
     !   in Marat's default SAM microphysics
     tlatqi = 0.

     ! initialize these arrays
     micro_field = 0.
     cloudliq = 0. !bloss/qt: auxially cloud liquid water variable, analogous to qn in MICRO_SAM1MOM
     fluxbmk = 0.
     fluxtmk = 0.
     mkwle = 0.
     mkwsb = 0.
     mkadv = 0.
     mkdiff = 0.
     mklsadv = 0.
     mkstor = 0.
     mk0 = 0. !bloss: placeholder

     aer_nflux_avg_xy(:,:) = 0.
     aer_qflux_avg_xy(:,:) = 0.

    ! initialize flag arrays to all mass, no number, no precip
     flag_wmass = 1
     flag_number = 0
     flag_precip = 0
     flag_advect = 1
     flag_micro3Dout = 0

     ! by default, effective radii in microphysics will be used in radiation,
     !   though this can be changed in the namelist using douse_reff*
  compute_reffc = douse_reffc
  compute_reffl = douse_reffl
  compute_reffi = douse_reffi

     if(dorrtm_cloud_optics_from_effrad_LegacyOption) then
       !bloss(2016-02-09): If using legacy radiative treatment, make sure snow is
       !   not radiatively active.
       dosnow_radiatively_active = .false.
       if(masterproc) write(*,*) '*** Snow is not radiatively active when using legacy radiation ***'
     end if

     ! initialize fields useful for radiation
     reffc = 25.
     reffi = 25.

     CloudLiquidMassMixingRatio = 0.
     CloudLiquidGammaExponent = 0.
     CloudLiquidLambda = 0.
     CloudIceMassMixingRatio = 0.
     SnowMassMixingRatio = 0.

     is_water_vapor(:) = .false.
     is_water_vapor(iqv) = .true.

     ! set up stuff for cloud radar simulator output (uses QUICKBEAM)
     if(doreflectivity_cloudradar) then
       call cloudradar_init( )
     end if

  end if

end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the 
!   beginning of each run, initial or restart:
subroutine micro_init()

  implicit none
  
  real, dimension(nzm) :: qc0, rh0, Na_accum, qa_accum
  real :: tmp_pgam, tmp_lambda

  real, external :: satadj_water, qsatw

  !berner variables
  real :: tmpRH, tmpx1,tmpx2, POC_length, arg, tmpNa, tmpqa, tmpqv
  integer :: k, m, n, it, jt, nx_Nc1, nx_Nc2, nx_trans, tmp_ind

  integer :: i, j, kinv, ierr
  real :: tmp_max, tmp_check, dry_ratio

  real :: InitialMode1MeanMass, InitialMode2MeanMass

  !bloss/qt: with the new dototalwater option, fill in the flag arrays element by element.
  flag_wmass(:) = 0
  flag_precip(:) = 0
  flag_number(:) = 0

  flag_wmass(iqv) = 1

  if(.NOT.dototalwater) then
    flag_wmass(iqcl) = 1 ! liquid water mass
  end if

  if(dopredictNc) then
!bloss(flag_#): flag_number only affect the scaling of the 3D output by rho
!bloss            so that it is not necessary to use it here since we want
!bloss            the output as a mixing ratio
!bloss(flag_#)    flag_number(incl) = 1 ! liquid water number
  end if

  if(doprecip) then
    flag_wmass(iqr) = 1 ! rain mass
    flag_precip(iqr) = 1 ! rain as precip
!bloss(flag_#)    flag_number(inr) = 1 ! rain number
    flag_precip(inr) = 1 ! rain number as precip
  end if

  if(doicemicro) then
    flag_wmass(iqci) = 1 ! cloud ice mass
!bloss(flag_#)    flag_number(inci) = 1 ! cloud ice number

    flag_wmass(iqs) = 1 ! snow mass
    flag_precip(iqs) = 1 ! snow as precip
!bloss(flag_#)    flag_number(ins) = 1 ! snow number
    flag_precip(ins) = 1 ! snow number as precip

    if(dograupel) then
      flag_wmass(iqg) = 1 ! graupel mass
      flag_precip(iqg) = 1 ! graupel as precip
!bloss(flag_#)      flag_number(ing) = 1 ! graupel number
      flag_precip(ing) = 1 ! graupel number as precip
    end if
  end if

  if(doprogaerosol) then
!bloss(flag_#)    flag_number(inacc) = 1 ! dry+cloud aerosol number mixing ratio
!bloss(flag_#)    flag_number(inad) = 1 ! dry aerosol number mixing ratio

    ! do not advect and diffuse these diagnostic variables.
    flag_advect(inad) = 0
    flag_advect(iqad) = 0
    flag_advect(iqaw) = 0
  end if

!!$  ! initialize flag arrays
!!$  if(dopredictNc) then
!!$     ! Cloud droplet number concentration is a prognostic variable
!!$     if(doicemicro) then
!!$        if(dograupel) then
!!$          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns, qg, Ng
!!$           flag_wmass  = (/1,0,1,0,1,0,1,0,1,0/)
!!$           flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
!!$           flag_number = (/0,1,0,1,0,1,0,1,0,1/)
!!$        else
!!$          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns
!!$           flag_wmass  = (/1,0,1,0,1,0,1,0/)
!!$           flag_precip = (/0,0,1,1,0,0,1,1/)
!!$           flag_number = (/0,1,0,1,0,1,0,1/)
!!$        end if
!!$     else
!!$        if(doprecip) then
!!$          !bloss/qt: qt, Nc, qr, Nr
!!$           flag_wmass  = (/1,0,1,0/)
!!$           flag_precip = (/0,0,1,1/)
!!$           flag_number = (/0,1,0,1/)
!!$        else
!!$          !bloss/qt: qt, Nc
!!$           flag_wmass  = (/1,0/)
!!$           flag_precip = (/0,0/)
!!$           flag_number = (/0,1/)
!!$        end if
!!$     end if
!!$  else
!!$     ! Cloud droplet number concentration is NOT a prognostic variable
!!$     if(doicemicro) then
!!$        if(dograupel) then
!!$          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns, qg, Ng
!!$           flag_wmass  = (/1,1,0,1,0,1,0,1,0/)
!!$           flag_precip = (/0,1,1,0,0,1,1,1,1/)
!!$           flag_number = (/0,0,1,0,1,0,1,0,1/)
!!$        else
!!$          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns
!!$           flag_wmass  = (/1,1,0,1,0,1,0/)
!!$           flag_precip = (/0,1,1,0,0,1,1/)
!!$           flag_number = (/0,0,1,0,1,0,1/)
!!$        end if
!!$     else
!!$        if(doprecip) then
!!$          !bloss/qt: qt, qr, Nr
!!$           flag_wmass  = (/1,1,0/)
!!$           flag_precip = (/0,1,1/)
!!$           flag_number = (/0,0,1/)
!!$        else
!!$          !bloss/qt: only total water variable is needed for no-precip, 
!!$          !            fixed droplet number, warm cloud and no cloud simulations.
!!$           flag_wmass  = (/1/)
!!$           flag_precip = (/0/)
!!$           flag_number = (/0/)
!!$        end if
!!$     end if
!!$  end if

  ! output all microphysical fields to 3D output files if using more than
  !   just docloud.  Otherwise, rely on basic SAM outputs
  if(docloud.AND.(doprecip.OR.dopredictNc)) then
     flag_micro3Dout = 1
  end if

  ! initialize factor for latent heat
  lfac(:) = 1. ! use one as default for number species
  lfac(iqv) = lcond
  if((.NOT.dototalwater).AND.docloud) lfac(iqcl) = lcond
  if(doprecip) lfac(iqr) = lcond
  if(doicemicro) then
     lfac(iqci) = lsub
     lfac(iqs) = lsub
     if(dograupel) lfac(iqg) = lsub
  end if

  call graupel_init() ! call initialization routine within mphys module

  shiptrack_timeflag = .true.

  if((doprogaerosol.AND.doshiptrack2D).AND.(time.LT.shiptrack_time0)) then
    shiptrack_timeflag = .false.
  end if

  if(nrestart.eq.0) then

 ! compute initial profiles of liquid water - M.K.
      call satadj_liquid(nzm,tabs0,q0,qc0,pres*100.)

     ! initialize microphysical quantities
     if(dototalwater) q0 = q0 + qc0
     do k = 1,nzm
       micro_field(:,:,k,iqv) = q0(k)
        if(.NOT.dototalwater) micro_field(:,:,k,iqcl) = qc0(k)
        cloudliq(:,:,k) = qc0(k)
        tabs(:,:,k) = tabs0(k)
        !bloss: approx initialization of effective radius based on 
        !  Hugh's formula.  Here, I'm taking the ratio of the gamma functions
        !  to be about two when they're all pulled inside the cube root.  
        !  Not perfect, but should be a reasonable approximation, I think.
        if (qc0(k).gt.0.) then
          if(dofix_pgam) then
            tmp_pgam = pgam_fixed
          else
            tmp_pgam=0.0005714*(Nc0*RHO(K))+0.2714
            tmp_pgam = MAX(2.,MIN(10.,1./(tmp_pgam**2)-1.))
          end if
          
          tmp_lambda = ( (3.14159*1000./6.)*1.e6*Nc0 &
               *(tmp_pgam+3.)*(tmp_pgam+2.)*(tmp_pgam+1.) / qc0(k) )**(1./3.)
          tmp_lambda = MAX((tmp_pgam+1.)/60.e-6,MIN((tmp_pgam+1.)/1.e-6, &
               tmp_lambda))

          CloudLiquidGammaExponent(:,:,k) = tmp_pgam
          CloudLiquidLambda(:,:,k) = tmp_lambda
          CloudLiquidMassMixingRatio(:,:,k) = qc0(k)

          reffc(:,:,k) = 1.e6 *(tmp_pgam+3.)/tmp_lambda/2.
          if(masterproc) write(*,*) 'Experimental reffc initialization: ', reffc(1,1,k)
        else
          reffc(:,:,k) = 25.
        end if
        reffi(:,:,k) = 25.
        reffs(:,:,k) = 25.


      end do
     if(dopredictNc) then ! initialize concentration somehow...
       do k = 1,nzm
         !bloss(2019-07-25): Only initialize non-zero N_cloud in cloudy grid cells.
         !   Used to be q0.gt.0. which is true everywhere
         if(qc0(k).gt.0.) then
            micro_field(:,:,k,incl) = 1.e6*Nc0 ! choose to make the number mixing ratio equal to Nc0
         end if
       end do
     end if

     !===== Aerosol initialization =========
     !bloss(2018-02): Patching in Berner aerosol initialization code here
     
     if(doprogaerosol) then

       !brnr/bloss: Mean aerosol particle mass
       InitialMode1MeanMass = rho_aerosol*4.*pi/3.*(aer_rm1**3)*EXP(9.*(log(aer_sig1)**2)/2.)
       InitialMode2MeanMass = rho_aerosol*4.*pi/3.*(aer_rm2**3)*EXP(9.*(log(aer_sig1)**2)/2.)

       !bloss: Restructure aerosol initialization
       !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
       do k = 1,nzm
         do j = 1,ny
           do i = 1,nx
             ! initialize cloud aerosol as having mean mass of mode 1
             micro_field(i,j,k,inacc) = aer_n1 !bloss: nacc = Nc + NAd
             micro_field(i,j,k,iqacc) = InitialMode1MeanMass*massfactor*micro_field(i,j,k,inacc) !bloss: qacc = QAw + QAd
           end do
         end do
       end do

       if(.NOT.doprecoff) micro_field(:,:,:,iqar) = 0.

       select case(aerinitmode)
       case(0)

         ! Use sounding information from aerosol input
         !   This currently only works with SCAM IOP netcdf forcings
         !   with accumulation mode aerosol number and mass mixing ratios
         !   provided in the variables Na_accum and qa_accum.
         if((.NOT.allocated(AccumAerosolNumber_snd)).AND.(.NOT.allocated(AccumAerosolMass_snd))) then
           write(*,*) 'Error in MICRO_M2005_PA: initial conditions for accumulation mode'
           write(*,*) '  aerosol when using aerinitmode==0 should be specified in '
           write(*,*) '  IOP netcdf file using the variables Na_accum and qa_accum.'
           write(*,*) '********* Model stopping ... *************'
           call task_abort()
         end if

         if((MINVAL(AccumAerosolNumber_snd).LT.0.).OR.MINVAL(AccumAerosolMass_snd).LT.0.) then
           write(*,*) 'Error in MICRO_M2005_PA: Input initial condition for accumulation'
           write(*,*) '  mode aerosol number and/or mass are not positive definite'
           write(*,*) '  MAX/MIN(AccumAerosolNumber_snd) = ', MAXVAL(AccumAerosolNumber_snd), &
                MINVAL(AccumAerosolNumber_snd)
           write(*,*) '  MAX/MIN(AccumAerosolMass_snd) = ', MAXVAL(AccumAerosolMass_snd), &
                MINVAL(AccumAerosolMass_snd)
           write(*,*) 'If these numbers are -9999, they were likely missing values in the IOP netcdf file'
           write(*,*) '********* Model stopping ... *************'
           call task_abort()
         end if

         !borrow code from forcing.f90
         call InterpolateFromForcings(nsnd,nzsnd,daysnd,zsnd,psnd,AccumAerosolNumber_snd, &
              nzm,day,z,pres,Na_accum,.true.)
         call InterpolateFromForcings(nsnd,nzsnd,daysnd,zsnd,psnd,AccumAerosolMass_snd, &
              nzm,day,z,pres,qa_accum,.true.)

         if(masterproc) then
         write(*,*) 'In MICRO_M2005_PA, using aerinitmode==0.  This means that '
         write(*,*) '  the accumulation mode aerosol is initialized using soundings '
         write(*,*) '  from the IOP netcdf forcing file.  Initial profile below ...'
         write(*,*) 
         do k = 1,nzm
           write(*,845) z(k), Na_accum(k), qa_accum(k), &
                ( qa_accum(k)/Na_accum(k) / EXP(9.*(log(aer_sig1)**2)/2.) &
                / ( rho_aerosol*4.*pi/3. ) )**(1./3.) 
845        format('Accum Mode Aerosol: z, N, q, D = ', F10.2,3E14.6)
         end do
         end if

         !bloss: Restructure aerosol initialization
         !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
         do k = 1,nzm
           do j = 1,ny
             do i = 1,nx
               ! Nacc = NAd + NC (accumulation mode = dry aerosol + wet (in-cloud-droplet) aerosol).
               micro_field(i,j,k,inacc) = Na_accum(k)
               micro_field(i,j,k,iqacc) = qa_accum(k)
             end do
           end do
         end do

       case(1)

         ! Use aerosol mode 2 (aer_rm2, aer_n2, aer_sig2) to initialize aerosol above the inversion

         ! find inversion
         !bloss: identify inversion as the height where
         !   -d(RH)/dz*d(sl)/dz 
         ! is maximum.  The usual criteria of the mean height of maximum
         ! theta_l gradient might not be so robust with mean profiles, so
         ! here we include the RH lapse rate as well.
         kinv = 0
         do k = 1,nzm
           rh0(k) = qv0(k)/qsatw(tabs0(k),pres(k))
         end do

         tmp_max = -1.
         do k = 1,nzm-1
           tmp_check = - (rh0(k+1)-rh0(k)) * (t0(k+1)-t0(k)) &
                / (z(k+1)-z(k))**2
           if(tmp_check.gt.tmp_max) then
             kinv = k
             tmp_max = tmp_check !mwyant: need to update tmp_max!
           end if
         end do
         if(masterproc) write(*,826) z(kinv)
826      format('Initial inversion height is ',F8.1, &
              ' meters.  This is used in the aerosol initialization.')

         ! set FT to aerosol size spec [aer_rm2,aer_n2,aer_sig2] aer_sig2 should equal aer_sig1
         do k = kinv+1,nzm
           !bloss(2018-02): In cloud-free FT, set dry aerosol concentration to aer_n1
           ! If any cloud layers exist, partition aer_n2 and aerosol mass between cloud and dry aerosol
           ! dry_ratio is the fraction of total aerosols that are not cloud droplets
           do j = 1,ny
             do i = 1,nx
               !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
               micro_field(i,j,k,inacc) = aer_n2
               micro_field(i,j,k,iqacc) = InitialMode2MeanMass*massfactor*micro_field(i,j,k,inacc)
             end do
           end do
         end do
         
         if(ldebug) then
           write(*,*)
           write(*,*) 'aer_n1, aer_n2 = ', aer_n1, aer_n2
           write(*,*) 'aer_rm1, aer_rm2 = ', aer_rm1, aer_rm2
           write(*,*) 'aer_sig1, aer_sig2 = ', aer_sig1, aer_sig2

           write(*,*)
           do k = 1,kinv
             write(*,827) k, z(k), rh0(k), micro_field(1,1,k,inacc), micro_field(1,1,k,iqacc), micro_field(1,1,k,incl)
827          format('k,z,relh,nacc,qacc,ncl = ', I4,F8.2,F8.4,3E10.2)
           end do
           write(*,*) 'Inversion'
           do k = kinv+1,nzm
             write(*,827) k, z(k), rh0(k), micro_field(1,1,k,inacc), micro_field(1,1,k,iqacc), micro_field(1,1,k,incl)
           end do
!Bloss           STOP 'in micro_init'

         end if

       end select
     end if !doprogaerosol       

     ! initialize microphysical quantities

     if(doprogaerosol.AND.(doPOCinit_Na.OR.doPOCinit_qt)) then

       call task_rank_to_index(rank,it,jt)

       nx_Nc1 = int(floor(OVC_length/dx))
       nx_trans = int(ceiling(POC_trans_length/dx))
       nx_Nc2 = int(nx_gl - 2*(nx_trans + nx_Nc1))

       POCperturb_Na = POCperturb_Na*1.E6
       POCperturb_Qt = POCperturb_Qt*1.e-3

       POC_length = real(nx_Nc2*dx)

       tmpx1 = (real(nx_Nc1) - 0.5)*dx
       tmpx2 = (real(nx_Nc1 + nx_trans + nx_Nc2) - 0.5)*dx

!bloss(2018-02): RHOA --> rho_aerosol
!bloss(2019-07): pi/6 --> 4.*pi/3. because mass ~ 4/3*pi*r^3
!!$       dumM2 = (POCperturb_Na)*RHOA*pi/6.*(aer_rm1**3)*EXP(9.*(log(aer_sig1)**2)/2.) 
       dumM2 = (POCperturb_Na)*rho_aerosol*4.*pi/3.*(aer_rm1**3)*EXP(9.*(log(aer_sig1)**2)/2.) 

       do k = 1,kinv ! only initialize up to inversion height
         ! Note that qv and qcl are initialized above
         do m = 1,nx

           tmp_ind = modulo(it + m,nx_gl)

           if( ( tmp_ind > nx_Nc1) .AND. (tmp_ind <= (nx_Nc1 + nx_trans)) ) then
             arg = 1/POC_trans_length*((tmpx1-OVC_length) + dx*(tmp_ind - nx_Nc1))
             !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
             micro_field(m,:,k,inacc) = micro_field(m,:,k,inacc) + (-1.*POCperturb_Na)/2.*(cos(arg*pi)+1) + POCperturb_Na
             micro_field(m,:,k,iqacc) = micro_field(m,:,k,iqacc) + (-1.*dumM2)/2.*(cos(arg*pi)+1) + dumM2
             micro_field(m,:,k,iqv) = micro_field(m,:,k,iqv) + (-1.*POCperturb_qt)/2.*(cos(arg*pi)+1) + POCperturb_qt
           else if( ( tmp_ind > (nx_Nc1 + nx_trans)) .AND. ( tmp_ind <= (nx_Nc1 + nx_trans + nx_Nc2)) ) then
             !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
             micro_field(m,:,k,inacc) = micro_field(m,:,k,inacc) + POCperturb_Na
             micro_field(m,:,k,iqacc) = micro_field(m,:,k,iqacc) + dumM2
             micro_field(m,:,k,iqv) = micro_field(m,:,k,iqv) + POCperturb_qt
           else if( ( tmp_ind > (nx_Nc1 + nx_trans + nx_Nc2)) .AND. ( tmp_ind <= (nx_Nc1 + 2*nx_trans + nx_Nc2)) ) then
             arg = 1/POC_trans_length*((tmpx2 - OVC_length - POC_trans_length - POC_length) + dx*(tmp_ind - (nx_Nc1+nx_trans+nx_Nc2)))
             !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
             micro_field(m,:,k,inacc) = micro_field(m,:,k,inacc) + (-1.*POCperturb_Na)/2.*(cos(arg*pi+pi)+1) + POCperturb_Na
             micro_field(m,:,k,iqacc) = micro_field(m,:,k,iqacc) + (-1.*dumM2)/2.*(cos(arg*pi+pi)+1) + dumM2
             micro_field(m,:,k,iqv) = micro_field(m,:,k,iqv) + (-1.*POCperturb_qt)/2.*(cos(arg*pi+pi)+1) + POCperturb_qt
           end if

         end do

         if (qc0(k).GT.0.) then
           do j = 1,ny
             do i = 1,nx
               micro_field(i,j,k,incl) = min(micro_field(i,j,k,inacc), 1.e6*Nc0) !removed factor of rho from Nc0 

             end do
           end do
         end if
       end do !k = 1,nzm

       !PLACEHOLDER FOR FUTURE WORK:
       !   Perturb microphysical variables in Na and/or qt to trigger POC formation

     end if

     if(doprogaerosol) then
       do k = 1,nzm
         do j = 1,ny
           do i = 1,nx
             if(micro_field(i,j,k,incl).gt.EPSILON(1.)*micro_field(i,j,k,inacc)) then
               !bloss(2020-11): Where cloud exists, make a diagnostic paritioning of the
               !   aerosol mass here in case the initial condition is output
               micro_field(i,j,k,inad) = MAX(0., micro_field(i,j,k,inacc) - micro_field(i,j,k,incl) )

               micro_field(i,j,k,iqad) = micro_field(i,j,k,iqacc) &
                    *DryAerosolMassFraction( micro_field(i,j,k,inad)/micro_field(i,j,k,inacc), &
                    aer_sig1 , 2 )
               micro_field(i,j,k,iqaw) = MAX(0., micro_field(i,j,k,iqacc) - micro_field(i,j,k,iqad) )
             else
               ! all aerosol is dry
               micro_field(i,j,k,inad) = micro_field(i,j,k,inacc)
               micro_field(i,j,k,iqad) = micro_field(i,j,k,iqacc)
               micro_field(i,j,k,iqaw) = 0.
             end if
           end do
         end do
       end do
     end if

   end if

   if(docloud) call micro_diagnose()   ! leave this here

   !------------------------------------------------------------------
   !  Initialize lookup table for M2011 scavenging if scheme active
   !  This is done whether the run is fresh or restarted
   if(do_m2011_scavenge.AND.doprogaerosol) then
     call t_startf('scav_init')
     call memory('allocate')
     call init_scavenging
     call t_stopf('scav_init')

     if(.NOT.isallocatedSCAV3D) then
       !bloss(2020-10): set up 3D array of scavenging tendencies 
       !  This allows the tendencies to be stored and re-used since
       !  the scavenging routine is not called every icycle.
       ALLOCATE(scav3d(nx,ny,nzm,4),scavname(4),scavlongname(4),scavunits(4), &
            STAT=ierr)
       if(ierr.ne.0) then
         write(*,*) 'Error in allocating scavenging outputs'
         call task_abort()
       end if
       scav3d(:,:,:,:) = 0.
       isallocatedSCAV3D = .true.
     end if

     if(do_scav_3d_output) then
       !bloss(2020-10): Add 3d outputs of scavenging tendencies
       nfields3D_micro = nfields3D_micro + 4 ! write out 3D scavenging tendencies

       ! set up output names and units
       n = 1
       scavname(n) = 'SCVTQAR'
       scavlongname(n) = 'Aerosol mass tendency due to interestitial scavenging by rain (Source of QAr, Sink of QAd)'
       scavunits(n) = 'kg/kg/s'
       n = n + 1
       scavname(n) = 'SCVTQAW'
       scavlongname(n) = 'Aerosol mass tendency due to interestitial scavenging by cloud (Source of QAw, Sink of QAd)'
       scavunits(n) = 'kg/kg/s'
       n = n + 1
       scavname(n) = 'SCVTNADR'
       scavlongname(n) = 'Aerosol number tendency due to interestitial scavenging by rain (Sink of NAd)'
       scavunits(n) = '#/kg/s'
       n = n + 1
       scavname(n) = 'SCVTNADC'
       scavlongname(n) = 'Aerosol number tendency due to interestitial scavenging by cloud (Sink of NAd)'
       scavunits(n) = '#/kg/s'
     end if
   end if

  ! set up names, units and scales for these microphysical quantities
  if(dototalwater) then
    mkname(iqv) = 'QTO'
    mklongname(iqv) = 'TOTAL WATER (VAPOR + CLOUD LIQUID)'
    mkunits(iqv) = 'g/kg'
    mkoutputscale(iqv) = 1.e3

  else
    mkname(iqv) = 'QV'
    mklongname(iqv) = 'WATER VAPOR'
    mkunits(iqv) = 'g/kg'
    mkoutputscale(iqv) = 1.e3

    if(docloud) then
      mkname(iqcl) = 'QC'
      mklongname(iqcl) = 'CLOUD WATER'
      mkunits(iqcl) = 'g/kg'
      mkoutputscale(iqcl) = 1.e3
    end if
  end if

  if(docloud.AND.dopredictNc) then
    mkname(incl) = 'NC'
    mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
    mkunits(incl) = '#/mg'
    mkoutputscale(incl) = 1.e-6
  end if

 if(doprogaerosol) then
    mkname(iqacc) = 'QAcc'
    mklongname(iqacc) = 'ACCUMULATION MODE AEROSOL MASS (=QAd+QAw, dry+wet)'
    mkunits(iqacc) = 'g/kg'
    mkoutputscale(iqacc) = 1.e3/massfactor
    
    mkname(inacc) = 'NAcc'
    mklongname(inacc) = 'ACCUMULATION MODE AEROSOL NUMBER (=NAd+NC, dry+NC)'
    mkunits(inacc) = '#/mg'
    mkoutputscale(inacc) = 1.e-6
    
    mkname(iqad) = 'QAd'
    mklongname(iqad) = 'DRY AEROSOL MASS'
    mkunits(iqad) = 'g/kg'
    mkoutputscale(iqad) = 1.e3/massfactor
    
    mkname(inad) = 'NAd'
    mklongname(inad) = 'DRY AEROSOL NUMBER CONCENTRATION'
    mkunits(inad) = '#/mg'
    mkoutputscale(inad) = 1.e-6
    
    mkname(iqaw) = 'QAw'
    mklongname(iqaw) = 'WET AEROSOL MASS'
    mkunits(iqaw) = 'g/kg'
    mkoutputscale(iqaw) = 1.e3/massfactor

    if(.NOT.doprecoff) then
       mkname(iqar) = 'QAr'
       mklongname(iqar) = 'RAIN AEROSOL MASS'
       mkunits(iqar) = 'g/kg'
       mkoutputscale(iqar) = 1.e3/massfactor
    end if
 end if
  
  if(doprecip) then
     mkname(iqr) = 'QR'
     mklongname(iqr) = 'RAIN'
     mkunits(iqr) = 'g/kg'
     mkoutputscale(iqr) = 1.e3

     mkname(inr) = 'NR'
     mklongname(inr) = 'RAIN NUMBER CONCENTRATION'
     mkunits(inr) = '#/mg'
     mkoutputscale(inr) = 1.e-6
  end if

  if(doicemicro) then
     mkname(iqci) = 'QI'
     mklongname(iqci) = 'CLOUD ICE'
     mkunits(iqci) = 'g/kg'
     mkoutputscale(iqci) = 1.e3

     mkname(inci) = 'NI'
     mklongname(inci) = 'CLOUD ICE NUMBER CONCENTRATION'
     mkunits(inci) = '#/mg'
     mkoutputscale(inci) = 1.e-6

     mkname(iqs) = 'QS'
     mklongname(iqs) = 'SNOW'
     mkunits(iqs) = 'g/kg'
     mkoutputscale(iqs) = 1.e3

     mkname(ins) = 'NS'
     mklongname(ins) = 'SNOW NUMBER CONCENTRATION'
     mkunits(ins) = '#/mg'
     mkoutputscale(ins) = 1.e-6

     if(dograupel) then
       if(dohail) then
         ! dense ice is hail
        mkname(iqg) = 'QH'
        mklongname(iqg) = 'HAIL'
        mkunits(iqg) = 'g/kg'
        mkoutputscale(iqg) = 1.e3

        mkname(ing) = 'NH'
        mklongname(ing) = 'HAIL NUMBER CONCENTRATION'
        mkunits(ing) = '#/mg'
        mkoutputscale(ing) = 1.e-6
      else
         ! dense ice is graupel
        mkname(iqg) = 'QG'
        mklongname(iqg) = 'GRAUPEL'
        mkunits(iqg) = 'g/kg'
        mkoutputscale(iqg) = 1.e3

        mkname(ing) = 'NG'
        mklongname(ing) = 'GRAUPEL NUMBER CONCENTRATION'
        mkunits(ing) = '#/mg'
        mkoutputscale(ing) = 1.e-6
      end if
     end if
end if

   if(mod(nstatfrq,2).eq.0) then
     nskip_quickbeam = nstatfrq/2 ! default is to call twice per statistics output.
   else
     nskip_quickbeam = nstatfrq ! if nstatfrq is odd, call once per stat output.
   end if

  if(do_chunk_mkbudget) then
    ! water vapor budget
    mkbudget_name(1) = 'NAT'; mkbudget_longname(1) = 'total accumulation-mode aerosol number (NAd+NC+NR)'; mkbudget_units(1) = '#/kg'
    mkbudget_name(2) = 'QAT'; mkbudget_longname(2) = 'total accumulation-mode aerosol mass (QAd+QAw+QAr)'; mkbudget_units(2) = 'kg/kg'
    mkbudget_name(3) = TRIM(mkname(iqv)); mkbudget_longname(3) = TRIM(mklongname(iqv)); mkbudget_units(3) = 'kg/kg'
    mkbudget_name(4) = TRIM(mkname(iqr)); mkbudget_longname(4) = TRIM(mklongname(iqr)); mkbudget_units(4) = 'kg/kg'
    if(.NOT.dototalwater) then
      mkbudget_name(5) = TRIM(mkname(iqcl)); mkbudget_longname(5) = TRIM(mklongname(iqcl)); mkbudget_units(5) = 'kg/kg'
    end if

    ! extra budget terms for rain evaporation
    n = 1
    mkbudget_extra_name(n) = 'PRC'
    mkbudget_extra_longname(n) = 'Tendency of QC (cloud liquid mass mixing ratio) due to autoconversion'
    mkbudget_extra_units(n) = 'kg/kg/s'
    n = n + 1
    mkbudget_extra_name(n) = 'NPRC'
    mkbudget_extra_longname(n) = 'Tendency of NC (cloud liquid number mixing ratio) due to autoconversion'
    mkbudget_extra_units(n) = '#/kg/s'
    n = n + 1
    mkbudget_extra_name(n) = 'NPRC1'
    mkbudget_extra_longname(n) = 'Tendency of NR (rain number mixing ratio) due to autoconversion'
    mkbudget_extra_units(n) = '#/kg/s'
    n = n + 1
    mkbudget_extra_name(n) = 'PRA'
    mkbudget_extra_longname(n) = 'Tendency of QR (rain mass mixing ratio) due to accretion of cloud droplets by rain'
    mkbudget_extra_units(n) = 'kg/kg/s'
    n = n + 1
    mkbudget_extra_name(n) = 'NPRA'
    mkbudget_extra_longname(n) = 'Tendency of NC (cloud liquid number mixing ratio) due to accretion of cloud droplets by rain'
    mkbudget_extra_units(n) = '#/kg/s'
    n = n + 1
    mkbudget_extra_name(n) = 'PRE'
    mkbudget_extra_longname(n) = 'Tendency of QR (rain mass mixing ratio) due to rain evaporation'
    mkbudget_extra_units(n) = 'kg/kg/s'
    n = n + 1
    mkbudget_extra_name(n) = 'NSUBR'
    mkbudget_extra_longname(n) = 'Tendency of NR (rain number mixing ratio) due to rain evaporation'
    mkbudget_extra_units(n) = '#/kg/s'

    if(do_m2011_scavenge) then
      n = n + 1
      mkbudget_extra_name(n) = 'SCVTQAR'
      mkbudget_extra_longname(n) = 'Aerosol mass tendency due to interestitial scavenging by rain (Source of QAr, Sink of QAd)'
      mkbudget_extra_units(n) = 'kg/kg/s'
      n = n + 1
      mkbudget_extra_name(n) = 'SCVTQAW'
      mkbudget_extra_longname(n) = 'Aerosol mass tendency due to interestitial scavenging by cloud (Source of QAw, Sink of QAd)'
      mkbudget_extra_units(n) = 'kg/kg/s'
      n = n + 1
      mkbudget_extra_name(n) = 'SCVTNADR'
      mkbudget_extra_longname(n) = 'Aerosol number tendency due to interestitial scavenging by rain (Sink of NAd)'
      mkbudget_extra_units(n) = '#/kg/s'
      n = n + 1
      mkbudget_extra_name(n) = 'SCVTNADC'
      mkbudget_extra_longname(n) = 'Aerosol number tendency due to interestitial scavenging by cloud (Sink of NAd)'
      mkbudget_extra_units(n) = '#/kg/s'
    end if

    if(n.ne.n_mkbudget_extra) then
      write(*,*) 'Error in setting up mkbudget_extra output fields in MICRO_M2005_PA/microphysics.f90'
      write(*,*) '=== Mismatch between n_mkbudget_extra and number of variables set up for output'
    end if

  end if

end subroutine micro_init

!----------------------------------------------------------------------
subroutine micro_finalize()
  implicit none
  integer :: ierr

   if(do_m2011_scavenge.AND.doprogaerosol) then
      call memory('deallocate')
   end if

  if(isallocatedMICRO) then
     ! allocate microphysical variables
    deallocate(micro_field, fluxbmk, fluxtmk, &
          reffc, reffr, reffi, reffs, &
          CloudLiquidMassMixingRatio, CloudLiquidGammaExponent, &
          CloudLiquidLambda, CloudIceMassMixingRatio, SnowMassMixingRatio, &
          mkwle, mkwsb, mkadv, mkdiff, mklsadv, mkstor, &
          stend, mtend, &
          mtend3d, mfrac, trtau, micro_proc_rates, mksed, tmtend, &
          cloudliq, tmtend3d, flag_micro3Dout, flag_wmass, flag_precip, &
          flag_number, lfac, mkname, mklongname, mkunits, mkoutputscale, &
          mk0, is_water_vapor, scvtndqadclstat,scvtndqadrstat, &
          scvtndnadclstat,scvtndnadrstat, STAT=ierr)
     if(ierr.ne.0) then
        write(*,*) 'Failed to deallocate microphysical arrays on proc ', rank
     end if
   end if
   
 end subroutine micro_finalize

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

integer :: i,j
real, dimension(nx, ny) :: Nflux_acc, Nflux_ait 

! shipv2 stuff
real, dimension(nx, ny) :: ship_flux
real :: t_ship, x0_ship, y0_ship, x1_ship, y1_ship, x2_ship, y2_ship, &
     tmpr, tmp1(1), tmp2(1), tmpx(nx), tmpy(ny)
integer :: it, jt, n


fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero

fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
fluxtmk(:,:,index_water_vapor) = fluxtq(:,:) ! top of domain qv flux

if(doprogaerosol) then
  select case(aerfluxmode)
    case(1) !fixed flux
      fluxbmk(:,:,iqacc) = fluxQAd0*massfactor*(1/RHOW(1))
      fluxbmk(:,:,inacc) = fluxNAd0*(1/RHOW(1))
    
    case(2)  !flux based on wind speed
       !bloss(2019-07-20): updated treatment from Wyant that allows
       ! aerosol size to be specified along with number flux
       !  parameter 
       Nflux_acc(:,:) = sflux_nacc_coef*whitecap_coef*(u10arr(:,:))**3.41  ! units are #/m2/s (Berner 2013 eq(5))
       fluxbmk(:,:,inacc) = Nflux_acc(:,:)/RHOW(1)  ! units are #/kg m/s
       fluxbmk(:,:,iqacc) = 1.3333 * pi * (sflux_racc*1.e-6)**3 * & 
          rho_aerosol * Nflux_acc(:,:) * &
          exp(4.5 * log(aer_sig1)**2)/RHOW(1) ! units are kg/m2/s / [rho air]

!bloss      do i=1,nx
!bloss         do j=1,ny
!bloss            fluxbmk(i,j,iqad) = fluxQAd0*massfactor*(1/RHOW(1))*(3.8e-6)*(u10arr(i,j))**3.41 
!bloss            fluxbmk(i,j,inad) = fluxNAd0*(1/RHOW(1))*(3.8e-6)*(u10arr(i,j))**3.41
!bloss         end do
!bloss      end do

    case default
      fluxbmk(:,:,inacc) = 0.
      fluxbmk(:,:,iqacc) = 0.
    end select

    if(ldebug.AND.masterproc) then
      write(*,*) 'Shipv2 = ', doshipv2
      write(*,*) 'shipv2_time0 = ', shipv2_time0, ' day = ', day, ' shipv2_timef = ', shipv2_timef
    end if

    ! add ship plume when doshipv2==.true.
    if(doshipv2) then
      if ((day.ge.shipv2_time0).AND.(day.le.shipv2_timef)) then
        if(masterproc.AND.ldebug) then
          write(*,*) 'Adding ship emissions, day = ', day
        end if

        select case (shipv2_mode)
        case (0)

          t_ship = 86400*(day - shipv2_time0)
          x0_ship = 0.5*dx*nx_gl - ug*t_ship
          y0_ship = - vg*t_ship

          !account for periodicity, allow for ten passages through domain at most
          if(x0_ship.gt.dx*nx_gl) then
            do n = 1,10
              if(x0_ship.gt.dx*nx_gl) x0_ship = x0_ship - dx*nx_gl
            end do
          end if

          if(x0_ship.lt.0.) then
            do n = 1,10
              if(x0_ship.lt.0.) x0_ship = x0_ship + dx*nx_gl
            end do
          end if

          if(y0_ship.gt.dy*ny_gl) then
            do n = 1,10
              if(y0_ship.gt.dy*ny_gl) y0_ship = y0_ship - dy*ny_gl
            end do
          end if

          if(y0_ship.lt.0.) then
            do n = 1,10
              if(y0_ship.lt.0.) y0_ship = y0_ship + dy*ny_gl
            end do
          end if

          if(masterproc) write(*,*) 'x0_ship, y0_ship = ', x0_ship, y0_ship

!!$        x1_ship = x0_ship + SIGN(dx*float(nx_gl), ug)
!!$        y1_ship = y0_ship + SIGN(dy*float(ny_gl), vg)
!!$
!!$        x2_ship = x0_ship + 2*SIGN(dx*float(nx_gl), ug)
!!$        y2_ship = y0_ship + 2*SIGN(dy*float(ny_gl), vg)

          ! figure out x and y coordinates on current processor
          call task_rank_to_index(rank,it,jt)
          do i = 1,nx
            tmpx(i) = dx*( float(i+it) - 0.5)
          end do
          do j = 1,ny
            tmpy(j) = dy*( float(j+jt) - 0.5)
          end do

          ship_flux(:,:) = 0.
          do j = 1,ny
            do i = 1,nx
!!$            tmpr = SQRT( &
!!$                   MIN(ABS(tmpx(i)-x0_ship),ABS(tmpx(i)-x1_ship),ABS(tmpx(i)-x2_ship))**2 &
!!$                 + MIN(ABS(tmpy(j)-y0_ship),ABS(tmpy(j)-y1_ship),ABS(tmpy(j)-y2_ship))**2 )
              tmpr = SQRT( (tmpx(i)-x0_ship)**2 + (tmpy(j)-y0_ship)**2 )
              if(tmpr.lt.shipv2_radius) then
                ship_flux(i,j) = 0.5*(1. + cos(pi * MIN(1., tmpr / shipv2_radius ) ) )
              end if
            end do
          end do
          tmp1(1) = SUM(ship_flux(1:nx,1:ny))
          if(dompi) then
            call task_sum_real(tmp1,tmp2,1)
            tmp1(1) = tmp2(1)
          end if
          if(ldebug) then
            write(*,*) 'dx/dy/sum = ', rank, dx, dy, tmp1(1), maxval(ship_flux(:,:))
            write(*,*) 'x*_ship= ', rank, x0_ship !, x1_ship, x2_ship
            write(*,*) 'y*_ship = ', rank, y0_ship !, y1_ship, y2_ship
            write(*,*) 'min/max x = ', rank, minval(tmpx(:)), maxval(tmpx(:))
            write(*,*) 'min/max y = ', rank, minval(tmpy(:)), maxval(tmpy(:))
          end if
          if(tmp1(1).gt.EPSILON(1.)) then
            ship_flux(:,:) = ship_flux(:,:) / tmp1(1) / dx / dy
          else
            ! ship has left the domain.
            ship_flux(:,:) = 0.
          end if

          fluxbmk(1:nx,1:ny,inacc) = fluxbmk(1:nx,1:ny,inacc) + ship_flux(1:nx,1:ny) * shipv2_aer_n / rhow(1)

          fluxbmk(1:nx,1:ny,iqacc) = fluxbmk(1:nx,1:ny,iqacc) + ship_flux(1:nx,1:ny) &
               * 1.3333 * pi * (shipv2_aer_r)**3 &
               * rho_aerosol * shipv2_aer_n  &
               * exp(4.5 * log(aer_sig1)**2)/RHOW(1)


        case default
          write(*,*) 'No flux at surface due to ship track -- wrong value of shipv2_mode'
        end select

      end if ! if within shipv2 plume release times
    end if ! if(doshipv2)

    ! accumulate aerosol surface fluxes so that we can capture fluxes from transient ship tracks
    aer_nflux_avg_xy(:,:) = aer_nflux_avg_xy(:,:) + fluxbmk(:,:,inacc)*dtfactor
    aer_qflux_avg_xy(:,:) = aer_qflux_avg_xy(:,:) + fluxbmk(:,:,iqacc)*dtfactor

  end if ! if doprogaerosol

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., that is  all the microphysics processes except for the spatial transport happen.

! IMPORTANT: You need to use the thermodynamic constants like specific heat, or
! specific heat of condensation, gas constant, etc, the same as in file params.f90
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg) 
! It should not be changed during all of your point microphysical processes!

subroutine micro_proc()

implicit none

real, dimension(nzm) :: &
     tmpqcl, tmpqci, tmpqr, tmpqacc, tmpqad, tmpqaw, tmpqar, tmpqs, tmpqg, tmpqv, &
     tmpncl, tmpnci, tmpnr, tmpnacc, tmpnad, tmpns, tmpng, tmpzero, tmpcmd, tmpcmddry, tmprhw, &
     tmpnuc, tmpnur, tmpw, tmpwsub, tmppres, tmpdz, tmptabs, &
     tmtend1d, &
     mtendqcl, mtendqci, mtendqr, mtendqad, mtendqaw, mtendqar, mtendqs, mtendqg, mtendqv, &
     mtendncl, mtendnci, mtendnr, mtendnad, mtendns, mtendng,  &
     stendqcl, stendqci, stendqr, stendqad, stendqaw, stendqar, stendqs, stendqg, stendqv, &
     stendncl, stendnci, stendnr, stendnad, stendns, stendng,  &
     scvtndqadcl,scvtndqadr,scvtndnadcl,scvtndnadr, &
     tmpec3d, &
     effg1d, effr1d, effs1d, effc1d, effi1d, &
     tmp_cl_pgam, tmp_cl_lambda
    
real, dimension(nzm,nmicro_fields) :: stend1d, mtend1d
real :: tmpc, tmpr, tmpi, tmps, tmpg
integer :: i1, i2, j1, j2, i, j, k, m, n

real(8) :: tmp_total, tmptot

logical :: do_accumulate_process_rates

!bloss: extra outputs for chunk-averaged microphysical process rates
real :: proc_extra(nzm,n_mkbudget_extra)

!brnr: variables used for fixed nc
real :: tmpRH, tmpq, tmp, tgrad_max(nx,ny),tgrad
integer :: k_inv_ind(nx,ny), it, jt, kinv

integer :: kc, kb

real :: rdry(1), rwet(1), hygro_arr(1), relhum_arr(1)
integer :: npoints = 1
real :: ddry

!bloss: cloudradar arrays
integer :: nprof, ngate
real :: esat1, qsat1
real*8, dimension(nx,nzm) :: hgt_matrix, p_matrix, t_matrix, rh_matrix, & ! inputs
     Ze_non, Ze_ray, h_atten_to_vol, g_atten_to_vol, dBZe, &
     g_to_vol_in, g_to_vol_out
real*8, dimension(nhclass,nx,nzm) :: hm_matrix, re_matrix, Np_matrix

logical, parameter :: do_new_scavenge = .false.

real :: scav_factor, scav_mass, scav_number
real :: dry_aerosol_mass_before, dry_aerosol_number_before

real, external :: qsatw

call t_startf ('micro_proc')

if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
      do i=1,nx
         precsfc(i,j)=0.
      end do
   end do
   do k=1,nzm
      precflux(k) = 0.
   end do
end if

if(do_chunked_energy_budgets) then
  if(mod(nstep-1,nsaveMSE).eq.0.and.icycle.eq.1) then
    ! initialize variables that will accumulate surface precipitation as a function of x,y
    prec_accum(:,:) = 0.
    prec_ice_accum(:,:) = 0.
    
    ! initialize variables that will accumulate 3D tendency due to sedimentation
    qtot_sed(:,:,:) = 0.
    qice_sed(:,:,:) = 0.

    if(do_chunk_mkbudget) then
      ! initialize mkbudget sedimentation terms
      mkbudget_sed(:,:,:,:) = 0.
    end if

    if(do_mkbudget_extra) then
      ! initialize mkbudget sedimentation terms
      mkbudget_extra(:,:,:,:) = 0.
    end if
  end if ! if(mod(nstep-1,nsaveMSE).eq.0.and.icycle.eq.1) 
end if ! if(do_chunked_energy_budgets)

if(dostatis) then ! initialize arrays for statistics
   mfrac(:,:) = 0.
   mtend(:,:) = 0.
   trtau(:,:) = 0.
   qpfall(:)=0.
   tlat(:) = 0.
   tmtend3d(:,:,:) = 0.

   micro_proc_rates(:,:) = 0.
   do_accumulate_process_rates = dostatis.AND.do_output_micro_process_rates

   if(do_m2011_scavenge.AND.doprogaerosol) then
     scvtndqadclstat(:) = 0.
     scvtndqadrstat(:) = 0.
     scvtndnadclstat(:) = 0.
     scvtndnadrstat(:) = 0.
   end if

end if
stend(:,:) = 0.
mksed(:,:) = 0.


dostatis_quickbeam = .false.
! only call quickbeam every nskip_quickbeam statistics step
if(doreflectivity_cloudradar) then
if(dostatis.AND.(mod(mod(nstep,nstat),nskip_quickbeam*nstatis).eq.0)) then
  dostatis_quickbeam = .true.
  factor_quickbeam = float(nskip_quickbeam)
  if(masterproc) write(*,*) 'Calling quickbeam this statistics step'
else
  factor_quickbeam = 0.
  if(masterproc.AND.dostatis) write(*,*) 'Skipping call of quickbeam this statistics step'
end if
end if
!!$if(doprecip) total_water_prec = total_water_prec + total_water()
 

mtend3d(:,:,:,:) = 0. !rezero micro3d fields 

do j = 1,ny
   do i = 1,nx

      ! zero out mixing ratios of microphysical species
      tmpqv(:) = 0.
      tmpqcl(:) = 0.
      tmpncl(:) = 0.
      tmpnuc(:) = 0. !gamma exponent
      tmpqr(:) = 0.
      tmpnr(:) = 0.
      tmpnur(:)= 0. !gamma exponent
      tmpqacc(:) = 0.
      tmpnacc(:) = 0.
      tmpqad(:) = 0.
      tmpnad(:) = 0.
      tmpcmd(:) = 0. !used for aerosol scavenging
      tmpcmddry(:) = 0.
      tmpqaw(:) = 0.
      tmpqar(:) = 0.
      tmpqci(:) = 0.
      tmpnci(:) = 0.
      tmpqs(:) = 0.
      tmpns(:) = 0.
      tmpqg(:) = 0.
      tmpng(:) = 0.

      tmp_cl_pgam(:) = 0.
      tmp_cl_lambda(:) = 0.

      stend1d(:,:) = 0.

      if(dototalwater) then
        ! get microphysical quantities in this grid column
        tmpqv(:) = micro_field(i,j,:,iqv) !bloss/qt: This is total water (qv+qcl)
        tmpqcl(:) = 0.
        !bloss/qt: compute cloud liquid below from saturation adjustment.
      else
        tmpqv(:) = micro_field(i,j,:,iqv) 
        tmpqcl(:) = micro_field(i,j,:,iqcl)
      end if

      if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)
      if(doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
         tmpnr(:) = micro_field(i,j,:,inr)
      end if

      if(doicemicro) then
         tmpqci(:) = micro_field(i,j,:,iqci)
         tmpnci(:) = micro_field(i,j,:,inci)
         tmpqs(:) = micro_field(i,j,:,iqs)
         tmpns(:) = micro_field(i,j,:,ins)
         if(dograupel) then
            tmpqg(:) = micro_field(i,j,:,iqg)
            tmpng(:) = micro_field(i,j,:,ing)
         end if
      end if

      !bloss(2020-11): Move below ice variables to match order in micro_field
      if(doprogaerosol) then
         !bloss(2020-11): Modify for new dry+cloud aerosol prognostic mass/number
         tmpqacc(:) = micro_field(i,j,:,iqacc)/massfactor
         tmpnacc(:) = micro_field(i,j,:,inacc)

         if(.NOT.doprecoff) tmpqar(:) = micro_field(i,j,:,iqar)/massfactor
         
         !bloss(2020-11): Partition total (dry+wet) aerosol mass into
         !   dry aerosol mass and wet (in-cloud-droplet) aerosol mass
         call PartitionAerosolMass( nzm, tmpnacc, tmpqacc, aer_sig1, &
              tmpncl, tmpnad, tmpqad, tmpqaw)

      end if

      ! get absolute temperature in this column
      !bloss/qt: before saturation adjustment for liquid,
      !          this is Tcl = T - (L/Cp)*qcl (the cloud liquid water temperature)
      tmptabs(:) = t(i,j,:)  &           ! liquid water-ice static energy over Cp
           - gamaz(:) &                                   ! potential energy
           + fac_cond * (tmpqcl(:) + tmpqr(:)) &    ! bloss/qt: liquid latent energy due to rain only if dototalwater==.true.
           + fac_sub  * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy

      tmpdz = adz(:)*dz
!      tmpw = 0.5*(w(i,j,1:nzm) + w(i,j,2:nz))  ! MK: changed for stretched grids 
      tmpw = ((zi(2:nz)-z(1:nzm))*w(i,j,1:nzm)+ &
             (z(1:nzm)-zi(1:nzm))*w(i,j,2:nz))/(zi(2:nz)-zi(1:nzm))
      tmpwsub = 0.

      tmppres(:) = 100.*pres(1:nzm)

      if(dototalwater) then
        !bloss/qt: saturation adjustment to compute cloud liquid water content.
        !          Note: tmpqv holds qv+qcl on input, qv on output.
        !                tmptabs hold T-(L/Cp)*qcl on input, T on output.
        !                tmpqcl hold qcl on output.
        !                tmppres is unchanged on output, should be in Pa.
        call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
      end if

      i1 = 1 ! dummy variables used by WRF convention in subroutine call
      i2 = 1
      j1 = 1
      j2 = 1

      mtendqv = 0.
      mtendqcl = 0.
      mtendqr = 0.
      mtendqad = 0.
      mtendqaw = 0.
      mtendqar = 0.
      mtendqci = 0.
      mtendqs = 0.
      mtendqg = 0.
      mtendncl = 0.
      mtendnr = 0.
      mtendnad = 0.
      mtendnci = 0.
      mtendns = 0.
      mtendng = 0.

      tmtend1d = 0.

      sfcpcp = 0.
      sfcicepcp = 0.

      effc1d(:) = 10. ! default liquid and ice effective radii
      effr1d(:) = 30. ! default rain effective radii, this is a swag as this value is calculated, anyhow
      effi1d(:) = 75.

      if(do_mkbudget_extra) proc_extra(:,:) = 0.

      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/o sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           mtendqcl,mtendqci,mtendqs,mtendqr,mtendqad,mtendqaw,mtendqar,&
           mtendncl,mtendnci,mtendns,mtendnr,mtendnad,          &
           tmpqcl,tmpqci,tmpqs,tmpqr,tmpqad,tmpqaw,tmpqar, &
           tmpncl,tmpnci,tmpns,tmpnr,tmpnad,tmpnuc,tmpnur,        &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho,tmpdz,tmpw,tmpwsub, &
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
           stendqr,stendqad,stendnad,stendqaw,stendqar,stendqci,stendqs,stendqcl, &
           tmp_cl_pgam, tmp_cl_lambda, &
           micro_proc_rates,do_accumulate_process_rates, &
           nmicro_fields, stend1d, & !bloss(2020-11): Make sedimentation tendencies available for all species
           proc_extra, n_mkbudget_extra, do_mkbudget_extra)
 
     ! update microphysical quantities in this grid column
      if(doprecip) then
         total_water_prec = total_water_prec + sfcpcp

         ! take care of surface precipitation
         precsfc(i,j) = precsfc(i,j) + sfcpcp/dz
         prec_xy(i,j) = prec_xy(i,j) + sfcpcp/dz

         ! update rain
         micro_field(i,j,:,iqr) = tmpqr(:)
         micro_field(i,j,:,inr) = tmpnr(:)

         if(do_chunked_energy_budgets) then
           prec_accum(i,j) = prec_accum(i,j) + sfcpcp/dz
           qtot_sed(i,j,:) = qtot_sed(i,j,:) &
                + dtn*( stendqcl(:) + stendqr(:)) 

           if(doicemicro) then
             prec_ice_accum(i,j) = prec_ice_accum(i,j) + sfcicepcp/dz
             qtot_sed(i,j,:) = qtot_sed(i,j,:) &
                  + dtn*( stendqci(:) + stendqs(:) + stendqg(:) )
             qice_sed(i,j,:) = qice_sed(i,j,:) &
                  + dtn*( stendqci(:) + stendqs(:) + stendqg(:) )
           end if

          if(do_chunk_mkbudget) then
            ! store the sedimentation tendencies in mkbudget_sed
            do n = 1,n_mkbudget
              ! standard water
              if(flag_mkbudget(iqv,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqv)
              
              if(iqcl.gt.0) then
                if(flag_mkbudget(iqcl,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqcl)
              end if

              if(iqr.gt.0) then
                if(flag_mkbudget(iqr,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqr)
              end if

              if(incl.gt.0) then
                if(flag_mkbudget(incl,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,incl)
              end if

              if(inr.gt.0) then
                if(flag_mkbudget(inr,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,inr)
              end if

              if(doicemicro) then
                if(flag_mkbudget(iqci,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqci)
                if(flag_mkbudget(iqs,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqs)
                if(flag_mkbudget(iqg,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqg)
              end if

              if(doprogaerosol) then
                if(flag_mkbudget(inacc,n).gt.0) then
                  mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,inad) + dtn*stend1d(:,incl)
                end if
                if(flag_mkbudget(iqad,n).gt.0) then
;                  mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqad) + dtn*stend1d(:,iqaw)
                end if
                if(flag_mkbudget(iqar,n).gt.0) mkbudget_sed(i,j,:,n) = mkbudget_sed(i,j,:,n) + dtn*stend1d(:,iqar)
              end if
            end do

            if(do_mkbudget_extra) then
              ! additional microphysical tendencies for output in chunk-averaged budgets
              do n = 1,n_mkbudget_extra
                do k = 1,nzm
                  mkbudget_extra(i,j,k,n) =  mkbudget_extra(i,j,k,n) + dtn*proc_extra(k,n)
                end do
              end do
            end if
            
          end if
        end if
      else
         ! add rain to cloud
         tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
         tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

         ! zero out rain 
         tmpqr(:) = 0.
         tmpnr(:) = 0.

         ! add rain tendencies to cloud
         stendqcl(:) = stendqcl(:) + stendqr(:)
         mtendqcl(:) = mtendqcl(:) + mtendqr(:)
         mtendncl(:) = mtendncl(:) + mtendnr(:)

         ! zero out rain tendencies
         stendqr(:) = 0.
         mtendqr(:) = 0.
         mtendnr(:) = 0.
      end if

      if(dototalwater) then
        !bloss/qt: update total water and cloud liquid.
        !          Note: update of total water moved to after if(doprecip),
        !                  since no precip moves rain --> cloud liq.
        micro_field(i,j,:,iqv) = tmpqv(:) + tmpqcl(:) !bloss/qt: total water
      else
        micro_field(i,j,:,iqv) = tmpqv(:) ! water vapor
        micro_field(i,j,:,iqcl) = tmpqcl(:) ! cloud liquid water mass mixing ratio
      end if
      cloudliq(i,j,:) = tmpqcl(:) !bloss/qt: auxilliary cloud liquid water variable
      if(dopredictNc) micro_field(i,j,:,incl) = tmpncl(:)

      reffc(i,j,:) = effc1d(:)
      CloudLiquidMassMixingRatio(i,j,:) = tmpqcl(:)
      CloudLiquidGammaExponent(i,j,:) = tmp_cl_pgam(:)
      CloudLiquidLambda(i,j,:) = tmp_cl_lambda(:)
      reffr(i,j,:) = effr1d(:)

      if(doprogaerosol) then
        !bloss(2020-11): prognostic variablses for wet+dry aerosol
        tmpqacc(:) = tmpqad(:)+tmpqaw(:)
        tmpnacc(:) = tmpnad(:)+tmpncl(:)

        !bloss(2020-11): Partition total (dry+wet) aerosol mass into
         !   dry aerosol mass and wet (in-cloud-droplet) aerosol mass
         call PartitionAerosolMass( nzm, tmpnacc, tmpqacc, aer_sig1, &
              tmpncl, tmpnad, tmpqad, tmpqaw)
        
        !bloss(2020-11): update prognostic aerosol variables
        micro_field(i,j,:,iqacc) = tmpqacc(:)*massfactor
        micro_field(i,j,:,inacc) = tmpnacc(:)
        micro_field(i,j,:,incl) = tmpncl(:) ! in case ncl was limited to nacc inside PartitionAerosolMass

        !bloss(2020-11): update diagnostic aerosol variables
        micro_field(i,j,:,inad) = tmpnad(:) ! dry aerosol number
        micro_field(i,j,:,iqad) = tmpqad(:)*massfactor ! dry aerosol mass
        micro_field(i,j,:,iqaw) = tmpqaw(:)*massfactor ! wet (in-cloud-droplet) aerosol mass

        if(.NOT.doprecoff) micro_field(i,j,:,iqar) = tmpqar(:)*massfactor ! in-rain-drop aerosol mass


        if((.NOT.shiptrack_timeflag).AND.(time.GE.shiptrack_time0)) then
          call task_rank_to_index(rank,it,jt)
          if((it.LE.(nx_gl/2)).AND.((it+nx).GE.(nx_gl/2))) then
!bloss(2018-02): RHOA --> rho_aerosol
!bloss(2019-07): pi/6 --> 4.*pi/3. because mass ~ 4/3*pi*r^3
!!$            dumM1 = aer_ship_n*RHOA*pi/6.*(aer_rm1**3)*EXP(9.*(log(aer_sig1)**2)/2.)
            dumM1 = aer_ship_n*rho_aerosol*4.*pi/3.*(aer_rm1**3)*EXP(9.*(log(aer_sig1)**2)/2.)
            do k = 1,nzm
              if (zi(k+1).LT.100) then
                micro_field(((nx_gl/2)-it+1),:,k,iqacc) = dumM1
                micro_field(((nx_gl/2)-it+1),:,k,inacc) = aer_ship_n
              end if
            end do
            shiptrack_timeflag = .true.
          end if
        end if
      end if

      if(doicemicro) then
         micro_field(i,j,:,iqci) = tmpqci(:)
         micro_field(i,j,:,inci) = tmpnci(:)
         micro_field(i,j,:,iqs) = tmpqs(:)
         micro_field(i,j,:,ins) = tmpns(:)
         if(dograupel) then
            micro_field(i,j,:,iqg) = tmpqg(:)
            micro_field(i,j,:,ing) = tmpng(:)
         end if
         reffi(i,j,:) = effi1d(:)  
         reffs(i,j,:) = effs1d(:)  
         CloudIceMassMixingRatio(i,j,:) = tmpqci(:)
         SnowMassMixingRatio(i,j,:) = tmpqs(:)
      end if

      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================
      
      !=====================================================
      ! Call M2011 scavenging if applicable
      if(do_m2011_scavenge.AND.doprogaerosol) then
          
          call t_startf ('micro_scav')
          dum = 0.
          dum2 = 0.
          dum3 = 0.
          dum4 = 0.
          dum5 = 0.

          !bloss(2020-10): Compute scavenging tendencies at every time step.
          !bloss          if(((mod(nstep,nscv).EQ.0).AND.(icycle.EQ.1)).OR.((nstep.EQ.1).AND.(icycle.EQ.1))) then
           
            aer_sig_arr(1) = aer_sig1
            scvtndqadcl(:) = 0.
            scvtndqadr(:) = 0.
            scvtndnadcl(:) = 0.
            scvtndnadr(:) = 0.

            !bloss(TODO): Convert the m2011_scavenging() subroutine
            !  into two pure functions that would each operate on arguments
            !   at a single grid point: m2011_scavenging_cloud() and
            !   m2011_scavenging_rain() which return vectors of
            !   length 2 with dn/dt and dq/dt for the dry aerosol.
            !  That way, we could loop over k and only populate the
            !  arguments when qcl(k)>qxeps for cloud or qr(k)>qxeps
            !  for rain.

            if(dototalwater) then
              tmpqcl(:) = cloudliq(i,j,:)
            else
              tmpqcl(:) = micro_field(i,j,:,iqcl)
            end if
            if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)

            ! get aerosol number masses
            tmpnad(:) = micro_field(i,j,:,inad)
            tmpqad(:) = micro_field(i,j,:,iqad)
            tmpqaw(:) = micro_field(i,j,:,iqaw)

            if(doprecip) then
              tmpqr(:) = micro_field(i,j,:,iqr)
              tmpnr(:) = micro_field(i,j,:,inr)
            end if

            if(doprecip.AND.doprogaerosol) then
              tmpqar(:) = micro_field(i,j,:,iqar)/massfactor
            end if
           
            tmpqv(:) = micro_field(i,j,:,iqv)
            tmptabs(:) = t(i,j,:)  &           ! liquid water-ice static energy over Cp
                 - gamaz(:) &                                   ! potential energy
                 + fac_cond * (tmpqr(:) + tmpqcl(:))     ! liquid latent energy         

            if(doicemicro) then
              tmptabs(:) = tmptabs(:)  &           ! liquid water-ice static energy over Cp
                   + fac_sub * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy         
            end if

            tmppres(:) = 100.*pres(1:nzm)
            tmpzero(:) = 0.
            
            if(.NOT.isinitializedSCAV) then
              ! This code will be invoked at the start of a simulation, either
              !  fresh or restarted.  Using this at restart may cause restarted runs
              !  to not be bit-for-bit identical to non-restarted runs.  I view this
              !  as acceptable since we may save a lot of computation time by not calling
              !  the scavenging at every step.
              iscv = 1
              isinitializedSCAV = .true.
            end if

            if(isinitializedSCAV) then
              if(icycle.eq.1) iscv = iscv - 1
              if((nscv.eq.0).OR.(iscv.eq.0)) then
                !bloss(2020-11): Only compute scavenging for cloud every nscv time steps.
                !    If nscv=0, compute it every substep (icycle) of every timestep
                !    If nscv=1, compute it only for the first substep (icycle==1) of every timestep
                !    If nscv>1, compute it only for the first substep (icycle==1) of every nscv_th timestep
                !    Compute it every time for rain since rain falls faster and is more transient.

                ! set the hygroscopicity from value in micro_params.f90
                hygro_arr(1) = hygro


                do k = 1,nzm
                  ! compute RH wrt liquid
                  tmprhw(k) = tmpqv(k)/qsatw(tmptabs(k),pres(k))
                  if(DoFixInCloudSuperSatForScavenging) then
                    if ((tmprhw(k).GT.0.99).AND.(tmpqcl(k).GT.1.e-5)) then
!!$                    if ((tmprhw(k).GT.0.99).AND.(tmpqcl(k).GT.0.0005)) then
                      tmprhw(k) = 1. + InCloudSuperSatPercentForScavenging/100.
                      !bloss                    tmprhw(k)= 1.005 !in cloud RH hack
                    end if
                  else
                    write(*,*) 'Predicted in-cloud supersaturation is not enabled at this time.'
                    call task_abort()
                  end if

                  !bloss(2019-07-22): check for qad>0 and Nad>0
                  if ((tmpnad(k).gt.0.).AND.(tmpqad(k).gt.0.)) then
                    !bloss(2019-07-22): NOTE: exponential factor already includes 1/3 power
                    ddry = EXP(-1.5*LOG(aer_sig1)**2.)*(6.*tmpqad(k)/&
                         (pi*rho_aerosol*tmpnad(k)))**(1./3.)

                    ! code from Matt Wyant for computing wet radius
                    !   Note: aerosol swells to sizes bigger than dry diameter as RH --> 1.
                    relhum_arr(1) = tmprhw(k)
                    rdry(1) = ddry/2.
                  
                    call modal_aero_kohler(rdry, hygro_arr, relhum_arr, rwet(:), npoints)
                    if (rwet(1).lt.5.e-10) then
                      rwet(1) = 5.e-10
                    end if
                    tmpcmd(k) = rwet(1)*2.
                    tmpcmddry(k) = ddry
                    if((nstep.eq.1).AND.(icycle.eq.1).AND.masterproc) then
                      write(*,942) k, tmpqcl(k), tmprhw(k), tmpcmddry(k), tmpcmd(k)
                      942 format('k/QCL/RHW/DryDiam/WetDiam = ',I4,4e12.4)
                    end if
                  else
                    tmpcmd(k) = 1.e-9 ! min value from Matt Wyant
                    tmpcmddry(k) = 1.e-9
                  end if
                end do

                !bloss: Only use this option for testing to determine how much the
                !  use of wet diameter affects scavenging.
                if(DoScavengeAerosolDryDiameter) tmpcmd(:) = tmpcmddry(:)

                ! compute scavenging due to cloud
                if (doscavcloud2m.eq..true.) then
                  call scav_cloud_2m(tmpnad,tmpqad,tmpcmd,aer_sig_arr,&
                       tmptabs,tmppres,tmpqcl,tmpncl,tmpnuc,tmpec3d,dodissip3d,&
                       1,1,nzm,1,dtn)
                else
                  if(MAXVAL(tmpqcl(1:nzm)).gt.qxeps) then
                    !bloss(2020-10): Only call scavenging routine if cloud is present.
                    !From MWyant: Add dry diameter as an argument.
                    call m2011_scavenging(tmpnad,tmpqad,tmpcmd,tmpcmddry,aer_sig_arr,&
                         tmptabs,tmppres,tmpqv,tmprhw,tmpzero,&
                         tmpqcl,tmpzero,tmpzero,tmpzero,tmpzero,tmpncl,tmpzero,tmpzero,&
                         tmpzero,tmpzero,.true.,.false.,.false.,.false.,.false.,&
                         1,1,nzm,1,dtn)
                  end if !if(MAX(QCL)>1e-5 kg/kg)
                end if
                scvtndqadcl(:) = (micro_field(i,j,:,iqad)/massfactor-tmpqad(:))/dtn !positive value
                scvtndnadcl(:) = (micro_field(i,j,:,inad)-tmpnad(:))/dtn            !positive value

                iscv = nscv ! reset counter

              else ! in this case, iscv > 0
                ! Use cloud scavenging tendencies from last computed value
                scvtndqadcl(1:nzm) = scav3d(i,j,1:nzm,2)
                scvtndnadcl(1:nzm) = scav3d(i,j,1:nzm,4)

              end if !if((iscv==0)|(nscv==0))

            end if !if(isinitializedSCAV)

      ! compute scavenging due to rain
            tmpqad(:) = micro_field(i,j,:,iqad)/massfactor
            tmpnad(:) = micro_field(i,j,:,inad)
            
            !bloss(2020-10): Only call scavenging routine if rain is present in the column
            if(MAXVAL(tmpqr(1:nzm)).gt.qxeps) then

            do k = 1, nzm
              !bloss(2019-07-22): check for qad>0 and Nad>0
              if ((tmpnad(k).gt.0.).AND.(tmpqad(k).gt.0.)) then
                tmpcmd(k) = EXP(-1.5*LOG(aer_sig1)**2.)*(6.*tmpqad(k)/&
                     (pi*rho_aerosol*tmpnad(k)))**(1./3.)
              else 
                tmpcmd(k) = 0.
              end if
            end do

            !From MWyant: Add dry diameter as an argument.
            call m2011_scavenging(tmpnad,tmpqad,tmpcmd,tmpcmddry, aer_sig_arr,&
                 tmptabs,tmppres,tmpqv,tmpzero,tmpzero,&
                 tmpzero,tmpqr,tmpzero,tmpzero,tmpzero,tmpzero,tmpnr,tmpzero,&
                 tmpzero,tmpzero,.false.,.true.,.false.,.false.,.false.,&
                 1,1,nzm,1,dtn)
        
          end if !if(MAX(QR)>1e-5 kg/kg)

          scvtndqadr(:) = (micro_field(i,j,:,iqad)/massfactor-tmpqad(:))/dtn !positive value
          scvtndnadr(:) = (micro_field(i,j,:,inad)-tmpnad(:))/dtn            !positive value

!bloss: Call scavenging at every step        
!bloss        end if

          ! refresh aerosol number/mass
          tmpnacc(:) = micro_field(i,j,:,inacc)
          tmpqacc(:) = micro_field(i,j,:,iqacc)
          tmpnad(:) = micro_field(i,j,:,inad)
          tmpqad(:) = micro_field(i,j,:,iqad)
          tmpqar(:) = micro_field(i,j,:,iqar)

        !bloss(2020-10): Only limit tendencies when they are applied.
        do k=1,nzm

          ! note that all of these tendencies should be non-negative
          !   Make sure of that so that the limiting does not divide by zero unintentionally.
          scav_mass = MAX(0., (scvtndqadr(k) + scvtndqadcl(k))*dtn )
          scav_number = MAX(0., (scvtndnadr(k) + scvtndnadcl(k))*dtn )

          dry_aerosol_mass_before = MAX(0., micro_field(i,j,k,iqad) )
          dry_aerosol_number_before = MAX(0., micro_field(i,j,k,inad) )

          !bloss: Restructure limiting 
          if ( (scav_mass*massfactor.gt.dry_aerosol_mass_before) &
               .OR. (scav_number.gt.dry_aerosol_number_before) ) THEN

            scav_factor = MIN(1., &
                 dry_aerosol_mass_before/(scav_mass*massfactor), &
                 dry_aerosol_number_before/scav_number ) 

            scvtndqadr(k)  = scvtndqadr(k)*scav_factor
            scvtndqadcl(k) = scvtndqadcl(k)*scav_factor
            scvtndnadr(k)  = scvtndnadr(k)*scav_factor
            scvtndnadcl(k) = scvtndnadcl(k)*scav_factor
          end if

          ! update aerosol mass and number
          ! note: we do not need to update dry and wet aerosol quantities
          !   separately, because we will diagnose these from the total
          !   aerosol mass and the fraction of dry aerosols after scavenging.

          ! dry+wet aerosol mass is only affected by rain scavenging
          tmpqacc(k) = MAX(0., tmpqacc(k) - dtn*scvtndqadr(k) )
          tmpqar(k)  = MAX(0., tmpqar(k)  + dtn*scvtndqadr(k) )

          ! update aerosol number
          ! dry+wet aerosol number is affected by both cloud and rain scavenging
          tmpnacc(k) = MAX(0., tmpnacc(k) - dtn*scvtndnadcl(k)  - dtn*scvtndnadr(k) )
        end do

        ! Now, since we have changed the number and mass of wet and dry aerosol,
        !   we need to re-partition the aerosol mass among dry and wet aerosols.
        !   This is necessary to be consistent with the assumption that the 
        !   lognormal size distribution of accumulation mode aerosol is divided
        !   into smaller dry aerosols and larger wet aerosols.
        call PartitionAerosolMass( nzm, tmpnacc, tmpqacc, aer_sig1, &
             tmpncl, tmpnad, tmpqad, tmpqaw)

        !bloss(2020-11): update prognostic aerosol variables
        micro_field(i,j,:,iqacc) = tmpqacc(:)*massfactor
        micro_field(i,j,:,inacc) = tmpnacc(:)
        micro_field(i,j,:,incl) = tmpncl(:) ! in case ncl was limited to nacc inside PartitionAerosolMass

        !bloss(2020-11): update diagnostic aerosol variables
        micro_field(i,j,:,inad) = tmpnad(:) ! dry aerosol number
        micro_field(i,j,:,iqad) = tmpqad(:)*massfactor ! dry aerosol mass
        micro_field(i,j,:,iqaw) = tmpqaw(:)*massfactor ! wet (in-cloud-droplet) aerosol mass

        if(dostatis) then
          ! accumulate tendencies so that statistics will
          !   hold horizonally-averaged tendencies.
          scvtndqadrstat(:) = scvtndqadrstat(:) + scvtndqadr(:)
          scvtndqadclstat(:) = scvtndqadclstat(:) + scvtndqadcl(:)
          scvtndnadrstat(:) = scvtndnadrstat(:) + scvtndnadr(:)
          scvtndnadclstat(:) = scvtndnadclstat(:) + scvtndnadcl(:)
        end if

        if(do_mkbudget_extra) then
          ! scavenging tendencies for aerosol mass
          do k = 1,nzm
            mkbudget_extra(i,j,k,n_mkbudget_extra-3) =  &
                 mkbudget_extra(i,j,k,n_mkbudget_extra-3) + dtn*scvtndqadr(k)
            mkbudget_extra(i,j,k,n_mkbudget_extra-2) =  &
                 mkbudget_extra(i,j,k,n_mkbudget_extra-2) + dtn*scvtndqadcl(k)

            mkbudget_extra(i,j,k,n_mkbudget_extra-1) = &
                 mkbudget_extra(i,j,k,n_mkbudget_extra-1) + dtn*scvtndnadr(k)
            mkbudget_extra(i,j,k,n_mkbudget_extra)   =  &
                 mkbudget_extra(i,j,k,n_mkbudget_extra)   + dtn*scvtndnadcl(k)
          end do
        end if

        if(do_scav_3d_output) then
          !bloss(2020-10): store scavenging tendencies in 3D array 
          !  in case we need to use them again, since scavenging is not
          !  called every icycle
          scav3d(i,j,1:nzm,1) = scvtndqadr(1:nzm)
          scav3d(i,j,1:nzm,2) = scvtndqadcl(1:nzm)
          scav3d(i,j,1:nzm,3) = scvtndnadr(1:nzm)
          scav3d(i,j,1:nzm,4) = scvtndnadcl(1:nzm)
        end if

         call t_stopf ('micro_scav')
      end if !do_m2011_scavenge.AND.doprogaerosol
      
      if(dostatis) then
        if(dototalwater) then
          !bloss/qt: total water microphysical tendency includes qv and qcl
          mtend(:,iqv) = mtend(:,iqv) + mtendqv + mtendqcl
        else
          ! separate tendencies for vapor and cloud liquid mass
          mtend(:,iqv) = mtend(:,iqv) + mtendqv 
          mtend(:,iqcl) = mtend(:,iqcl) + mtendqcl
        end if

         if(dopredictNc) mtend(:,incl) = mtend(:,incl) + mtendncl
         if(doprecip) then
            mtend(:,iqr) = mtend(:,iqr) + mtendqr
            mtend(:,inr) = mtend(:,inr) + mtendnr
         end if
         
         if(doprogaerosol) then
           !bloss(2020-11): Combined dry+wet aerosol
           mtend(:,iqacc) = mtend(:,iqacc) + mtendqad + mtendqaw
           mtend(:,inacc) = mtend(:,inacc) + mtendnad + mtendncl

            mtend(:,iqad) = mtend(:,iqad) + mtendqad
            mtend(:,iqaw) = mtend(:,iqaw) + mtendqaw
            mtend(:,inad) = mtend(:,inad) + mtendnad
            if(.NOT.doprecoff) mtend(:,iqar) = mtend(:,iqar) + mtendqar
         endif

         if(doicemicro) then
            mtend(:,iqci) = mtend(:,iqci) + mtendqci
            mtend(:,inci) = mtend(:,inci) + mtendnci

            mtend(:,iqs) = mtend(:,iqs) + mtendqs
            mtend(:,ins) = mtend(:,ins) + mtendns

            if(dograupel) then
               mtend(:,iqg) = mtend(:,iqg) + mtendqg
               mtend(:,ing) = mtend(:,ing) + mtendng
            end if
         end if

         do n = 1,nmicro_fields
            do k = 1,nzm
               if(micro_field(i,j,k,n).ge.1.e-6) mfrac(k,n) = mfrac(k,n)+1.
            end do
         end do

!!$         if (doprogaerosol.AND.do_m2011_scavenge.AND.do_micro3Dout) then
!!$            do n = 1,30
!!$               mtend3d(i,j,:,n) = mtendaux(:,n)
!!$            end do
!!$            mtend3d(i,j,:,31) = scvtndqadclstat(:)
!!$            mtend3d(i,j,:,32) = scvtndqadrstat(:)
!!$            mtend3d(i,j,:,33) = scvtndnadclstat(:)
!!$            mtend3d(i,j,:,34) = scvtndnadrstat(:)
!!$         end if

         !bloss(2018-02): This should be 1.5, rather than 1.8.
         ! approximate optical depth = 0.0015*lwp/effrad
         !  integrated up to level at which output
         tmpc = 0.
         tmpr = 0.
         tmpi = 0.
         tmps = 0.
         tmpg = 0.

         do k = 1,nzm
            tmpc = tmpc + 0.0015*rho(k)*dz*adz(k)*tmpqcl(k)/(1.e-20+1.e-6*effc1d(k))
            tmpr = tmpr + 0.0015*rho(k)*dz*adz(k)*tmpqr(k)/(1.e-20+1.e-6*effr1d(k))
            !bloss/qt: put cloud liquid optical depth in trtau(:,iqv)
            trtau(k,iqv) = trtau(k,iqv) + tmpc
            if(doprecip) trtau(k,iqr) = trtau(k,iqr) + tmpr

            if(doicemicro) then
               tmpi = tmpi + 0.0015*rho(k)*dz*adz(k)*tmpqci(k)/(1.e-20+1.e-6*effi1d(k))
               tmps = tmps + 0.0015*rho(k)*dz*adz(k)*tmpqs(k)/(1.e-20+1.e-6*effs1d(k))

               trtau(k,iqci) = trtau(k,iqci) + tmpi
               trtau(k,iqs) = trtau(k,iqs) + tmps

               if(dograupel) then
                  tmpg = tmpg + 0.0015*rho(k)*dz*adz(k)*tmpqg(k)/(1.e-20+1.e-6*effg1d(k))
                  trtau(k,iqg) = trtau(k,iqg) + tmpg
               end if
            end if
         end do


         tlat(1:nzm) = tlat(1:nzm) &
              - dtn*fac_cond*(stendqcl+stendqr) &
              - dtn*fac_sub*(stendqci+stendqs+stendqg)
         qpfall(1:nzm) = qpfall(1:nzm) + dtn*(stendqr+stendqs+stendqg)

         !bloss: temperature tendency (sensible heating) due to phase changes
         tmtend3d(i,j,1:nzm) = tmtend1d(1:nzm)

      end if ! dostatis

       if(doreflectivity_cloudradar.AND. &
            (dostatis_quickbeam.OR. (mod(nstep,nsave3D).eq.0.AND.icycle.eq.ncycle)) ) then

         call t_startf ('micro_quickbeam')

         hgt_matrix(i,1:nzm) = z(1:nzm)*1.e-3 ! in km

         hm_matrix(1,i,1:nzm) = tmpqcl(1:nzm)*1e3 ! mixing ratio in g/kg
         hm_matrix(2,i,1:nzm) = tmpqr(1:nzm)*1e3

         re_matrix(1,i,1:nzm) = effc1d(1:nzm) ! in microns
         re_matrix(2,i,1:nzm) = effr1d(1:nzm)

         if(doicemicro) then
           hm_matrix(3,i,1:nzm) = tmpqci(1:nzm)*1e3 ! cloud ice
           hm_matrix(4,i,1:nzm) = tmpqs(1:nzm)*1e3 ! snow
           hm_matrix(5,i,1:nzm) = tmpqg(1:nzm)*1e3 ! graupel

           re_matrix(3,i,1:nzm) = effi1d(1:nzm)
           re_matrix(4,i,1:nzm) = effs1d(1:nzm)
           re_matrix(5,i,1:nzm) = effg1d(1:nzm)

         end if

         ! use Hugh's effective radii instead of droplet number
         Np_matrix(:,:,:) = -1.

         p_matrix(i,1:nzm) = pres(1:nzm) ! in hPa
         t_matrix(i,1:nzm) = tmptabs(1:nzm) ! in deg K

         do k = 1,nzm
           esat1 = MIN(0.99*pres(k), polysvp(tmptabs(k),0) ) ! esat is in Pa, with fix from Hugh for low pressure
           qsat1 = 0.622*esat1/ (100.*pres(k) - esat1)
           rh_matrix(i,k) = MIN(100., 100.*tmpqv(k)/qsat1)
         end do

         call t_stopf ('micro_quickbeam')

       end if ! if(doreflectivity_cloudradar.AND.dostatis_quickbeam) 
      if(dototalwater) then
        ! since iqv includes both vapor and cloud liquid, add qcl sedimentation tendency here.
        stend(:,iqv) = stend(:,iqv) + stendqcl !bloss/qt: iqcl --> iqv
      else
        ! since cloud liquid is separate here, add qcl sedimentation tendency here.
        stend(:,iqcl) = stend(:,iqcl) + stendqcl 
      end if

      if(doprecip) then
         stend(:,iqr) = stend(:,iqr) + stendqr
      end if
      
      if(doprogaerosol) then
        !bloss(2020-11): Save sedimentation tendencies for wet+dry
        !                  aerosol number/mass and rain aerosol mass
         stend(:,iqacc) = stend(:,iqacc) + stend1d(:,iqad) + stend1d(:,iqaw) !bloss(2020-11)
         stend(:,inacc) = stend(:,inacc) + stend1d(:,inad) + stend1d(:,incl) !bloss(2020-11)
         if(.NOT.doprecoff) stend(:,iqar) = stend(:,iqar) + stend1d(:,iqar) !bloss(2020-11)

!bloss         stend(:,iqad) = stend(:,iqad) + stendqad
!bloss         stend(:,iqaw) = stend(:,iqaw) + stendqaw
!bloss         stend(:,inad) = stend(:,inad) + stendnad
      end if

      if(doicemicro) then
         stend(:,iqci) = stend(:,iqci) + stendqci
         stend(:,inci) = stend(:,inci) + stend1d(:,inci) !bloss(2020-11)
         stend(:,iqs) = stend(:,iqs) + stendqs
         stend(:,ins) = stend(:,ins) + stend1d(:,ins) !bloss(2020-11)
         if(dograupel) then
           stend(:,iqg) = stend(:,iqg) + stendqg
           stend(:,ing) = stend(:,ing) + stend1d(:,ing) !bloss(2020-11)
         end if
      end if

   end do ! i = 1,nx

   if(doreflectivity_cloudradar.AND. &
        (dostatis_quickbeam .OR. (mod(nstep,nsave3D).eq.0.AND.icycle.eq.ncycle)) ) then

     call t_startf ('micro_quickbeam')

     !bloss: cloud radar reflectivity computation.
     !       Call once per row to allow some vectorization.
     nprof = nx
     ngate = nzm
     if(masterproc) then
       write(*,*) 'Calling Radar Simulator'
       write(*,*) 'Max/min(Np_matrix) = ', MAXVAL(Np_matrix), MINVAL(Np_Matrix)
     end if
     call radar_simulator( &
          hp_cloudradar, & ! structure that holds radar parameters, description of drop/particle size distn
          nprof,ngate, & ! # of columns, # of levels
          missing_value_cloudradar, & ! like it sounds
          hgt_matrix, hm_matrix, re_matrix, Np_matrix, &
          p_matrix, t_matrix, rh_matrix, &
          Ze_non,Ze_ray,h_atten_to_vol,g_atten_to_vol,dBZe, &
          g_to_vol_in, g_to_vol_out)

     do k = 1,nzm
       dBZ_cloudradar(1:nx,j,k) = dBZe(1:nx,k)
     end do

     call t_stopf ('micro_quickbeam')

   end if ! if(doreflectivity_cloudradar.AND.dostatis)

end do ! j = 1,ny

! back sedimentation flux out from sedimentation tendencies
if(dototalwater) then
   tmpc = 0.
   do k = 1,nzm
      m = nz-k
      tmpc = tmpc + stend(m,iqv)*rho(m)*dz*adz(m)  !bloss/qt: iqcl --> iqv
      mksed(m,iqv) = tmpc
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqv)*dtn/dz
else
   !bloss(2019-07): We need to include this so that the cloud droplet
   !  sedimentation will be accounted for when dototalwater==.false.
   tmpc = 0.
   do k = 1,nzm
      m = nz-k
      tmpc = tmpc + stend(m,iqcl)*rho(m)*dz*adz(m)
      mksed(m,iqcl) = tmpc
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqcl)*dtn/dz
end if

if(doprecip) then
   tmpr = 0.
   do k = 1,nzm
      m = nz-k
      tmpr = tmpr + stend(m,iqr)*rho(m)*dz*adz(m)
      mksed(m,iqr) = tmpr
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqr)*dtn/dz
end if

if(doicemicro) then
   tmpi = 0.
   tmps = 0.
   tmpg = 0.
   do k = 1,nzm
      m = nz-k
      tmpi = tmpi + stend(m,iqci)*rho(m)*dz*adz(m)
      tmps = tmps + stend(m,iqs)*rho(m)*dz*adz(m)
      mksed(m,iqci) = tmpi
      mksed(m,iqs) = tmps
      if(dograupel) then
        tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
        mksed(m,iqg) = tmpg
      end if
   end do
   precflux(1:nzm) = precflux(1:nzm) &
        - (mksed(:,iqci) + mksed(:,iqs))*dtn/dz
   if(dograupel) precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqg)*dtn/dz
end if

if (any(micro_field.ne.micro_field)) loc_nan = 1

if(dompi) then
   call task_max_integer(loc_nan,glb_nan,1)
else
   glb_nan = loc_nan
end if

if (glb_nan(1).eq.1) then
  print *,"NaN detected!!!"
  call write_fields3D()
  call task_abort()
end if

!!$if(doprecip) total_water_prec = total_water_prec - total_water()

if (docloud)  call micro_diagnose()   ! leave this line here

call t_stopf ('micro_proc')

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
!  This is the pace where the microphysics field that SAM actually cares about
!  are diagnosed.

subroutine micro_diagnose()

real omn, omp, tmp(nzm,nmicro_fields)
integer i,j,k,n

if(dototalwater) then
  ! water vapor = total water - cloud liquid
  qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqv) &
       - cloudliq(1:nx,1:ny,1:nzm)

  ! cloud liquid water
  qcl(1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm)
else
  qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqv)
  qcl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqcl)
end if

! rain water
if(doprecip) qpl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqr)

! cloud ice 
if(doicemicro) then
   qci(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqci)

   if(dograupel) then
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs) &
           + micro_field(1:nx,1:ny,1:nzm,iqg)
   else
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs)
   end if
end if

! save dmoain-mean profiles of each microphysical field
do n = 1,nmicro_fields
  do k = 1,nzm
    mk0(k,n) = SUM(micro_field(1:nx,1:ny,k,n))/real(nx*ny)
  end do
end do
if(dompi) then
  call task_sum_real(mk0,tmp,nzm*nmicro_fields)
  mk0 = tmp/real(nsubdomains)
end if

end subroutine micro_diagnose

!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need supply functions to compute terminal velocity for all of your 
! precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and water vapor (single values, not arrays). Var1 and var2 
! are some microphysics variables like water content and concentration.
! Don't change the number of arguments or their meaning!

!!$real function term_vel_qr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_qr
!!$
!!$real function term_vel_Nr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_Nr
!!$
!!$real function term_vel_qs(qs,ns,tabs,rho)
!!$! .......  
!!$end function term_vel_qs

! etc.

!----------------------------------------------------------------------
!!! compute sedimentation 
!
!  The perpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles
subroutine micro_precip_fall()

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for hour hypothetical case, there is no mixed hydrometeor, so omega is not actually used.

integer hydro_type
real omega(nx,ny,nzm) 

integer i,j,k

return ! do not need this routine -- sedimentation done in m2005micro.

!!$! Initialize arrays that accumulate surface precipitation flux
!!$
!!$ if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
!!$   do j=1,ny
!!$    do i=1,nx
!!$     precsfc(i,j)=0.
!!$    end do
!!$   end do
!!$   do k=1,nzm
!!$    precflux(k) = 0.
!!$   end do
!!$ end if
!!$
!!$ do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
!!$    qpfall(k)=0.
!!$    tlat(k) = 0.
!!$ end do
!!$   
!!$! Compute sedimentation of falling variables:
!!$
!!$ hydro_type=0
!!$ call precip_fall(qr, term_vel_qr, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Nr, term_vel_Nr, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qs, term_vel_qs, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ns, term_vel_Ns, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qg, term_vel_qg, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ng, term_vel_Ng, hydro_type, omega)
!!$


end subroutine micro_precip_fall

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
  implicit none
  integer :: k

  ! print out min/max values of all microphysical variables
  do k=1,nmicro_fields
     call fminmax_print(trim(mkname(k))//':', &
          micro_field(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  end do

end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics that will be outputted
!!  to *.stat statistics file

subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,microcount)


character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,microcount, n, ii, jj, ncond

character*8 name
character*80 longname
character*10 units

microcount = 0

name = 'QTFLUX'
longname = 'Total (resolved + subgrid) total water (vapor+cloud) flux'
units = 'W/m2'
call add_to_namelist(count,microcount,name,longname,units,0)

do n = 1,nmicro_fields
  if(dototalwater.OR.(n.ne.iqv)) then
    ! add mean value of microphysical field to statistics
    !   EXCEPT for water vapor (added in statistics.f90)
    name = trim(mkname(n))
    longname = trim(mklongname(n))
    units = trim(mkunits(n))
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if
  if(n.eq.iqv) then
      ! add variance of ONLY total water (vapor + cloud liq) field to statistics
      !   cloud water variance and cloud ice variance
      !   already output in statistics.f90
      name = trim(mkname(n))//'2'
      longname = 'Variance of '//trim(mklongname(n))
      units = '('//trim(mkunits(n))//')^2'
      call add_to_namelist(count,microcount,name,longname,units,0)
   end if

   if(flag_advect(n).eq.1) then
     ! add vertical advective tendency
     name = trim(mkname(n))//'ADV'
     longname = 'Tendency of '//trim(mklongname(n))// &
          ' due to resolved vertical advection'
     units = trim(mkunits(n))//'/day'
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add vertical diffusive tendency
     name = trim(mkname(n))//'DIFF'
     longname = 'Tendency of '//trim(mklongname(n))// &
          ' due to vertical SGS transport'
     units = trim(mkunits(n))//'/day'
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add tendency due to large-scale vertical advection
     name = trim(mkname(n))//'LSADV'
     longname = 'Tendency of '//trim(mklongname(n))// &
          ' due to large-scale vertical advection'
     units = trim(mkunits(n))//'/day'
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add tendency due to microphysical processes
     name = trim(mkname(n))//'MPHY'
     longname = 'Tendency of '//trim(mklongname(n))// &
          ' due to sedimentation and microphysical processes'
     units = trim(mkunits(n))//'/day'
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add vertical diffusive tendency
     name = trim(mkname(n))//'SED'
     longname = 'Tendency of '//trim(mklongname(n))//' due to sedimentation'
     units = trim(mkunits(n))//'/day'
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add storage terms
     name = trim(mkname(n))//'STRG'
     longname = 'Storage of '//trim(mklongname(n))//''
     units = trim(mkunits(n))//'/day'
     call add_to_namelist(count,microcount,name,longname,units,0)

     if(flag_wmass(n).gt.0) then
       ! fluxes output in W/m2 for mass mixing ratios
       units = 'W/m2'
     else
       ! fluxes output in #/m2/s for number concentrations
       units = '#/m2/s'
     end if

     ! add flux of microphysical fields to scalar
     name = trim(mkname(n))//'FLXR'
     longname = 'Resolved flux of '//trim(mklongname(n))
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add subgrid flux of microphysical fields to scalar
     name = trim(mkname(n))//'FLXS'
     longname = 'Subgrid flux of '//trim(mklongname(n))
     call add_to_namelist(count,microcount,name,longname,units,0)

     ! add sedimentation flux of microphysical fields to scalar
     name = trim(mkname(n))//'SDFLX'
     longname = 'Sedimentation flux of '//trim(mklongname(n))
     call add_to_namelist(count,microcount,name,longname,units,0)

   end if

   if((flag_wmass(n).gt.0).and.(n.ne.iqv)) then
      ! add area fraction of microphysical field to statistics
      name = trim(mkname(n))//'FRAC'
      longname = trim(mklongname(n))//' FRACTION'
      units = '1'
      call add_to_namelist(count,microcount,name,longname,units,0)

      ! add approximate optical depth of hydrometeor fields
      name = 'TAU'//trim(mkname(n))
      longname = 'Approx optical depth of '//trim(mklongname(n))
      units = '1'
      call add_to_namelist(count,microcount,name,longname,units,0)

!bloss (Apr 09): Eliminate this output.  It is unreliable when
!            hydrometeor fractions are variable across processors
!            or in time.  You can still compute this from 
!            TAU* and Q* values in the statistics file.
!bloss      ! add approximate optical depth of hydrometeor fields
!bloss      name = 'EFFR'//trim(mkname(n))
!bloss      longname = 'Effective radius of '//trim(mklongname(n))
!bloss      units = 'microns'
!bloss      call add_to_namelist(count,microcount,name,longname,units,0)

      ! add field which can be used to recover mean effective radius.
      name = trim(mkname(n))//'OEFFR'
      longname = 'Mixing ratio of '//trim(mklongname(n)) &
           //' over effective radius, EFFR = ' &
           //trim(mkname(n))//'/'//trim(mkname(n))//'OEFFR'
      units = 'g/kg/microns'
      call add_to_namelist(count,microcount,name,longname,units,0)
   end if

end do

if(dototalwater) then
  !bloss/qt: add output for cloud liquid water (not included explicitly in 
  !  total water formulation).
  call add_to_namelist(count,microcount,'QC', &
       'Cloud liquid water mass mixing ratio', 'g/kg',0)

  ! add approximate optical depth of cloud liquid water
  name = 'TAUQC'
  longname = 'Approx optical depth of cloud liquid water'
  units = '1'
  call add_to_namelist(count,microcount,name,longname,units,0)

  ! add field which can be used to recover mean effective radius.
  name = 'QCOEFFR'
  longname = 'Mixing ratio of QC'// &
       ' over effective radius, EFFR = QC/QCOEFFR'
  units = 'g/kg/microns'
  call add_to_namelist(count,microcount,name,longname,units,0)

  do ncond = 1,ncondavg
    !bloss/qt: add conditional averages for water vapor and cloud liquid water
    call add_to_namelist(count,microcount,'QV' // TRIM(condavgname(ncond)), &
         'Water vapor mixing ratio in ' // TRIM(condavglongname(ncond)),'kg/kg',ncond)
    call add_to_namelist(count,microcount,'QC' // TRIM(condavgname(ncond)), &
         'Cloud liquid water mixing ratio in ' // TRIM(condavglongname(ncond)),'kg/kg',ncond)
  end do

else

  !bloss/qt: Can't be computed reliably in total water formulation.
  ! add temperature tendency (sensible energy) tendency due to mphys
  call add_to_namelist(count,microcount,'QLAT', &
       'Sensible energy tendency due to phase changes', 'K/day',0)

  do ncond = 1,ncondavg
    !bloss/qt: Can't be computed reliably in total water formulation.
    call add_to_namelist(count,microcount,'QLAT' // TRIM(condavgname(ncond)), &
         'Sensible energy tendency due to phase changes in ' // TRIM(condavglongname(ncond)), &
         'K/day',ncond)
  end do

end if

do ncond = 1,ncondavg
  ! add conditional averages of hydrometeor fields
   do n = 1,nmicro_fields
      call add_to_namelist(count,microcount,trim(mkname(n)) // TRIM(condavgname(ncond)), &
           trim(mklongname(n)) // ' in ' // TRIM(condavglongname(ncond)), &
           trim(mkunits(n)),ncond)
   end do
end do

if(do_output_micro_process_rates) then
  ! warm cloud process rates
  do n = 1,nmicro_process_rates_warm_mass
    call add_to_namelist(count,microcount, &
         trim(micro_process_rate_names_warm_mass(n)), &
         trim(micro_process_rate_longnames_warm_mass(n)), &
         'g/kg/day',0)
  end do
  do n = 1,nmicro_process_rates_warm_number
    call add_to_namelist(count,microcount, &
         trim(micro_process_rate_names_warm_number(n)), &
         trim(micro_process_rate_longnames_warm_number(n)), &
         '#/mg/day',0)
  end do

  if(doprogaerosol) then
    ! prognostic aerosol process rates
    do n = 1,nmicro_process_rates_progaer_mass
      call add_to_namelist(count,microcount, &
           trim(micro_process_rate_names_progaer_mass(n)), &
           trim(micro_process_rate_longnames_progaer_mass(n)), &
           'g/kg/day',0)
    end do
    do n = 1,nmicro_process_rates_progaer_number
      call add_to_namelist(count,microcount, &
           trim(micro_process_rate_names_progaer_number(n)), &
           trim(micro_process_rate_longnames_progaer_number(n)), &
           '#/mg/day',0)
    end do
  end if

  if(doicemicro) then
    ! cold cloud process rates
    do n = 1,nmicro_process_rates_cold_mass
      call add_to_namelist(count,microcount, &
           trim(micro_process_rate_names_cold_mass(n)), &
           trim(micro_process_rate_longnames_cold_mass(n)), &
           'g/kg/day',0)
    end do
    do n = 1,nmicro_process_rates_cold_number
      call add_to_namelist(count,microcount, &
           trim(micro_process_rate_names_cold_number(n)), &
           trim(micro_process_rate_longnames_cold_number(n)), &
           '#/mg/day',0)
    end do
  end if

end if

if (doprogaerosol.AND.do_m2011_scavenge) then
  ! set up outputs for scavenging tendencies
  name = 'SCVTQAR'
  longname = 'Aerosol mass tendency due to interestitial scavenging by rain (Source of QAr, Sink of QAd)'
  units = 'g/kg/day'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'SCVTQAW'
  longname = 'Aerosol mass tendency due to interestitial scavenging by cloud (Source of QAw, Sink of QAd)'
  units = 'g/kg/day'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'SCVTNADR'
  longname = 'Aerosol number tendency due to interestitial scavenging by rain (Sink of NAd)'
  units = '#/mg/day'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'SCVTNADC'
  longname = 'Aerosol number tendency due to interestitial scavenging by cloud (Sink of NAd)'
  units = '#/mg/day'
  call add_to_namelist(count,microcount,name,longname,units,0)
end if

if(doreflectivity_cloudradar) then
  ! add histogram bins to list of statistics that will be output.
  do n = 1,nbins_cloudradar-1
    binedges_cloudradar(n) = float(min_dBZbin_cloudradar + (n-1)*binwidth_cloudradar)
  end do
  call histogram_hbuf_init(count,microcount,'CLRAD', &
       'Area fraction with cloud radar', 'dBZ', &
       nbins_cloudradar, binedges_cloudradar, binname_cloudradar, &
       .true.) ! the last argument indicates that bins are named by dBZ value.

  ! profile of cloud fraction (area fraction w/ cloud radar dBZ>-40)
  name = 'CLDCRM40'
  longname = 'Area fraction with cloud radar dBZ > -40'
  units = ''
  call add_to_namelist(count,microcount,name,longname,units,0)

  ! profile of cloud fraction (area fraction w/ cloud radar dBZ>-30)
  name = 'CLDCRM30'
  longname = 'Area fraction with cloud radar dBZ > -30'
  units = ''
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if(masterproc) then
   write(*,*) 'Added ', microcount, ' arrays to statistics for M2005 microphysics'
end if

end subroutine micro_hbuf_init

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!! Note that only the fields declared in micro_hbuf_init() are allowed to
! be collected

subroutine micro_statistics()

use hbuffer, only: hbuf_put

real, dimension(nzm) :: tr0, tr2, frac30, frac40

real, dimension(nzm) :: zeros 

real tmp(2), factor_xy
integer i,j,k,m, n, ii, jj, nn, ncond, ioffset

call t_startf ('micro_statistics')

factor_xy = 1./float(nx*ny)

zeros(:) = 0.

do n = 1,nmicro_fields
   do k = 1,nzm
      tmp(1) = dz
      tmp(2) = dz/dtn
      tr0(k) = SUM(micro_field(1:nx,1:ny,k,n))
      tr2(k) = SUM(micro_field(1:nx,1:ny,k,n)*micro_field(1:nx,1:ny,k,n))

      mkwle(k,n) = mkwle(k,n)*tmp(2)!*lfac(n) ! resolved flux changed from W/m^2
      mkwsb(k,n) = mkwsb(k,n)*tmp(1)!*lfac(n) ! subgrid flux changed from W/m^2
      mksed(k,n) = mksed(k,n)!*lfac(n) ! sedimentation flux changed from W/m^2

   end do

!bloss(2019-07): !!!! OUTPUT NUMBER AS #/mg, NOT #/cm3 !!!!
!bloss   if(flag_number(n).eq.1) then
!bloss      ! remove factor of rho from number concentrations
!bloss      tr0(:) = tr0(:)*rho(:)
!bloss      tr2(:) = tr2(:)*rho(:)**2
!bloss      mkadv(1:nzm,n) = mkadv(1:nzm,n)*rho(:)
!bloss      mkdiff(1:nzm,n) = mkdiff(1:nzm,n)*rho(:)
!bloss      mtend(1:nzm,n) = mtend(1:nzm,n)*rho(:)
!bloss      stend(1:nzm,n) = stend(1:nzm,n)*rho(:)
!bloss      mklsadv(1:nzm,n) = mklsadv(1:nzm,n)*rho(:)
!bloss
!bloss   end if

!bloss/qt: output all microphysical fields
   if(dototalwater.OR.(n.ne.iqv)) then
     ! mean microphysical field
     call hbuf_put(trim(mkname(n)),tr0,mkoutputscale(n)*factor_xy)
   end if
  if(n.eq.iqv) then
      ! variance of microphysical field,  only for QTO (qv+qcl)
      call hbuf_put(trim(mkname(n))//'2',tr2,mkoutputscale(n)**2*factor_xy)
   end if

   if(flag_advect(n).eq.1) then
     ! do not rescale fluxes
     call hbuf_put(trim(mkname(n))//'FLXR',mkwle(1,n),factor_xy)
     call hbuf_put(trim(mkname(n))//'FLXS',mkwsb(1,n),factor_xy)
     call hbuf_put(trim(mkname(n))//'SDFLX',mksed(1,n),factor_xy)

     ! tendencies
     call hbuf_put(trim(mkname(n))//'ADV', &
          mkadv(:,n),mkoutputscale(n)*factor_xy*86400./dtn)
     call hbuf_put(trim(mkname(n))//'DIFF', &
          mkdiff(:,n),mkoutputscale(n)*factor_xy*86400./dtn)
     call hbuf_put(trim(mkname(n))//'LSADV', &
          mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
     call hbuf_put(trim(mkname(n))//'MPHY', &
          mtend(:,n),mkoutputscale(n)*factor_xy*86400.)
     call hbuf_put(trim(mkname(n))//'SED', &
          stend(:,n),mkoutputscale(n)*factor_xy*86400.)


     if(mod(nstep,nstat).eq.0) then
       mkstor(1:nzm,n) = (tr0(1:nzm) - mkstor(1:nzm,n))/dt/float(nstatis)
       call hbuf_put(trim(mkname(n))//'STRG', &
            mkstor(:,n),mkoutputscale(n)*factor_xy*86400.)
     else
       call hbuf_put(trim(mkname(n))//'STRG', &
            zeros(:),mkoutputscale(n)*factor_xy*86400.)
     end if

     if(mod(nstep,nstat).eq.0) then
       mkstor(1:nzm,n) = tr0(1:nzm)
     end if
   end if


   if((flag_wmass(n).gt.0).and.(n.ne.iqv)) then
      ! fractional area of microphysical field > 1.e-6
      call hbuf_put(trim(mkname(n))//'FRAC',mfrac(1,n),factor_xy)

      ! approx optical depth
      call hbuf_put('TAU'//trim(mkname(n)),trtau(:,n),factor_xy)

      !bloss (Apr 09):  This measure of effective radius is unreliable if the 
      !          microphysical fraction is not roughly uniform across
      !          the processors in a MPI run.  As a result, I am
      !          removing it from the outputs.  It is reliable if computed from
      !          the quantities TAU* and Q* in the output file.
!bloss      ! effective radius
!bloss      tr2(:) = 0.
!bloss      if(trtau(1,n).gt.0.) then
!bloss         tr2(1) = 1.e6*0.0018*rho(1)*dz*adz(1)*tr0(1)/trtau(1,n)
!bloss      end if

!bloss      do k = 2,nzm
!bloss         if(trtau(k,n).gt.trtau(k-1,n)) then
!bloss            tr2(k) = 1.e6*0.0018*rho(k)*dz*adz(k)*tr0(k)/(trtau(k,n)-trtau(k-1,n))
!bloss         end if
!bloss      end do
!bloss      call hbuf_put('EFFR'//trim(mkname(n)),tr2,1.)

      !bloss (Apr 09): Make an alternate statistic that can be used
      ! to easily compute the mean effective radius in a consistent
      ! way from optical depth.  This quantity Q*OEFFR is essentially
      ! the layer optical depth scaled in such a way that
      !
      !    EFFR = <Q*> / <Q*OEFFR>
      !
      ! where <.> is a time- and horizontal average.
      tr2(:) = 0.
      tr2(1) = trtau(1,n) / (1.e6*0.0018*rho(1)*dz*adz(1)*1.e-3)
      do k = 2,nzm
            tr2(k) = (trtau(k,n)-trtau(k-1,n)) / (1.e6*0.0018*rho(k)*dz*adz(k)*1.e-3) 
      end do
      call hbuf_put(trim(mkname(n))//'OEFFR',tr2,factor_xy)
      
   end if

   do ncond = 1,ncondavg
      do k = 1,nzm
         tr0(k) = SUM(micro_field(1:nx,1:ny,k,n)*condavg_mask(1:nx,1:ny,k,ncond))
      end do
!bloss (!!! #/mg NOT #/cm3 !!!)
!bloss      if(flag_number(n).eq.1) tr0(:) = tr0(:)*rho(:) ! remove factor of rho from number concentrations
      call hbuf_put(TRIM(mkname(n)) // TRIM(condavgname(ncond)), &
           tr0,mkoutputscale(n))
   end do

end do

if(dototalwater) then
  !bloss/qt: in total water formulation, fluxes of qv and qcl computed together.
  tr0(:) = mkwle(1:nzm,iqv) + mkwsb(1:nzm,iqv) ! qv + qcl tendencies
else
  tr0(:) = mkwle(1:nzm,iqv) + mkwsb(1:nzm,iqv) &
       + mkwle(1:nzm,iqcl) + mkwsb(1:nzm,iqcl)
end if

if(doicemicro) then
   tr0(:) = tr0(:) + mkwle(1:nzm,iqci) + mkwsb(1:nzm,iqci)
end if
call hbuf_put('QTFLUX',tr0,factor_xy*lcond) !bloss/diamond(2020-06): Add factor of Lv for output in W/m2

if(dototalwater) then
  !bloss/qt: add separate output for cloud liquid water
  !           and approx cloud liquid optical depth.
  do k = 1,nzm
    tr0(k) = SUM(cloudliq(1:nx,1:ny,k))
  end do
  call hbuf_put('QC',tr0,factor_xy*mkoutputscale(iqv))

  !bloss/qt: add separate conditional averages for cloud liquid water and vapor.
  do ncond = 1,ncondavg
    do k = 1,nzm
      tr0(k) = SUM(cloudliq(1:nx,1:ny,k)*condavg_mask(1:nx,1:ny,k,ncond))
    end do
    call hbuf_put('QC' // TRIM(condavgname(ncond)),tr0,mkoutputscale(iqv))
    do k = 1,nzm
      tr0(k) = SUM((micro_field(1:nx,1:ny,k,iqv)-cloudliq(1:nx,1:ny,k))*condavg_mask(1:nx,1:ny,k,ncond))
    end do
    call hbuf_put('QV' // TRIM(condavgname(ncond)),tr0,mkoutputscale(iqv))
  end do

else
  ! since vapor and cloud liquid mass are separate prognostic variables,
  !   we can report the latent heating tendency.  This can not be done
  !   in the total water formulation because one cannot distinguish between
  !   cloud liquid water tendencies due to advection and those due to phase changes.
  do k = 1,nzm
    tr0(k) = SUM(tmtend3d(1:nx,1:ny,k))
  end do
  call hbuf_put('QLAT',tr0,factor_xy*86400.)

  do ncond = 1,ncondavg
    do k = 1,nzm
      tr0(k) = SUM(tmtend3d(1:nx,1:ny,k)*condavg_mask(1:nx,1:ny,k,ncond))
    end do
    call hbuf_put('QLAT' // TRIM(condavgname(ncond)),tr0,86400.)
  end do
end if

call hbuf_put('TAUQC',trtau(:,iqv),factor_xy)

!bloss (Apr 09): Make an alternate statistic that can be used
! to easily compute the mean effective radius in a consistent
! way from optical depth.  This quantity Q*OEFFR is essentially
! the layer optical depth scaled in such a way that
!
!    EFFR = <Q*> / <Q*OEFFR>
!
! where <.> is a time- and horizontal average.
tr2(:) = 0.
tr2(1) = trtau(1,iqv) / (1.e6*0.0018*rho(1)*dz*adz(1)*1.e-3)
do k = 2,nzm
  tr2(k) = (trtau(k,iqv)-trtau(k-1,iqv)) / (1.e6*0.0018*rho(k)*dz*adz(k)*1.e-3) 
end do
call hbuf_put('QCOEFFR',tr2,factor_xy)

if(dopredictNc) then
  nn = 0.
  tmp(1)=0.
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      if(qcl(i,j,k).gt.0.) then
         tmp(1) = tmp(1) + micro_field(i,j,k,incl)*1.e-6
         nn = nn + 1
       end if
    end do
   end do      
  end do
  if (nn.gt.0) ncmn = ncmn + tmp(1)/dble(nn)
  ! RP - why not
  ! 1.e6*sum(micro_field(1:nx,1:ny,1:nzm,incl), mask = micro_field(1:nx,1:ny,1:nzm,incl) > 0.) / & 
  !      count(micro_field(1:nx,1:ny,1:nzm,incl) > 0.)
else
  ncmn = Nc0
end if
if(doprecip) then
  nn = 0.
  tmp(1)=0.
  do k=1,nzm
   do j=1,ny
    do i=1,nx 
      if(micro_field(i,j,k,iqr).gt.0.) then 
         tmp(1) = tmp(1) + micro_field(i,j,k,inr)*1.e-6
         nn = nn + 1
       end if
    end do
   end do
  end do
  if (nn.gt.0) then
      nrainy = nrainy + 1
      nrmn = nrmn + tmp(1)/dble(nn)
  end if
else
  nrmn = 0.
end if

if(do_output_micro_process_rates) then
  ioffset = 0

  do n = 1,nmicro_process_rates_warm_mass
    tr0(1:nzm) = micro_proc_rates(1:nzm,ioffset+n)
    call hbuf_put(TRIM(micro_process_rate_names_warm_mass(n)),tr0, &
         factor_xy*86400.*1.e3) ! kg/kg/s --> g/kg/day
  end do
  ioffset = ioffset + nmicro_process_rates_warm_mass

  do n = 1,nmicro_process_rates_warm_number
    tr0(1:nzm) = micro_proc_rates(1:nzm,ioffset+n)
    call hbuf_put(TRIM(micro_process_rate_names_warm_number(n)),tr0, &
         factor_xy*86400.*1.e-6) ! #/kg/s --> #/mg/day
  end do
  ioffset = ioffset + nmicro_process_rates_warm_number


  if(doprogaerosol) then
    do n = 1,nmicro_process_rates_progaer_mass
      tr0(1:nzm) = micro_proc_rates(1:nzm,ioffset+n)
      call hbuf_put(TRIM(micro_process_rate_names_progaer_mass(n)),tr0, &
           factor_xy*86400.*1.e3) ! kg/kg/s --> g/kg/day
    end do
    ioffset = ioffset + nmicro_process_rates_progaer_mass

    do n = 1,nmicro_process_rates_progaer_number
      tr0(1:nzm) = micro_proc_rates(1:nzm,ioffset+n)
      call hbuf_put(TRIM(micro_process_rate_names_progaer_number(n)),tr0, &
           factor_xy*86400.*1.e-6) ! #/kg/s --> #/mg/day
    end do
    ioffset = ioffset + nmicro_process_rates_progaer_number
  end if

  if(doicemicro) then
    do n = 1,nmicro_process_rates_cold_mass
      tr0(1:nzm) = micro_proc_rates(1:nzm,ioffset+n)
      call hbuf_put(TRIM(micro_process_rate_names_cold_mass(n)),tr0, &
           factor_xy*86400.*1.e3) ! kg/kg/s --> g/kg/day
    end do
    ioffset = ioffset + nmicro_process_rates_cold_mass

    do n = 1,nmicro_process_rates_cold_number
      tr0(1:nzm) = micro_proc_rates(1:nzm,ioffset+n)
      call hbuf_put(TRIM(micro_process_rate_names_cold_number(n)),tr0, &
           factor_xy*86400.*1.e-6) ! #/kg/s --> #/mg/day
    end do
    ioffset = ioffset + nmicro_process_rates_cold_number
  end if

end if

if (doprogaerosol.AND.do_m2011_scavenge) then
  ! output scavenging tendencies
  tr0(:) = scvtndqadclstat(:)
  call hbuf_put('SCVTQAW',tr0,factor_xy*86400.*1.e3) ! kg/kg/s --> g/kg/day

  tr0(:) = scvtndqadrstat(:)
  call hbuf_put('SCVTQAR',tr0,factor_xy*86400.*1.e3) ! kg/kg/s --> g/kg/day

  tr0(:) = scvtndnadclstat(:)
  call hbuf_put('SCVTNADC',tr0,factor_xy*86400.*1.e-6) ! #/kg/s --> #/mg/day

  tr0(:) = scvtndnadrstat(:)
  call hbuf_put('SCVTNADR',tr0,factor_xy*86400.*1.e-6) ! #/kg/s --> #/mg/day
end if

if(doreflectivity_cloudradar) then

  hist_cloudradar(:,:) = 0.
  frac30(:) = 0.
  frac40(:) = 0.

  if(dostatis_quickbeam) then

    !bloss: compute histogram of area fractions for cloud radar
    call compute_histogram(dBZ_cloudradar,nbins_cloudradar, binedges_cloudradar, &
         hist_cloudradar, 1, nx, 1, ny, 1, nzm, 1, nx, 1, ny, 1, nzm)

    do k = 1,nzm
      do j = 1,ny
        do i = 1,nx
          ! area fraction where cloud radar dBZ > -30, -40
          frac30(k) = frac30(k) + MAX(0.,SIGN(1.,dbZ_cloudradar(i,j,k)+30.))
          frac40(k) = frac40(k) + MAX(0.,SIGN(1.,dbZ_cloudradar(i,j,k)+40.))
        end do
      end do
    end do
    
  end if

  do n = 1,nbins_cloudradar
    call hbuf_put(TRIM(binname_cloudradar(n)),hist_cloudradar(1:nzm,n),factor_quickbeam*factor_xy)
  end do
  call hbuf_put('CLDCRM30',frac30,factor_quickbeam*factor_xy)
  call hbuf_put('CLDCRM40',frac40,factor_quickbeam*factor_xy)

end if ! if(doreflectivity_cloudradar)

call t_stopf ('micro_statistics')

end subroutine micro_statistics

!-----------------------------------------
subroutine micro_stat_2Dinit(ResetStorage)
  implicit none
  integer, intent(in) :: ResetStorage

  ! initialize microphysical 2D outputs as necessary

  if(ResetStorage.eq.1) then
    !bloss: If computing storage terms for individual hydrometeors,
    !         store initial profiles for computation of storage terms in budgets
    mkstor(1:nzm,:) = real(nx*ny)*mk0(1:nzm,:)
  end if

end subroutine micro_stat_2Dinit
!-----------------------------------------
subroutine micro_write_fields2D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1

  integer :: i, j, k, nn
  real(4) tmp(nx,ny,nzm)
  real :: coef
  character *80 long_name
  character *8 name
  character *10 units

  if(save2Davg) then
    write(*,*) ' No averaging of aerosol 2D outputs in MICRO_M2005_PA.  FIX THIS LATER IF NEEDED...'
    !    coef = 1./float(nsave2D)
    coef = 1.
  else
    coef = 1.
  end if

  if(doprogaerosol) then
    ! if isotopes are enabled, output surface fluxes (in water mass, not energy, units)
    !   of both standard isotope and heavy isotopologues.
    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do j=1,ny
      do i=1,nx
        tmp(i,j,1)=fluxbmk(i,j,inacc)*rhow(1)
      end do
    end do
    name='NAFLX'
    long_name='Instantaneous surface flux of aerosol number'
    units='/m2/s'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do j=1,ny
      do i=1,nx
        tmp(i,j,1)=aer_nflux_avg_xy(i,j)*rhow(1)/float(nsave2D)
      end do
    end do
    name='NAFLXAVG'
    long_name='Time-averaged surface flux of aerosol number'
    units='/m2/s'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)
    aer_nflux_avg_xy(:,:) = 0.

    ! if isotopes are enabled, output surface fluxes (in water mass, not energy, units)
    !   of both standard isotope and heavy isotopologues.
    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do j=1,ny
      do i=1,nx
        tmp(i,j,1)=fluxbmk(i,j,iqacc)*rhow(1)
      end do
    end do
    name='QAFLX'
    long_name='Instantaneous surface flux of aerosol mass'
    units='kg/m2/s'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do j=1,ny
      do i=1,nx
        tmp(i,j,1)=aer_qflux_avg_xy(i,j)*rhow(1)/float(nsave2D)
      end do
    end do
    name='QAFLXAVG'
    long_name='Time-averaged surface flux of aerosol mass'
    units='kg/m2/s'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)
    aer_qflux_avg_xy(:,:) = 0.

    ! if prognostic aerosols are enabled, output column-integrated quantities
    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*micro_field(i,j,k,inad)
        end do
      end do
    end do
    name='NAdPATH'
    long_name='Column dry aerosol (vertically integrated)'
    units='/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*micro_field(i,j,k,incl)
        end do
      end do
    end do
    name='NCPATH'
    long_name='Column CDNC (vertically integrated)'
    units='/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*qcl(i,j,k)
        end do
      end do
    end do
    name='QCPATH'
    long_name='Column cloud liquid mass (vertically integrated)'
    units='kg/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*micro_field(i,j,k,inr)
        end do
      end do
    end do
    name='NRPATH'
    long_name='Column rain number (vertically integrated)'
    units='/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*micro_field(i,j,k,iqr)
        end do
      end do
    end do
    name='QRPATH'
    long_name='Column rain mass (vertically integrated)'
    units='/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*qcl(i,j,k)*(micro_field(i,j,k,inacc)+micro_field(i,j,k,inr))
        end do
      end do
    end do
    name='NAQCPATH'
    long_name='Column integral of QC*NA (NAQCPATH/QCPATH=LWC-weighted total aerosol #, NAd+NC+NR)'
    units='kg/kg/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*qcl(i,j,k)*micro_field(i,j,k,incl)
        end do
      end do
    end do
    name='NCQCPATH'
    long_name='Column integral of QC*NC (useful for computing mass-weighted NC)'
    units='kg/kg/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    tmp(:,:,:) = 0.
    do k = 1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,1)=tmp(i,j,1)+rho(k)*dz*adz(k)*micro_field(i,j,k,inr)*micro_field(i,j,k,iqr)
        end do
      end do
    end do
    name='NRQRPATH'
    long_name='Column integral of QR*NR (useful for computing mass-weighted NR)'
    units='kg/kg/m2'
    call compress3D(tmp,nx,ny,1,name,long_name,units, &
         save2Dbin,dompi,rank,nsubdomains)


  end if


end subroutine micro_write_fields2D

!-----------------------------------------
subroutine micro_write_fields3D(nfields1)
  implicit none
  integer, intent(inout) :: nfields1
  character *80 long_name
  character *8 name
  character *10 units
  integer :: i, j, k, n, tens, ones
  real(4), dimension(nx,ny,nzm) :: tmp

  if(doreflectivity_cloudradar) then

    nfields1=nfields1+1
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tmp(i,j,k)=dBZ_cloudradar(i,j,k)
        end do
      end do
    end do
    name='dBZCLRAD'
    tens = floor(freq_cloudradar/10.)
    ones = floor(freq_cloudradar) - 10*tens
    long_name= char(48+tens) // char(48+ones) // &
         'GHz Cloud Radar Reflectivity'
    units='dBZ'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
         save3Dbin,dompi,rank,nsubdomains)

  end if

  if(do_m2011_scavenge.AND.doprogaerosol.AND.do_scav_3d_output) then
    !bloss(2020-10): Add 3d outputs of scavenging tendencies
    do n = 1,4
      nfields1=nfields1+1
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            tmp(i,j,k)=scav3d(i,j,k,n)
          end do
        end do
      end do
      name=TRIM(scavname(n))
      long_name= TRIM(scavlongname(n))
      units=TRIM(scavunits(n))
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains)
    end do
  end if

end subroutine micro_write_fields3D

!-----------------------------------------
subroutine satadj_liquid(nzm,tabs,qt,qc,pres)
  !bloss/qt: Utility routine based on cloud.f90 in 
  !  MICRO_SAM1MOM that was written by Marat Khairoutdinov.
  !  This routine performs a saturation adjustment for
  !  cloud liquid water only using a Newton method.
  !  While 20 iterations are allowed, most often this
  !  routine should exit in five iterations or less.
  !  Only a single calculation of the saturation vapor
  !  pressure is required in subsaturated air.

  implicit none

  integer, intent(in) :: nzm
  real, intent(inout), dimension(nzm) :: tabs ! absolute temperature, K
  real, intent(inout), dimension(nzm) :: qt  ! on input: qt; on output: qv
  real, intent(out), dimension(nzm) :: qc ! cloud liquid water, kg/kg
  real, intent(in), dimension(nzm) :: pres ! pressure, Pa

  real tabs1, dtabs, thresh, esat1, qsat1, fff, dfff
  integer k, niter

  integer, parameter :: maxiter = 20

  !bloss/qt: quick saturation adjustment to compute cloud liquid water content.
  do k = 1,nzm
    tabs1 = tabs(k) 
    esat1 = MIN(0.99*pres(k), polysvp(tabs1,0)) ! fix from Hugh for low pressure
    qsat1 = 0.622*esat1/ (pres(k) - esat1)
    qc(k) = 0. ! no cloud unless qt > qsat
    
    if (qt(k).gt.qsat1) then

      ! if unsaturated, nothing to do (i.e., qv=qt, T=Tl) --> just exit.
      ! if saturated, do saturation adjustment 
      !    (modeled after Marat's cloud.f90).

      ! generate initial guess based on above calculation of qsat
      dtabs = + fac_cond*MAX(0.,qt(k) - qsat1) &
           / ( 1. + lcond**2*qsat1/(cp*rv*tabs1**2) )
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(0.01, 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        ! borrow fix for low pressure from Hugh in module_mp_graupel.f90
        esat1 = MIN(0.99*pres(k), polysvp(tabs1,0))
        qsat1 = 0.622*esat1/ (pres(k) - esat1) ! saturation mixing ratio

        fff = tabs(k) - tabs1 + fac_cond*MAX(0.,qt(k) - qsat1)
        dfff = 1. + lcond**2*qsat1/(cp*rv*tabs1**2)
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      qc(k) = MAX( 0.,tabs1 - tabs(k) )/fac_cond ! cloud liquid mass mixing ratio
      qt(k) = qt(k) - qc(k) ! This now holds the water vapor mass mixing ratio.
      tabs(k) = tabs1 ! update temperature.
      
      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,nzm

end subroutine satadj_liquid

!========================================================================
subroutine cloudradar_init()
  ! initialize stuff for cloud radar simulator (QUICKBEAM)
  !   -- outputs include 3D snapshots of radar reflectivity as well
  !      as reflectivity histograms and profiles of area fractions 
  !      that exceed -40 and -30 dBZ.
  !   -- We call the QUICKBEAM to evaluate the reflectivity.  It was
  !      developed by John Haynes and Roger Marchand.
  implicit none
  integer :: ierr, n, i
  real :: delt, deltp

  ! inputs to radar_simulator_init
  integer, parameter :: nhclass_max = 5
  real, dimension(nhclass_max) :: hclass_dmin,hclass_dmax, &
       hclass_apm,hclass_bpm,hclass_rho, &
       hclass_p1,hclass_p2,hclass_p3
  integer,dimension(nhclass_max)  ::    hclass_type,hclass_phase
  logical     :: load_scale_LUTs_flag = .false.,update_scale_LUTs_flag = .false.
  character*240 :: LUT_file_name = './RUNDATA/quickbeam_lookuptable'

  ! check to see that nstatfrq is evenly divisible by nskip_quickbeam
  if(mod(nstatfrq,nskip_quickbeam).ne.0) then
    if(masterproc) then
      write(*,*) '*************************************************************'
      write(*,*) 'ERROR in MICRO_M2005_PA: '
      write(*,*) '        nskip_quickbeam must be a factor of nstatfrq'
      write(*,*) '*************************************************************'
    end if
    call task_abort()
  end if

  ! set up, output message and allocate variables
  nfields3D_micro = nfields3D_micro + 1
  nbins_cloudradar = 2 + FLOOR( &
       (max_dBZbin_cloudradar - min_dBZbin_cloudradar)/float(binwidth_cloudradar))
  if(masterproc) then
    write(*,*) '*************************************************************'
    write(*,*) 'MICRO_M2005_PA: Cloud radar reflectivity output enabled.'
999 format('             ', I4, ' bins from ', I4,' to ', I4, ' dBZ')
    write(*,999) nbins_cloudradar, min_dBZbin_cloudradar, max_dBZbin_cloudradar
    write(*,*) ' Uses quickbeam cloud radar simulator, v. 1.03'
    write(*,*) '   Copyright (c) 2006, J.M. Haynes.  All rights reserved.'
    write(*,*) '*************************************************************'
  end if
  ALLOCATE( &
       hist_cloudradar(nzm,nbins_cloudradar), &
       dBZ_cloudradar(nx,ny,nzm), &
       binedges_cloudradar(nbins_cloudradar-1), &
       binname_cloudradar(nbins_cloudradar), &
       STAT=ierr)
  if(ierr.ne.0) then
    write(*,*) 'Failed to allocate M2005 cloud radar reflectivity arrays on proc ', rank
    call task_abort()
  end if

  ! set up hydrometeor classes
  if(doicemicro) then
    nhclass = 5
  else
    nhclass = 2
  end if

  ! first, initialize arrays describing drop/particle size distributions
  !   to dummy values
  hclass_type(:) = -1
  hclass_phase(:) = -1
  hclass_dmin(:) = 0. ! unused for gamma/exponential distn.
  hclass_dmax(:) = 0.
  hclass_apm(:) = -1. ! optional way to specify ice density
  hclass_bpm(:) = -1.
  hclass_rho(:) = -1.
  hclass_p1(:) = -1. ! drop size distn parameters.
  hclass_p2(:) = -1.
  hclass_p3(:) = -1.

  ! Now, fill those arrays with meaningful values as necessary.
  !   put liquid in first two positions so that we can make the inputs smaller for 
  !   ice-free simulations

  ! cloud liquid water -- gamma distribution
  n = 1
  hclass_type(n) = 1
  hclass_phase(n) = 0
  hclass_rho(n) = 1000.
  hclass_p2(n) = 10. ! mean diameter.
  if(dofix_pgam) then
    hclass_p3(n) = pgam_fixed
  else
    ! PGAM formula in M2005 
    hclass_p3(n) = MAX(2., MIN(10., 1./(0.0005714*Nc0 + 0.2714)**2 - 1.) ) 
  end if
  ! subtract one by quickbeam convention
  hclass_p3(n) = hclass_p3(n) - 1.

  ! rain -- exponential distribution
  n = 2
  hclass_type(n) = 2
  hclass_phase(n) = 0
  hclass_rho(n) = 997.
  ! note: since we are supplying effective radii, this fixed value for lambda
  !         should not enter the reflectivity computations.
  hclass_p2(n) = 0.01 ! use geometric mid point of MAX/MIN range

  if(doicemicro) then

    ! cloud ice -- exponential distribution
    n = 3
    hclass_type(n) = 2
    hclass_phase(n) = 1
    hclass_rho(n) = 500.
    ! note: since we are supplying effective radii, this fixed value for lambda
    !         should not enter the reflectivity computations.
    hclass_p2(n) = 0.05 ! use geometric mid point of MAX/MIN range

    ! snow -- exponential distribution
    n = 4
    hclass_type(n) = 2
    hclass_phase(n) = 1
    hclass_rho(n) = 100.
    ! note: since we are supplying effective radii, this fixed value for lambda
    !         should not enter the reflectivity computations.
    hclass_p2(n) = 0.007 ! use geometric mid point of MAX/MIN range

    ! graupel/hail -- exponential distribution
    n = 5
    hclass_type(n) = 2
    hclass_phase(n) = 1
    if(dohail) then
      hclass_rho(n) = 900. ! hail
    else
      hclass_rho(n) = 400. ! graupel
    end if
    ! note: since we are supplying effective radii, this fixed value for lambda
    !         should not enter the reflectivity computations.
    hclass_p2(n) = 0.005 ! use geometric mid point of MAX/MIN range

  end if ! if(doicemicro)

  call radar_simulator_init( &
       freq_cloudradar, k2_cloudradar, & ! inputs
       usegasabs_cloudradar, doray_cloudradar, missing_value_cloudradar, &
       nhclass, & 
       hclass_type,hclass_phase, &
       hclass_dmin,hclass_dmax, &
       hclass_apm,hclass_bpm,hclass_rho, &
       hclass_p1,hclass_p2,hclass_p3, &
       load_scale_LUTs_flag,update_scale_LUTs_flag,LUT_file_name, &
       hp_cloudradar ) ! output

end subroutine cloudradar_init

!========================================================================
subroutine twodigit_integer_to_string(itmp,name,longname)
  implicit none
  integer, intent(in) :: itmp
  character(LEN=3), intent(out) :: name, longname

  integer :: tens, ones

  tens = floor(float(ABS(itmp))/10.)
  ones = ABS(itmp) - 10*tens

  if(itmp.lt.0) then
    name = 'M' // char(48+tens) // char(48+ones)  ! no negative sign
    ! in netcdf variable names
    longname = '-' // char(48+tens) // char(48+ones)  ! negative sign here
  else
    name = char(48+tens) // char(48+ones)
    longname = char(48+tens) // char(48+ones)
  end if
end subroutine twodigit_integer_to_string

!========================================================================
! auxilliary routine to set up a bunch of output fields for cloud radar histograms
subroutine histogram_hbuf_init(count, trcount, shortstring, longstring, &
     hist_unit, nbin, binedges, binnames, dobinedges_in_names)
  implicit none

  integer, intent(inout) :: count, trcount
  character(LEN=*), intent(in) :: shortstring
  character(LEN=*), intent(in) :: longstring
  character(LEN=*), intent(in) :: hist_unit ! units of thing being binned
  integer, intent(in) :: nbin
  real, dimension(nbin-1), intent(in) :: binedges
  CHARACTER(LEN=8), DIMENSION(nbin), intent(out) :: binnames
  logical, intent(in) :: dobinedges_in_names

  character*8 name
  character*80 longname
  character*10 units

  character*3 :: num_string, numlo_string
  character*3 :: num_longstring, numlo_longstring

  integer :: n

  ! long string
  call twodigit_integer_to_string(NINT(binedges(1)), num_string, num_longstring)
  longname = longstring // ' ' // hist_unit // ' < ' // num_longstring

  ! short string
  if(dobinedges_in_names) then
    name = shortstring // 'LO'
  else
    call twodigit_integer_to_string(1, num_string, num_longstring)
    name = shortstring // num_string
  end if

  units = ''
  call add_to_namelist(count,trcount,name,longname,units,0)
  binnames(1) = TRIM(name)

  do n = 2,nbin-1
    ! long string
    call twodigit_integer_to_string(NINT(binedges(n-1)), numlo_string, numlo_longstring)
    call twodigit_integer_to_string(NINT(binedges(n)), num_string, num_longstring)
    longname = longstring // ' ' // numlo_longstring // ' < ' // hist_unit &
         // ' < ' // num_longstring

    ! short string
    if(.NOT.dobinedges_in_names) then
      ! put bin number in name
      call twodigit_integer_to_string(n, num_string, num_longstring)
    end if
    name = shortstring // num_string
    units = ''
    call add_to_namelist(count,trcount,name,longname,units,0)
    binnames(n) = TRIM(name)
  end do

  ! long string 
  longname = longstring // ' ' // hist_unit // ' > ' // num_longstring

  if(dobinedges_in_names) then
    name = shortstring // 'HI'
  else
    call twodigit_integer_to_string(nbin, num_string, num_longstring)
    name = shortstring // num_string
  end if

  units = ''
  call add_to_namelist(count,trcount,name,longname,units,0)
  binnames(nbin) = TRIM(name)

end subroutine histogram_hbuf_init

! auxilliary routine for processing cloud radar output into bins
subroutine compute_histogram( field, nbin, binedges, hist, &
     imin, imax, jmin, jmax, kmin, kmax, i1, i2, j1, j2, k1, k2)

  implicit none

  integer, intent(in) :: imin, imax, jmin, jmax, kmin, kmax, &
       i1, i2, j1, j2, k1, k2, &
       nbin
  real, intent(in) :: field(imin:imax, jmin:jmax, kmin:kmax)
  real, intent(in) :: binedges(nbin-1)
  real, intent(out) :: hist(k1:k2,nbin)

  integer :: k, n
  real :: tmp(i1:i2, j1:j2)

  hist(k1:k2,1:nbin) = 0.

  do k = k1, k2

    tmp(:,:) = 0.
    WHERE (field(i1:i2,j1:j2,k).LE.binedges(1)) tmp(i1:i2,j1:j2) = 1.
    hist(k,1) = hist(k,1) + SUM(tmp)

    do n = 2,nbin-1
      tmp(:,:) = 0.
      WHERE (field(i1:i2,j1:j2,k).LE.binedges(n) &
           .AND. field(i1:i2,j1:j2,k).GT.binedges(n-1)) &
           tmp(i1:i2,j1:j2) = 1.
      hist(k,n) = hist(k,n) + SUM(tmp)
    end do

    tmp(:,:) = 0.
    WHERE (field(i1:i2,j1:j2,k).GT.binedges(nbin-1)) tmp(i1:i2,j1:j2) = 1.
    hist(k,nbin) = hist(k,nbin) + SUM(tmp)

  end do

end subroutine compute_histogram

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

logical function micro_provides_reffc()
  micro_provides_reffc = douse_reffc
end function micro_provides_reffc

logical function micro_provides_reffi()
  micro_provides_reffi = douse_reffi
end function micro_provides_reffi

function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffl() ! cloud water plus drizzle
  real, dimension(nx,ny,nzm) :: Get_reffl
  Get_reffl = 2.5
  where((qpl+qcl).GT.0.)
    Get_reffl = (qpl+qcl)/(qpl/reffr+qcl/reffc) 
  end where
end function Get_reffl

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi

function Get_nca() ! aerosol
  real, pointer, dimension(:,:,:) :: Get_nca
  Get_nca = 0.
end function Get_nca

end module microphysics



