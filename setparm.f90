
subroutine setparm

!       initialize parameters:

use vars
!use micro_params
use params
use microphysics, only: micro_setparm
use sgs, only: sgs_setparm
use movies, only : irecc
use instrument_diagnostics, only: zero_instr_diag
use grid, only: dompi
implicit none

integer icondavg, ierr, ios_uw, ios_kuang, ios_missing_namelist, place_holder, dowtg_num

real ttheta_tot

NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
                dolongwave, doshortwave, dosgs, dz, doconstdz, &
                docoriolis, docoriolisz, dosurface, dolargescale, doradforcing, &
								fluxt0,fluxq0,tau0,tabs_s,z0,nelapse, dt, dx, dy,  &
                fcor, ug, vg, nstop, caseid, case_restart,caseid_restart, &
								nstat, nstatfrq, nprint, nrestart, doradsimple, &
								nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
								donudging_uv, donudging_tq, &
                donudging_t, donudging_q, tauls,tautqls,&
                nudging_uv_z1, nudging_uv_z2, nudging_t_z1, nudging_t_z2, &
                nudging_q_z1, nudging_q_z2, dofplane, &
                dovariable_tauz,  variable_tauz_offset_above_inversion, & !bloss(2019-10)
                variable_tauz_thickness_of_onset, variable_tauz_minimum_height, & !bloss(2019-10)
                nudging_t_zramp, nudging_q_zramp, & !bloss(2019-10)
								timelargescale, longitude0, latitude0, day0, nrad, &
								CEM,LES,OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, soil_wetness, &
                doensemble, nensemble, dowallx, dowally, &
                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, &
                docolumn, save2Dbin, save2Davg, save3Dbin, &
                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
	        			doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
                domodis, domisr, dodynamicocean, ocean_type, delta_sst, &
                depth_slab_ocean, Szero, deltaS, timesimpleocean, &
								dosolarconstant, solar_constant, zenith_angle, rundatadir, &
                dotracers, output_sep, perturb_type, &
                doSAMconditionals, dosatupdnconditionals, &
                doscamiopdata, iopfile, dozero_out_day0, &
                use_scam_reference_sounding, set_obs_to_reference_sounding, &
                use_scam_reference_sounding_for_wtg, use_scam_initial_sounding, & !bloss(2018-11-30/2019-02-22)
                nstatmom, nstatmomstart, nstatmomend, savemomsep, savemombin, &
                nmovie, nmoviestart, nmovieend, nrestart_skip, &
                bubble_x0,bubble_y0,bubble_z0,bubble_radius_hor, &
                bubble_radius_ver,bubble_dtemp,bubble_dq, dosmoke, dossthomo, &
                rad3Dout, nxco2, dosimfilesout, notracegases, &
                ncycle_max, ncycle_min, ncycle0, cfl_safety_factor, &
                do_chunked_energy_budgets, do_chunked_momentum_budgets, &
                nchunk_x_gl, nchunk_y_gl, nsaveMSE

NAMELIST /UWOPTIONS/ rad_simple_fluxdiv1, &
     rad_simple_fluxdiv2, rad_simple_kappa, &
     dofixdivg, do_linear_subsidence, divg_ls, divg_lapse, &
     nelapsemin, doExtrapolate_UpperBound, &
     doFixedWindSpeedForSurfaceFluxes, WindSpeedForFluxes, &
     doDerbyshire, Derbyshire_z1, Derbyshire_z2, Derbyshire_z3, &
     Derbyshire_RelH_Low, Derbyshire_RelH_High, &
     Derbyshire_theta0, Derbyshire_LapseRate, &
     Derbyshire_tau, compute_advection_everywhere, &
     dowtg_blossey_etal_JAMES2009, &
     dowtg_qnudge, itau_wtg_qnudge, &
     dowtg_tnudge, itau_wtg_tnudge,  taulz_wtg_tnudge, &
     tauz0_wtg_qnudge, taulz_wtg_qnudge, &
     am_wtg, am_wtg_exp, am_tscale, lambda_wtg, nudge_to_sounding_winds, &
     doSmoothDamping, zbot_SmoothDamping, tau_SmoothDamping, &
     doenforce_cgils_qfloor, ztop_qfloor, qfloor, tau_qfloor

! Options added by Kuang Lab at Harvard
NAMELIST /KUANG_OPTIONS/ dompiensemble, doradtendency, troptend, &
            dowtg_raymondzeng_QJRMS2005, dowtg_daleuetal_JAMES2015, dowtg_decomp2022, dowtgLBL, &
            boundstatic, tau_wtg, tauscale_h, tauscale_f, tau_tscale, dthetadz_min, &
            dooceantimeperturb, tabs_ptscale, tabs_pamp, tabs_pphase

!bloss: Create dummy namelist, so that we can figure out error code
!       for a mising namelist.  This lets us differentiate between
!       missing namelists and those with an error within the namelist.
NAMELIST /BNCUIODSBJCB/ place_holder

!----------------------------------
!  Read namelist variables from the standard input:
!------------

open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')
read (55,PARAMETERS,IOSTAT=ierr)
if (ierr.ne.0) then
  ! try to get a more useful error message
  rewind(55)
  read (55,PARAMETERS)
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in PARAMETERS namelist'
        call task_abort()
end if
close(55)


!----------------------------------
!  Read namelist for uw options from same prm file:
!------------
open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

!bloss: get error code for missing namelist (by giving the name for
!       a namelist that doesn't exist in the prm file).
read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
rewind(55) !note that one must rewind before searching for new namelists

!bloss: read in UWOPTIONS namelist
read (UNIT=55,NML=UWOPTIONS,IOSTAT=ios_uw)
if (ios_uw.ne.0) then
  if(masterproc) write(*,*) 'ios_missing_namelist = ', ios_missing_namelist
  if(masterproc) write(*,*) 'ios for UWOPTIONS = ', ios_uw
   !namelist error checking
   if(ios_uw.ne.ios_missing_namelist) then
     rewind(55) !note that one must rewind before searching for new namelists
     read (UNIT=55,NML=UWOPTIONS)
     if(masterproc) then
       write(*,*) '****** ERROR: bad specification in UWOPTIONS namelist'
     end if
      call task_abort()
   elseif(masterproc) then
      write(*,*) '****************************************************'
      write(*,*) '******** No UWOPTIONS namelist in prm file *********'
      write(*,*) '****************************************************'
   end if
end if
close(55)

!----------------------------------
!  Read namelist for kuang_lab options from same prm file:
!------------
open(55,file='./'//trim(case)//'/prm', status='old',form='formatted')

!bloss: get error code for missing namelist (by giving the name for
!       a namelist that doesn't exist in the prm file).
read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
rewind(55) !note that one must rewind before searching for new namelists

! MPI Ensemble Run: read in KUANG_OPTIONS namelist (Nathanael Wong, 2022)
read (UNIT=55,NML=KUANG_OPTIONS,IOSTAT=ios_kuang)
if (ios_kuang.ne.0) then
  if(masterproc) write(*,*) 'ios_missing_namelist = ', ios_missing_namelist
  if(masterproc) write(*,*) 'ios for KUANG_OPTIONS = ', ios_kuang
   !namelist error checking
   if(ios_kuang.ne.ios_missing_namelist) then
     rewind(55) !note that one must rewind before searching for new namelists
     read (UNIT=55,NML=KUANG_OPTIONS)
     if(masterproc) then
       write(*,*) '****** ERROR: bad specification in KUANG_OPTIONS namelist'
     end if
      call task_abort()
   elseif(masterproc) then
      write(*,*) '****************************************************'
      write(*,*) '****** No KUANG_OPTIONS namelist in prm file *******'
      write(*,*) '****************************************************'
   end if
end if
close(55)

! write namelist values out to file for documentation
if(masterproc) then
      open(UNIT=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.nml',&
            form='formatted')
      write (55,NML=PARAMETERS)
      write (55,NML=UWOPTIONS)
      write (55,NML=KUANG_OPTIONS)
      write(55,*)
      close(55)
end if

!------------------------------------
!  Set parameters


        ! Allow only special cases for separate output:

        output_sep = output_sep.and.RUN3D
        if(output_sep)  save2Dsep = .true.

	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif

        if(docoriolis.and..not.dofplane.or.doradlat) dowally=.true.

	if(ny.eq.1) dy=dx
	dtn = dt

	notopened2D = .true.
	notopened3D = .true.

        call zero_instr_diag() ! initialize instruments output
        call sgs_setparm() ! read in SGS options from prm file.
        call micro_setparm() ! read in microphysical options from prm file.

        if(dosmoke) then
           epsv=0.
        else
           epsv=0.61
        endif

        if((nstatmomend.gt.nstatmomstart).AND.(nstop.ge.nstatmomstart)) then
          ! if moment statistics will be produced, check to make sure that
          !   navgmom_x and _y are set properly.
          if(navgmom_x.lt.0.or.navgmom_y.lt.0) then
            if(masterproc) then
              write(*,*) ' ******* ERROR in statmom specifications ******'
              write(*,*) 'Moment statistics will not be produced unless'
              write(*,*) 'navgmom_x and navgmom_y >= 0.  Model is stopping...'
            end if
            call task_abort()
          end if
        end if

        if(doseawater) then
          salt_factor = 0.981
        else
          salt_factor = 1.
        end if


        if(tautqls.eq.99999999.) tautqls = tauls

        !===============================================================
        ! KUANG_LAB ADDITION

        if(dompiensemble.AND.dompi) then
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Using the Kuang_Lab Ensemble Run Method'
            write(*,*) '  This will turn off MPI in the model run, such that'
            write(*,*) '  each subdomain is run independently of each other.'
            write(*,*) '  However, MPI is turned on for saving of output and'
            write(*,*) '  restart files.'
            write(*,*) '*********************************************************'
          end if
        else if(dompiensemble.AND.(.NOT.dompi)) then
          dompiensemble = .false.
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Do not use the Kuang_Lab Ensemble Run Method'
            write(*,*) '  MPI is not called because number of processors = 1'
            write(*,*) '  Setting dompiensemble to FALSE'
            write(*,*) '*********************************************************'
          end if
        end if

        !===============================================================
        ! Weak Temperature Gradient Approximation Schemes

        dowtg_num = 0

        if(dowtg_blossey_etal_JAMES2009) then
          if(masterproc) write(*,*) 'WTG (based on BBW09 in JAMES) is being used'
          am_wtg = am_wtg/86400. ! convert from 1/d to 1/s.
          dowtg_num = dowtg_num + 1
        end if

        if(dowtg_raymondzeng_QJRMS2005) then
          if(masterproc) write(*,*) 'WTG (based on Raymond and Zeng [2005]) is being used'
          tau_wtg = tau_wtg * 3600. ! convert from units of hours to units of sec.
          tau_wtg = 1 / tau_wtg      ! convert from sec to sec^-1
          dowtg_num = dowtg_num + 1
        end if

        if(dowtg_daleuetal_JAMES2015) then
          if(masterproc) write(*,*) 'WTG (based on Daleu et al. [2015]) is being used'
          tau_wtg = tau_wtg * 3600. ! convert from units of hours to units of sec.
          tau_wtg = 1 / tau_wtg      ! convert from sec to sec^-1
          dowtg_num = dowtg_num + 1
        end if

        if(dowtg_decomp2022) then
          if(masterproc) write(*,*) 'WTG (Spectral Decomposition into half- and full-sine) is being used'
          tau_wtg = tau_wtg * 3600. ! convert from units of hours to units of sec.
          tau_wtg = 1 / tau_wtg      ! convert from sec to sec^-1

          if(tauscale_h.gt.1) tauscale_h = 1
          if(tauscale_f.gt.1) tauscale_f = 1
          if(tauscale_h.lt.0) tauscale_h = 0
          if(tauscale_f.lt.0) tauscale_f = 0

          dowtg_num = dowtg_num + 1
        end if

        if(dowtg_num.GT.1) then
          if(masterproc) then
            write(*,*) '********************************************************'
            write(*,*) '  More than one of the available WTG schemes has been'
            write(*,*) '  called.  Please select only one of these schemes to use.'
            write(*,*) '********************************************************'
          end if
          call task_abort()
        end if

        !===============================================================
        ! UW ADDITION

        if((.NOT.save2Dbin).OR.(.NOT.save3Dbin)) then
          if(masterproc) then
            write(*,*) '*********************************************************'
            write(*,*) '  Compressed integer 2D and 3D output is disabled in '
            write(*,*) '  the UW version of SAM.  Please set save2Dbin=.true. '
            write(*,*) '  and save3Dbin=.true. in the PARAMETERS namelist. '
            write(*,*) '  If you really want compressed integer output, talk'
            write(*,*) '  to Peter.'
            write(*,*) '*********************************************************'
          end if
          call task_barrier()
          call task_abort()
        end if



        !bloss: set up conditional averages
        ncondavg = 1 ! always output CLD conditional average
        if(doSAMconditionals) ncondavg = ncondavg + 2
        if(dosatupdnconditionals) ncondavg = ncondavg + 3
        if(allocated(condavg_factor)) then ! avoid double allocation when nrestart=2
          DEALLOCATE(condavg_factor,condavg_mask,condavgname,condavglongname)
        end if
        ALLOCATE(condavg_factor(nzm,ncondavg), & ! replaces old cloud_factor, core_factor
             condavg_mask(nx,ny,nzm,ncondavg), & ! nx x ny x nzm indicator arrays
             condavgname(ncondavg), & ! short names (e.g. CLD, COR, SATUP)
             condavglongname(ncondavg), & ! long names (e.g. cloud, core, saturated updraft)
             STAT=ierr)
        if(ierr.ne.0) then
             write(*,*) '**************************************************************************'
             write(*,*) 'ERROR: Could not allocate arrays for conditional statistics in setparm.f90'
             call task_abort()
        end if

        ! indicators that can be used to tell whether a particular average
        !   is present.  If >0, these give the index into the condavg arrays
        !   where this particular conditional average appears.
        icondavg_cld = -1
        icondavg_cor = -1
        icondavg_cordn = -1
        icondavg_satup = -1
        icondavg_satdn= -1
        icondavg_env = -1

        icondavg = 0
        icondavg = icondavg + 1
        condavgname(icondavg) = 'CLD'
        condavglongname(icondavg) = 'cloud'
        icondavg_cld = icondavg

        if(doSAMconditionals) then
           icondavg = icondavg + 1
           condavgname(icondavg) = 'COR'
           condavglongname(icondavg) = 'core'
           icondavg_cor = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'CDN'
           condavglongname(icondavg) = 'downdraft core'
           icondavg_cordn = icondavg
        end if

        if(dosatupdnconditionals) then
           icondavg = icondavg + 1
           condavgname(icondavg) = 'SUP'
           condavglongname(icondavg) = 'saturated updrafts'
           icondavg_satup = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'SDN'
           condavglongname(icondavg) = 'saturated downdrafts'
           icondavg_satdn = icondavg

           icondavg = icondavg + 1
           condavgname(icondavg) = 'ENV'
           condavglongname(icondavg) = 'unsaturated environment'
           icondavg_env = icondavg
        end if

        ! END UW ADDITIONS
        !===============================================================

        irecc = 1


end
