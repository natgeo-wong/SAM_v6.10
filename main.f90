program crm

!       Main module.

use vars
use hbuffer
use params, only: dompiensemble
use microphysics
use sgs
use tracers
use movies, only: init_movies
use mse, only: init_MSE_tendency, compute_and_increment_MSE_tendency, &
     init_momentum_tendency, compute_and_increment_momentum_tendency, &
     initializeMSE
use simple_ocean, only: sst_perturb
implicit none

integer k, icyc, nn, nstatsteps
double precision cputime, oldtime, init_time !bloss wallclocktime
double precision usrtime, systime
integer itmp1(1), oldstep
logical :: log_Exists
!-------------------------------------------------------------------
! determine the rank of the current task and of the neighbour's ranks

call task_init() 
!-------------------------------------------------------------------
!bloss: Make sure that the run will not restart if it does not stop cleanly
if(masterproc) then
  inquire(file='ReadyForRestart',exist=log_Exists)
  if(.NOT.log_Exists) then
    open(unit=47,file='ReadyForRestart',form='FORMATTED',status='NEW')
    close(47)
  end if

  open(unit=47,file='ReadyForRestart',form='FORMATTED',status='REPLACE')
  rewind(47)
  write(47,991) 
991 format('FALSE')
  close(47)
end if

!------------------------------------------------------------------
! print time, version, etc

if(masterproc) call header()	
!------------------------------------------------------------------
! Initialize timing library.  2nd arg 0 means disable, 1 means enable

   call t_setoptionf (1, 0)
   call t_initializef ()

   call t_startf ('total')
   call t_startf ('initialize')
!------------------------------------------------------------------
! Get initial time of job

call t_stampf(init_time,usrtime,systime)
!------------------------------------------------------------------

call init()     ! initialize some statistics arrays
call setparm()	! set all parameters and constants

!------------------------------------------------------------------
! Initialize or restart from the save-dataset:

if(nrestart.eq.0) then
   day=day0 
   call setgrid() ! initialize vertical grid structure
   call setdata() ! initialize all variables
elseif(nrestart.eq.1) then
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call sgs_init()
   call micro_init()  !initialize microphysics
   if(do_chunked_energy_budgets) call initializeMSE()
elseif(nrestart.eq.2) then  ! branch run
   call read_all()
   call setgrid() ! initialize vertical grid structure
   call diagnose()
   call setparm() ! overwrite the parameters
   call sgs_init()
   call micro_init()  !initialize microphysics
   if(do_chunked_energy_budgets) call initializeMSE()
   nstep = 0
   day0 = day
else
   print *,'Error: confused by value of NRESTART'
   call task_abort() 
endif

call init_movies()
call stat_2Dinit(1) ! argument of 1 means storage terms in stats are reset
call tracers_init() ! initialize tracers
call setforcing()
call printout()
call write_fields3D() ! output initial condition
!------------------------------------------------------------------
!  Initialize statistics buffer:

call hbuf_init()
	
!------------------------------------------------------------------
nstatis = nstat/nstatfrq
nstat = nstatis * nstatfrq
nstatsteps = 0
call t_stopf ('initialize')
!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------
call t_stampf(oldtime,usrtime,systime)
oldstep = nstep


do while(nstep.lt.nstop.and.nelapse.gt.0) 
        
  nstep = nstep + 1
  time = time + dt
  day = day0 + nstep*dt/86400.
  nelapse = nelapse - 1
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  ! MPI Ensemble run: turn off mpi entering each loop (Song Qiyu, 2022)
  if(dompiensemble) dompi = .false.

  call kurant()

  total_water_before = total_water()
  total_water_evap = 0.
  total_water_prec = 0.
  total_water_ls = 0.

  do icyc=1,ncycle

     icycle = icyc
     dtn = dt/ncycle
     dt3(na) = dtn
     dtfactor = dtn/dt

     if(mod(nstep,nstatis).eq.0.and.icycle.eq.ncycle) then
        nstatsteps = nstatsteps + 1
        dostatis = .true.
        if(masterproc) print *,'Collecting statistics...'
     else
        dostatis = .false.
     endif

     !bloss:make special statistics flag for radiation,since it's only updated at icycle==1.
     dostatisrad = .false.
     if(mod(nstep,nstatis).eq.0.and.icycle.eq.1) dostatisrad = .true.

!---------------------------------------------
!    Perturb ocean SST in time sinusoidally around mean SST

     if(.not.dosfcforcing.and.dooceantimeperturb) call sst_perturb()

!---------------------------------------------
!  	the Adams-Bashforth scheme in time

     call abcoefs()
 
!---------------------------------------------
!  	initialize stuff: 
	
     call zero()

     if(do_chunked_momentum_budgets) call init_momentum_tendency() !bloss: see mse.f90
!-----------------------------------------------------------
!       Buoyancy term:
	     
     call buoyancy()

     if(do_chunked_momentum_budgets) call compute_and_increment_momentum_tendency('buoy') !bloss: see mse.f90
!------------------------------------------------------------

     total_water_ls =  total_water_ls - total_water()

!------------------------------------------------------------
!  for mse.f90, store s_li, h_f and q_tot before forcing.
     if(do_chunked_energy_budgets) call init_MSE_tendency()
     if(do_chunked_momentum_budgets) call init_momentum_tendency() !bloss: see mse.f90

!------------------------------------------------------------
!       Large-scale and surface forcing:
     call forcing()

!------------------------------------------------------------
!  get increment from forcing and add it to chunk averages.
     if(do_chunked_energy_budgets) call compute_and_increment_MSE_tendency('lsf')
     if(do_chunked_momentum_budgets) call compute_and_increment_momentum_tendency('lsf') !bloss: see mse.f90

!------------------------------------------------------------
!  for mse.f90, store s_li, h_f and q_tot before nudging/damping/upperbound
     if(do_chunked_energy_budgets) call init_MSE_tendency()
     if(do_chunked_momentum_budgets) call init_momentum_tendency() !bloss: see mse.f90

!----------------------------------------------------------
!       Nadging to sounding:

     call nudging()

!----------------------------------------------------------
!   	spange-layer damping near the upper boundary:

     if(dodamping) call damping()

!------------------------------------------------------------
!  get increment from nudging/damping and add it to chunk averages.
     if(do_chunked_energy_budgets) call compute_and_increment_MSE_tendency('misc')
     if(do_chunked_momentum_budgets) call compute_and_increment_momentum_tendency('misc') !bloss: see mse.f90

!----------------------------------------------------------

     total_water_ls =  total_water_ls + total_water()

!---------------------------------------------------------
!   Ice fall-out
   
      if(docloud) then
          call ice_fall()
      end if

!----------------------------------------------------------
!     Update scalar boundaries after large-scale processes:

     call boundaries(3)

!---------------------------------------------------------
!     Update boundaries for velocities:

      call boundaries(0)

!-----------------------------------------------
!     surface fluxes:

     if(dosurface) call surface()
     call micro_flux() !bloss: moved here for consistency of isotopic fluxes

!-----------------------------------------------------------
!  SGS physics:

     if (dosgs) call sgs_proc()

!----------------------------------------------------------
!     Fill boundaries for SGS diagnostic fields:

     call boundaries(4)

     if(do_chunked_momentum_budgets) call init_momentum_tendency() !bloss: see mse.f90
!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!----------------------------------------------------------
!	SGS effects on momentum:

     if(dosgs) call sgs_mom()

     if(do_chunked_momentum_budgets) call compute_and_increment_momentum_tendency('eddy') !bloss: see mse.f90
     if(do_chunked_momentum_budgets) call init_momentum_tendency() !bloss: see mse.f90
!-----------------------------------------------------------
!       Coriolis force:
	     
     if(docoriolis) call coriolis()
	 
     if(do_chunked_momentum_budgets) call compute_and_increment_momentum_tendency('lsf') !bloss: see mse.f90
     if(do_chunked_momentum_budgets) call init_momentum_tendency() !bloss: see mse.f90
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

     if(do_chunked_momentum_budgets) call compute_and_increment_momentum_tendency('pgrad') !bloss: see mse.f90
!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
!  Note that at the end of the call, the velocities are in nondimensional form.
	 
     call adams()

!----------------------------------------------------------
!     Update boundaries for all prognostic scalar fields for advection:

     call boundaries(2)

!---------------------------------------------------------
!      advection of scalars :

     call advect_all_scalars()

!-----------------------------------------------------------
!    Convert velocity back from nondimensional form:

      call uvw()

!----------------------------------------------------------
!     Update boundaries for scalars to prepare for SGS effects:

     call boundaries(3)
   
!------------------------------------------------------------
!  for mse.f90, store s_li, h_f and q_tot before SGS transport
     if(do_chunked_energy_budgets) call init_MSE_tendency()

!---------------------------------------------------------
!      SGS effects on scalars :

     if (dosgs) call sgs_scalars()

!------------------------------------------------------------
!  get increment from advection/diffusion and add it to chunk averages.
     if(do_chunked_energy_budgets) call compute_and_increment_MSE_tendency('eddy')

!------------------------------------------------------------
!  for mse.f90, store s_li, h_f and q_tot before upperbound
     if(do_chunked_energy_budgets) call init_MSE_tendency()

!-----------------------------------------------------------
!       Handle upper boundary for scalars

     total_water_ls =  total_water_ls - total_water() !bloss: Include any water changes in large-scale diagnostic
     if(doupperbound) call upperbound()
     total_water_ls =  total_water_ls + total_water()

!------------------------------------------------------------
!  get increment from upperbound and add it to chunk averages.
     if(do_chunked_energy_budgets) call compute_and_increment_MSE_tendency('misc')

!-----------------------------------------------------------
!       Cloud condensation/evaporation and precipitation processes:

     if(do_chunked_energy_budgets) call init_MSE_tendency() ! for mse.f90, store fields before micro_proc()

      if(docloud.or.dosmoke) call micro_proc()

     if(do_chunked_energy_budgets) call compute_and_increment_MSE_tendency('mphy')

!----------------------------------------------------------
!  Tracers' physics:

      call tracers_physics()

!-----------------------------------------------------------
!	Radiation

      if(dolongwave.or.doshortwave) then 
	call radiation()     
      end if

!-----------------------------------------------------------
!    Compute diagnostic fields:

      call diagnose()

!----------------------------------------------------------

! Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle	
          
  total_water_after = total_water()
!----------------------------------------------------------
!  collect statistics, write save-file, etc.

   call stepout(nstatsteps)

   ! MPI Ensemble run: turn on mpi after each loop (Song Qiyu, 2022)
   if(dompiensemble) dompi = .true.
  
!----------------------------------------------------------
!----------------------------------------------------------

   ! bloss:  only stop when stat and 2D fields are output.
   !bloss(2016-09-07): if saving local MSE budgets, only stop when they are output
   if ( (nstep.lt.nstop).and.(mod(nstep,nstat).eq.0) .AND. &
        ( (.NOT.save2Davg) .OR. (mod(nstep,nsave2D).eq.0) ) .AND. &
        ( (.NOT.do_chunked_energy_budgets) .OR. (mod(nstep,nsaveMSE).eq.0) ) )then

     if(masterproc) then
       call t_stampf(cputime,usrtime,systime)
       write(*,999) cputime-oldtime, float(nstep-oldstep)*dt/3600.
999    format('CPU TIME = ',f12.4, ' OVER ', f8.2,' MODEL HOURS')
       oldtime = cputime
       oldstep = nstep
     end if

      if(dompi) then
        if(masterproc) itmp1 = nelapse
        call task_bcast_integer(0,itmp1,1)
        nelapse = itmp1(1)
      end if
   end if

end do ! main loop

!----------------------------------------------------------
!----------------------------------------------------------

   call t_stopf('total')
   if(masterproc) call t_prf(rank)

if(masterproc) write(*,*) 'Finished with SAM, exiting...'
call task_stop()

end program crm
