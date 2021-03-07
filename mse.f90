
! peters module mse in mse_budget.f90

!bloss(2016-08): This module for computing grid-point-by-grid-point
!  vertically-integrated FMSE budget quantities that was developed by
!  Matt Peters (circa 2006) will be adapted to give vertically-resolved
!  budgets for FMSE, s_li, q_tot, u and v on a grid that is reduced in
!  the horizontal relative to the model grid.  This setup should 
!  provide both the tools for analyzing convective self-aggregation
!  throutgh budgets of FMSE (Bretherton et al, 2005, JAS) and 
!  FMSE variance (Wing & Emanuel, 2013, JAMES).  In addition, the
!  individual budgets of qtot and sli should ease the construction of
!  single-column model forcings.  The momentum budgets might be
!  useful for analysis of convective momentum transport.  At least,
!  the author hopes that they might...

!bloss(2020-10): Adding option do_chunk_mkbudget.  This outputs chunk
!  averages of two things: all variables in micro_field(:,:,:,:) AND
!  full budgets for individual or combined variables from micro_field().
!  In MICRO_*/microphysics.f90, a few variables must be included and set
!  to enable these budgets:
!
!       - do_chunk_mkbudget (logical): true for budgets to be
!           accumulated and output.
!       - n_mkbudget (integer): # of fields that get budgets computed
!       - flag_mkbudget(nmicro_fields, n_mkbudget), integer array:
!           each row of flag_mkbudget tells which fields within
!           micro_field(:,:,:,:) should be summed for that budget
!           output.  For example, the total water budget QTOT sums
!           all water masses.  The total aerosol budget in
!           MICRO_M2005_PA/ sums the dry aeroosl number (NAd), cloud
!           droplet number (NC) and rain number (NR).
!       - mkbudget_sed(nx,ny,nzm,n_mkbudget), real array: holds
!           sedimentation tendency for fields that have budgets computed
!
!  flag_mkbudget and mkbudget_sed only need to be allocated if
!  do_chunk_mkbudget==.true., but they must be declared in each
!  MICRO_*/microphysics.f90 for SAM to compile.
!
!  Examples: In MICRO_M2005_PA/, this facility will be used to output
!  budgets for aerosol number (dry aerosol NAd, cloud droplet number
!  NC, and rain number NR).  In the isotope-enabled Thompson
!  microphysics, the budgets of water vapor mass mixing ratio for the
!  standard and heavy isotopologues HDO and H218O will be output.

module mse

  use vars, only: u, v, w, t, rho, rhow, pres, p, tabs, qv, qcl, qci, qpl, qpi, &
       dudt, dvdt, dwdt, &
       tabs0, qv0, qn0, qp0, qcldliq0, qcldice0, qpcpliq0, qpcpice0, t0, u0, v0, &
       bet, fluxbq, fluxbt, sstxy, t00
  use grid, only: z, zi, dx, dy, dz, adz, dt, dtn, day0, RUN3D, &
       nx, ny, nzm, nx_gl, ny_gl, nsubdomains_x, nsubdomains_y, nsubdomains, &
       nstep, icycle, ncycle, case, caseid, dogzip3D, &
       dompi, rank, masterproc, nrestart, nstat, nstatis, &
       do_chunked_energy_budgets, do_chunked_momentum_budgets, &
       nchunk_x_gl, nchunk_y_gl, nsaveMSE, na, output_sep

  use params, only: fac_cond, fac_sub, fac_fus, lcond, lsub, lfus, cp, epsv, ug, vg
  use microphysics, only: micro_field, nmicro_fields, flag_wmass, &
       qtot_sed, qice_sed, prec_accum, prec_ice_accum, fluxbmk, &
       do_chunk_mkbudget, n_mkbudget, flag_mkbudget, mkbudget_sed, &
       mkbudget_name, mkbudget_longname, mkbudget_units, &
       mkname, mklongname, mkunits, mkoutputscale, &
       do_mkbudget_extra, n_mkbudget_extra, mkbudget_extra, &
       mkbudget_extra_name, mkbudget_extra_longname, mkbudget_extra_units
  use rad, only: swntxy, swntcxy, swnsxy, swnscxy, &
       lwntxy, lwntcxy, lwnsxy, lwnscxy, &
       qrad_lw, qradclr_lw, qrad_sw, qradclr_sw, &
       qrad, isAllocatedIndividualQrad

  implicit none

  private ! make module variables/subroutines private by default
  public :: initializeMSE, diagnoseMSE, increment_chunk_averages, &
       init_MSE_tendency, compute_and_increment_MSE_tendency, &
       init_momentum_tendency, compute_and_increment_momentum_tendency

  logical, public :: doMSE = .false.
  logical, public :: isAllocatedEnergyBudgets = .false.
  logical, public :: isAllocatedMomentumBudgets = .false.

  integer, public :: nchunk_x, nchunk_y ! number of chunks on a processor in each direction (must be integer!!)
  integer, public :: nx_chunk, ny_chunk ! number of grid points within a chunk in each direction (must be integer)

  integer :: navgMSE, nsaveMSEstart, nsaveMSEend, mseflag

  ! extra variables for computing tendencies elsewhere in the code
  real, dimension(:,:,:), allocatable :: sli_before, fmse_before, qtot_before
  real, dimension(:,:,:), allocatable :: u_before, v_before, w_before

  !  variables to output
  real, dimension(:,:,:), allocatable :: u_mse, v_mse, w_mse   ! m/s
  real, dimension(:,:,:), allocatable :: tabs_mse, tvprime_mse, h_mse, sli_mse ! K
  real, dimension(:,:,:), allocatable :: qtot_mse, qv_mse, qcl_mse ! kg/kg -- Note: qtot includes all water mass
  real, dimension(:,:,:), allocatable :: qci_mse, qpl_mse, qpi_mse ! kg/kg
  real, dimension(:,:,:), allocatable :: rho_mse                 ! [kg/m3]
  real, dimension(:,:,:), allocatable :: pp_mse                 ! [Pa]
  real, dimension(:,:,:), allocatable :: cld_mse, cldcumup_mse, cldcumdn_mse ! cloud fraction, cumulative cloud fraction (up and down)
  real, dimension(:,:,:,:), allocatable :: micro_field_mse

  !  fluxes to output -- computed from anomalies to domain-mean quantities
  real, dimension(:,:,:), allocatable :: u_flux_mse, v_flux_mse, w_flux_mse   ! kg/m/s2
  real, dimension(:,:,:), allocatable :: h_flux_mse, sli_flux_mse ! K
  real, dimension(:,:,:), allocatable :: qtot_flux_mse, qv_flux_mse, qcl_flux_mse ! kg/kg -- Note: qtot includes all water mass
  real, dimension(:,:,:), allocatable :: qci_flux_mse, qpl_flux_mse, qpi_flux_mse ! kg/kg
  real, dimension(:,:,:), allocatable :: buoy_flux_mse                 ! m2/s3

  !  vertically-resolved FMSE, SLI and QTOT budget terms for output in K/s (kg/kg/s for QTOT)
  real, dimension(:,:,:), allocatable :: h_stor_mse, sli_stor_mse, qtot_stor_mse ! K/d -- storage
  real, dimension(:,:,:), allocatable :: h_eddy_mse, sli_eddy_mse, qtot_eddy_mse ! K/d -- resolved advection + subgrid diffusion
  real, dimension(:,:,:), allocatable :: h_sed_mse, sli_sed_mse, qtot_sed_mse, qice_sed_mse ! K/d -- sedimentation
  real, dimension(:,:,:), allocatable :: h_lsf_mse, sli_lsf_mse, qtot_lsf_mse ! K/d -- large-scale forcing
  real, dimension(:,:,:), allocatable :: h_misc_mse, sli_misc_mse, qtot_misc_mse ! K/d -- nudging, damping and upper boundary

  real, dimension(:,:,:), allocatable :: qradclr_lw_mse, qrad_lw_mse ! K/d -- lw radiative htg, clear and full-sky
  real, dimension(:,:,:), allocatable :: qradclr_sw_mse, qrad_sw_mse ! K/d -- sw radiative htg, clear and full-sky

  !  vertically-resolved momentum budget terms for output in m/s2
  real, dimension(:,:,:), allocatable :: u_stor_mse, v_stor_mse, w_stor_mse  ! m/s2 -- storage
  real, dimension(:,:,:), allocatable :: u_eddy_mse, v_eddy_mse, w_eddy_mse  ! m/s2 -- advection + diffusion
  real, dimension(:,:,:), allocatable :: u_pgrad_mse, v_pgrad_mse, w_pgrad_mse  ! m/s2 -- pressure gradient
  real, dimension(:,:,:), allocatable :: u_lsf_mse, v_lsf_mse, w_lsf_mse  ! m/s2 -- large-scale forcing + coriolis
  real, dimension(:,:,:), allocatable :: u_misc_mse, v_misc_mse, w_misc_mse  ! m/s2 -- nudging, damping and upper boundary
  real, dimension(:,:,:), allocatable :: w_buoy_mse  ! m/s2 -- buoyancy

  ! 3D budgets for microphysical quantities
  real, allocatable, dimension(:,:,:,:) :: mkbudget_before, mkbudget_stor_mse, mkbudget_eddy_mse, mkbudget_sed_mse, &
       mkbudget_lsf_mse, mkbudget_misc_mse, mkbudget_mphy_mse, mkbudget_vars_mse, mkbudget_extra_mse
  real, allocatable, dimension(:,:,:) :: mkbudget_surfflux_mse

  !  2D fields to output for budgets
  !  All 2D output fields have units of W/m^2 except for sst=K
  real, dimension(:,:), allocatable :: prec_mse, prec_ice_mse, shf_mse, lhf_mse, sst_mse
  real, dimension(:,:), allocatable :: h_stor_2d_mse, sli_stor_2d_mse, qtot_stor_2d_mse
  real, dimension(:,:), allocatable :: h_eddy_2d_mse, sli_eddy_2d_mse, qtot_eddy_2d_mse
  real, dimension(:,:), allocatable :: h_sed_2d_mse, sli_sed_2d_mse, qtot_sed_2d_mse, qice_sed_2d_mse
  real, dimension(:,:), allocatable :: h_lsf_2d_mse, sli_lsf_2d_mse, qtot_lsf_2d_mse
  real, dimension(:,:), allocatable :: h_misc_2d_mse, sli_misc_2d_mse, qtot_misc_2d_mse
  real, dimension(:,:), allocatable :: lwns_mse, lwnt_mse, lwnsc_mse, lwntc_mse
  real, dimension(:,:), allocatable :: swns_mse, swnt_mse, swnsc_mse, swntc_mse
  real, dimension(:,:), allocatable :: qradclr_lw_2d_mse, qrad_lw_2d_mse ! radiative flux divergence
  real, dimension(:,:), allocatable :: qradclr_sw_2d_mse, qrad_sw_2d_mse 

  ! 1D field to output
  real, dimension(:), allocatable   :: pres_mse, z_mse, zi_mse, layermass_mse

  ! variables used by SAM
  real, dimension(:,:,:), allocatable :: oldsli_mse, oldhf_mse   ! J/kg
  real, dimension(:,:,:), allocatable :: div      ! divergence


  logical :: isAllocatedMSE = .false.

  integer, external :: lenstr

  !--------------------------------------------------------------------
  !  end variable definition.  now start subroutine definition
CONTAINS

subroutine init_MSE_tendency()
  !bloss: This routine records the current values of FMSE, SLI and QTOT
  !  as 3D arrays on the model grid.  These will be used to compute
  !  the increment of these values due to different processes.

  ! NOTE: Here, the negative of all quantities is stored because we are
  !   computing the difference of these quantities across each process
  !   to get the tendency.
  implicit none

  integer :: n, m

  sli_before(:,:,:) = - t(1:nx,1:ny,1:nzm)
  fmse_before(:,:,:) = - t(1:nx,1:ny,1:nzm) 
  qtot_before(:,:,:) = 0.
  do n = 1,nmicro_fields
    if(flag_wmass(n).gt.0) then
      fmse_before(:,:,:) = fmse_before(:,:,:) &
           - fac_cond*micro_field(1:nx,1:ny,1:nzm,n)
      qtot_before(:,:,:) = qtot_before(:,:,:) &
           - micro_field(1:nx,1:ny,1:nzm,n)
    end if
  end do

  if(do_chunk_mkbudget) then
    mkbudget_before(:,:,:,:) = 0.
    do n = 1,n_mkbudget
      do m = 1,nmicro_fields
        if(flag_mkbudget(m,n).gt.0) then
          mkbudget_before(:,:,:,n) = mkbudget_before(:,:,:,n) - micro_field(1:nx,1:ny,1:nzm,m)
        end if
      end do
    end do
  end if

end subroutine init_MSE_tendency

subroutine compute_and_increment_MSE_tendency(which_process)
  !bloss: This routine computes the increment of FMSE, SLI and QTOT
  !  due to a given process and then adds that increment to an array
  !  that accumulates
  implicit none

  character(LEN=*), intent(in) :: which_process

  integer :: n, m

  call t_startf('mse_increment')

  ! compute increment in FMSE, SLI and QTOT since init_MSE_tendency was called
  !   NOTE: We add this to the before value because that already contains the negative sign.
  sli_before(:,:,:) = sli_before(:,:,:) + t(1:nx,1:ny,1:nzm) 
  fmse_before(:,:,:) = fmse_before(:,:,:) + t(1:nx,1:ny,1:nzm)
  
  do n = 1,nmicro_fields
    if(flag_wmass(n).gt.0) then
      fmse_before(:,:,:) = fmse_before(:,:,:) + fac_cond*micro_field(1:nx,1:ny,1:nzm,n)
      qtot_before(:,:,:) = qtot_before(:,:,:) + micro_field(1:nx,1:ny,1:nzm,n)
    end if
  end do

  select case (TRIM(which_process))
    case('eddy')
       call increment_chunk_averages(h_eddy_mse,fmse_before,1.)
       call increment_chunk_averages(sli_eddy_mse,sli_before,1.)
       call increment_chunk_averages(qtot_eddy_mse,qtot_before,1.)
    case('lsf')
       call increment_chunk_averages(h_lsf_mse,fmse_before,1.)
       call increment_chunk_averages(sli_lsf_mse,sli_before,1.)
       call increment_chunk_averages(qtot_lsf_mse,qtot_before,1.)
    case('misc')
       call increment_chunk_averages(h_misc_mse,fmse_before,1.)
       call increment_chunk_averages(sli_misc_mse,sli_before,1.)
       call increment_chunk_averages(qtot_misc_mse,qtot_before,1.)
    case('mphy')
      ! do nothing because the microphysical tendency is zero for
      ! these budgets (though the sedimentation tendency can be
      ! nonzero and is accounted for).
     case default
       write(*,*) 'Error in mse.f90, compute_and_increment_MSE_tendency: bad process type.'
       write(*,*) 'Input process was ', which_process, '.  Available options: eddy, lsf, misc'
       write(*,*) 'Stopping ...'
       call task_abort()
     end select

     !bloss(2020-10): Handle microphysical tendencies separately
     if(do_chunk_mkbudget) then
       do n = 1,n_mkbudget
         do m = 1,nmicro_fields
           if(flag_mkbudget(m,n).gt.0) then
             mkbudget_before(:,:,:,n) = mkbudget_before(:,:,:,n) + micro_field(1:nx,1:ny,1:nzm,m)
           end if
         end do

         select case (TRIM(which_process))
         case('eddy')
           call increment_chunk_averages(mkbudget_eddy_mse(:,:,:,n),mkbudget_before(:,:,:,n),1.)
         case('lsf')
           call increment_chunk_averages(mkbudget_lsf_mse(:,:,:,n),mkbudget_before(:,:,:,n),1.)
         case('misc')
           call increment_chunk_averages(mkbudget_misc_mse(:,:,:,n),mkbudget_before(:,:,:,n),1.)
         case('mphy')
           call increment_chunk_averages(mkbudget_mphy_mse(:,:,:,n),mkbudget_before(:,:,:,n),1.)
         case default
           write(*,*) 'Error in mse.f90, compute_and_increment_MSE_tendency: bad process type.'
           write(*,*) 'Input process was ', which_process, '.  Available options: eddy, lsf, misc'
           write(*,*) 'Stopping ...'
           call task_abort()
         end select
       end do
     end if

  call t_stopf('mse_increment')

end subroutine compute_and_increment_MSE_tendency

subroutine init_momentum_tendency()
  !bloss: This routine records the current values of FMSE, SLI and QTOT
  !  as 3D arrays on the model grid.  These will be used to compute
  !  the increment of these values due to different processes.
  implicit none

  u_before(:,:,:) = dudt(1:nx,1:ny,1:nzm,na) 
  v_before(:,:,:) = dvdt(1:nx,1:ny,1:nzm,na)
  w_before(:,:,:) = 0.5*( dwdt(1:nx,1:ny,1:nzm,na) &
                         + dwdt(1:nx,1:ny,2:nzm+1,na) )

end subroutine init_momentum_tendency

subroutine compute_and_increment_momentum_tendency(which_process)
  !bloss: This routine computes the increment of FMSE, SLI and QTOT
  !  due to a given process and then adds that increment to an array
  !  that accumulates
  implicit none

  character(LEN=*), intent(in) :: which_process

  ! compute increment in dudt, dvdt, dwdt since init_momentum_tendency was called
  u_before(:,:,:) = dudt(1:nx,1:ny,1:nzm,na) - u_before(:,:,:)
  v_before(:,:,:) = dvdt(1:nx,1:ny,1:nzm,na) - v_before(:,:,:)
  w_before(:,:,:) = 0.5*( dwdt(1:nx,1:ny,1:nzm,na) &
                         + dwdt(1:nx,1:ny,2:nzm+1,na) ) &
                         - w_before(:,:,:)

  ! multiply increment in d*dt by dtn to get increment in u/v/w
  ! and add chunk averages to proper budget terms
  select case (TRIM(which_process))
    case('eddy')
       call increment_chunk_averages(u_eddy_mse,u_before,dtn)
       call increment_chunk_averages(v_eddy_mse,v_before,dtn)
       call increment_chunk_averages(w_eddy_mse,w_before,dtn)
    case('pgrad')
       call increment_chunk_averages(u_pgrad_mse,u_before,dtn)
       call increment_chunk_averages(v_pgrad_mse,v_before,dtn)
       call increment_chunk_averages(w_pgrad_mse,w_before,dtn)
    case('lsf')
       call increment_chunk_averages(u_lsf_mse,u_before,dtn)
       call increment_chunk_averages(v_lsf_mse,v_before,dtn)
       call increment_chunk_averages(w_lsf_mse,w_before,dtn)
    case('misc')
       call increment_chunk_averages(u_misc_mse,u_before,dtn)
       call increment_chunk_averages(v_misc_mse,v_before,dtn)
       call increment_chunk_averages(w_misc_mse,w_before,dtn)
    case('buoy')
       call increment_chunk_averages(w_buoy_mse,w_before,dtn)
     case default
       write(*,*) 'Error in mse.f90, compute_and_increment_momentum_tendency: bad process type.'
       write(*,*) 'Input process was ', which_process, '.  Available options: eddy, lsf, misc, pgrad, buoy'
       write(*,*) 'Stopping ...'
       call task_abort()
     end select

end subroutine compute_and_increment_momentum_tendency

subroutine divergence()
 ! computes the divergence from the velocity field
 !  This is called after the velocity field is updated to the next half
 !  step, and u is replaced with u*rho(k)*dtn/dx, and v is replaced with
 !  v*rho(k)*dtn/dy
 !  Don't worry about the w term since d/dz(w*scalar) integrates to zero
 !
 ! on exit, div = du/dx + dv/dy
 
 implicit none

 integer i, j, k, ic, jc
 real coef

 ! computes the divergence from the velocity field
  
 if(RUN3D) then
   do k=1,nzm
     coef = 1./dtn/rho(k)
     do j=1,ny
       jc=j+1 
       do i=1,nx
         ic=i+1
         div(i,j,k) = (u(ic,j,k)-u(i,j,k) + v(i,jc,k)-v(i,j,k))*coef
       end do
     end do
   end do

 else    ! 2D run
   j=1
   do k=1,nzm
     coef = dtn/rho(k)
     do i=1,nx
       ic=i+1
       div(i,j,k) = (u(ic,j,k)-u(i,j,k))*coef
     end do
   end do
 end if

end subroutine divergence

!-------------------------------------------------------------------
subroutine diagnoseMSE()
  ! updates the fields to output, and writes them if necessary

  implicit none

  integer :: i, j, k, n, m
  real :: coef, coef2, omn, omp
  real :: tmp(nx,ny,nzm), tmprhow(nx,ny,nzm), tmp2d(nx,ny)


  call t_startf('mse_diagnose')

  if(masterproc) write(*,*) 'Accumulating budget terms for FMSE/SLI/QTOT in mse.f90 ...'

  ! normalization factor: fraction of total accumulation interval divided by the number of 
  !   grid columns that contribute to a chunk average
  coef = float(nstatis)*dt /( float(nsaveMSE)*dt ) /float( nx_chunk*ny_chunk )
!!$  write(*,*) 'MSE weighting coef = ', coef

  ! store w for fluxes
  do k = 1,nzm
    tmprhow(:,:,k) = 0.5*rhow(k)*w(1:nx,1:ny,k) + 0.5*rhow(k+1)*w(1:nx,1:ny,k+1)
  end do

  ! liquid-ice static energy and its vertical flux
  tmp(:,:,:) = t(1:nx,1:ny,1:nzm)
  call increment_chunk_averages(sli_mse,tmp,coef)

  do k = 1,nzm
    tmp(:,:,k) = cp*tmprhow(:,:,k) *( t(1:nx,1:ny,k)-t0(k) ) ! multiply by rho*cp*w to get approx s_li flux in W/m2
  end do
  call increment_chunk_averages(sli_flux_mse,tmp,coef)

  ! total water and its vertical flux
  tmp(:,:,:) = 0.
  do n = 1,nmicro_fields
    if(flag_wmass(n).gt.0) then
      tmp(:,:,:) = tmp(:,:,:) + micro_field(1:nx,1:ny,1:nzm,n)
    end if
  end do
  call increment_chunk_averages(qtot_mse,tmp,coef)

  do k = 1,nzm
    tmp(:,:,k) = tmprhow(:,:,k) *( tmp(:,:,k) - qv0(k) - qn0(k) - qp0(k) )
  end do
  call increment_chunk_averages(qtot_flux_mse,tmp,coef)

  ! frozen moist static energy and its vertical flux
  h_mse = sli_mse + fac_cond*qtot_mse
  h_flux_mse = sli_flux_mse + lcond*qtot_flux_mse

  if(do_chunk_mkbudget) then
    !bloss(2020-10): output chunk-averaged versions of all micro_field() variables
    do n = 1,nmicro_fields
      tmp(:,:,:) = micro_field(1:nx,1:ny,1:nzm,n)
      call increment_chunk_averages(micro_field_mse(:,:,:,n),tmp,coef)
    end do

    do n = 1,n_mkbudget
      if(SUM(flag_mkbudget(:,n)).gt.1) then
        ! this budget quantity will not be included above, save it here for separate output
        tmp(:,:,:) = 0.
        do m = 1,nmicro_fields
          if(flag_mkbudget(m,n).gt.0) then
            tmp(:,:,:) = tmp(:,:,:) + micro_field(1:nx,1:ny,1:nzm,m)
          end if
        end do
        call increment_chunk_averages(mkbudget_vars_mse(:,:,:,n),tmp,coef)
      end if
    end do
            
  end if


  call increment_chunk_averages(qv_mse,qv,coef) ! water vapor
  do k = 1,nzm
    tmp(:,:,k) = lcond*tmprhow(:,:,k) *( qv(:,:,k) -qv0(k) ) ! multiply by rho*lcond*w to get approx s_li flux in W/m2
  end do
  call increment_chunk_averages(qv_flux_mse,tmp,coef)

  call increment_chunk_averages(qcl_mse,qcl,coef) ! cloud liquid
  do k = 1,nzm
    tmp(:,:,k) = lcond*tmprhow(:,:,k) *( qcl(:,:,k) -qcldliq0(k) )
  end do
  call increment_chunk_averages(qcl_flux_mse,tmp,coef)

  call increment_chunk_averages(qci_mse,qci,coef) ! cloud ice
  do k = 1,nzm
    tmp(:,:,k) = lsub*tmprhow(:,:,k) *( qci(:,:,k) -qcldice0(k) )
  end do
  call increment_chunk_averages(qci_flux_mse,tmp,coef)

  call increment_chunk_averages(qpl_mse,qpl,coef) ! precipitating liquid
  do k = 1,nzm
    tmp(:,:,k) = lcond*tmprhow(:,:,k) *( qpl(:,:,k) -qpcpliq0(k) )
  end do
  call increment_chunk_averages(qpl_flux_mse,tmp,coef)

  call increment_chunk_averages(qpi_mse,qpi,coef) ! precipitating ice
  do k = 1,nzm
    tmp(:,:,k) = lsub*tmprhow(:,:,k) *( qpi(:,:,k) -qpcpice0(k) )
  end do
  call increment_chunk_averages(qpi_flux_mse,tmp,coef)

  tmp(:,:,:) = tabs(:,:,:)
  call increment_chunk_averages(tabs_mse,tmp,coef) ! absolute temperature

  do k = 1,nzm
    tmp(:,:,k) = tabs0(k)*( epsv*(qv(:,:,k)-qv0(k)) &
                           - (qcl(:,:,k) + qci(:,:,k) - qn0(k) &
                             + qpl(:,:,k) + qpi(:,:,k) - qp0(k)) ) &
               + (tabs(:,:,k)-tabs0(k)) *( 1. + epsv*qv0(k) - qn0(k) - qp0(k) )
  end do
  call increment_chunk_averages(tvprime_mse,tmp,coef) ! density temperature anomaly
  do k = 1,nzm
    tmp(:,:,k) = bet(k)*tmprhow(:,:,k)*tmp(:,:,k)/rho(k) ! divide by rho to get w'b' in m2/s3
  end do
  call increment_chunk_averages(buoy_flux_mse,tmp,coef) ! buoyancy flux in m2/s3

  tmp(:,:,:) = ug + 0.5*u(1:nx,1:ny,1:nzm) + 0.5*u(2:nx+1,1:ny,1:nzm)
  call increment_chunk_averages(u_mse,tmp,coef)
  do k = 1,nzm
    tmp(:,:,k) = tmprhow(:,:,k)*( tmp(:,:,k) - u0(k) )
  end do
  call increment_chunk_averages(u_flux_mse,tmp,coef)

  if(RUN3D) then
    tmp(:,:,:) = vg + 0.5*v(1:nx,1:ny,1:nzm) + 0.5*v(1:nx,2:ny+1,1:nzm)
  else
    tmp(:,:,:) = vg + v(1:nx,1:ny,1:nzm) 
  end if
  call increment_chunk_averages(v_mse,tmp,coef)
  do k = 1,nzm
    tmp(:,:,k) = tmprhow(:,:,k)*( tmp(:,:,k) - v0(k) )
  end do
  call increment_chunk_averages(v_flux_mse,tmp,coef)

  tmp(:,:,:) = 0.5*w(1:nx,1:ny,1:nzm) + 0.5*w(1:nx,1:ny,2:nzm+1)
  call increment_chunk_averages(w_mse,tmp,coef)
  do k = 1,nzm
    tmp(:,:,k) = 0.5*rhow(k)*w(1:nx,1:ny,k)*w(1:nx,1:ny,k) &
         + 0.5*rhow(k+1)*w(1:nx,1:ny,k+1)*w(1:nx,1:ny,k+1)
  end do
  call increment_chunk_averages(w_flux_mse,tmp,coef)

  tmp(:,:,:) = p(1:nx,1:ny,1:nzm)
  call increment_chunk_averages(pp_mse,tmp,coef)

  ! cloud fraction: cloud liquid + cloud ice > 1e-5 kg/kg
  tmp(:,:,:) = MAX(0., SIGN(1., qcl(:,:,:) + qci(:,:,:) - 1.e-5) )
  call increment_chunk_averages(cld_mse,tmp,coef)
  
  ! cumulative cloud fraction, going upward from surface
  tmp(:,:,:) = MAX(0., SIGN(1., qcl(:,:,:) + qci(:,:,:) - 1.e-5) )
  do k = 2,nzm
    do j = 1,ny
      do i = 1,nx
        ! will set cloud fraction to one if cloud exists at/below this level
        tmp(i,j,k) = MAX( tmp(i,j,k), tmp(i,j,k-1) )
      end do
    end do
  end do
  call increment_chunk_averages(cldcumup_mse,tmp,coef)

  ! cumulative cloud fraction, going downward from the top of the model
  tmp(:,:,:) = MAX(0., SIGN(1., qcl(:,:,:) + qci(:,:,:) - 1.e-5) )
  do k = nzm-1,1,-1
    do j = 1,ny
      do i = 1,nx
        ! will set cloud fraction to one if cloud exists at/above this level
        tmp(i,j,k) = MAX( tmp(i,j,k), tmp(i,j,k+1) )
      end do
    end do
  end do
  call increment_chunk_averages(cldcumdn_mse,tmp,coef)


  do k = 1,nzm
    ! future-proof in case anelastic assumption is relaxed
    rho_mse(:,:,k) = rho_mse(:,:,k) + float(nx_chunk*ny_chunk)*rho(k)*coef 
  end do

  ! compute tendencies due to radiation
  call increment_chunk_averages(qrad_lw_mse,qrad_lw,coef)
  call increment_chunk_averages(qradclr_lw_mse,qradclr_lw,coef)
  call increment_chunk_averages(qrad_sw_mse,qrad_sw,coef)
  call increment_chunk_averages(qradclr_sw_mse,qradclr_sw,coef)

  ! 2D fluxes
  tmp2d(:,:) = cp*rhow(1)*fluxbt(1:nx,1:ny)
  call increment_chunk_averages2D(shf_mse,tmp2d,coef)

  tmp2d(:,:) = lcond*rhow(1)*fluxbq(1:nx,1:ny)
  call increment_chunk_averages2D(lhf_mse,tmp2d,coef)

  tmp2d(:,:) = sstxy(1:nx,1:ny) + t00 ! include offset for sst
  call increment_chunk_averages2D(sst_mse,tmp2d,coef)

  if(do_chunk_mkbudget) then
    ! First, save the surface fluxes associated with each microphysical
    !   budget term.
    do n = 1,n_mkbudget
      tmp2d(:,:) = 0.
      do m = 1,nmicro_fields
        if(flag_mkbudget(m,n).gt.0) then
          tmp2d(:,:) = tmp2d(:,:) + fluxbmk(1:nx,1:ny,m)
        end if
      end do
      call increment_chunk_averages2D(mkbudget_surfflux_mse(:,:,n),tmp2d,coef)
    end do
  end if

  ! accumulate radiation fluxes.
  call increment_chunk_averages2D( lwns_mse, lwnsxy, coef )
  call increment_chunk_averages2D( lwnt_mse, lwntxy, coef )
  call increment_chunk_averages2D( lwnsc_mse, lwnscxy, coef )
  call increment_chunk_averages2D( lwntc_mse, lwntcxy, coef )

  call increment_chunk_averages2D( swns_mse, swnsxy, coef )
  call increment_chunk_averages2D( swnt_mse, swntxy, coef )
  call increment_chunk_averages2D( swnsc_mse, swnscxy, coef )
  call increment_chunk_averages2D( swntc_mse, swntcxy, coef )

  ! all the fluxes are now accumulated.  if time to write the data, update
  ! the adv, dif fluxes, the storage terms, and write the fields
  if((icycle.eq.ncycle) .and. (mod(nstep,nsaveMSE).eq.0) ) then

    ! finish computing storage: these terms already hold the negative of the initial value
    coef2 = 1./float(nx_chunk*ny_chunk)/( float(nsaveMSE)*dt ) ! normalize sums by total time

    tmp(:,:,:) = t(1:nx,1:ny,1:nzm)
    call increment_chunk_averages(sli_stor_mse,tmp,coef2)

    tmp(:,:,:) = 0.
    do n = 1,nmicro_fields
      if(flag_wmass(n).gt.0) then
        tmp(:,:,:) = tmp(:,:,:) + micro_field(1:nx,1:ny,1:nzm,n)
      end if
    end do
    call increment_chunk_averages(qtot_stor_mse,tmp,coef2)

    h_stor_mse = sli_stor_mse + fac_cond*qtot_stor_mse

    if(do_chunk_mkbudget) then
      !bloss(2020-10): Storage terms for microphysical budgets
      do n = 1,n_mkbudget
        tmp(:,:,:) = 0.
        do m = 1,nmicro_fields
          if(flag_mkbudget(m,n).gt.0) then
            tmp(:,:,:) = tmp(:,:,:) + micro_field(1:nx,1:ny,1:nzm,m)
          end if
        end do
        call increment_chunk_averages(mkbudget_stor_mse(:,:,:,n),tmp,coef2)
      end do
    end if

    ! eddy, lsf and misc budget terms now contain the total increment due to each process
    !  multiplied by the number of columns in each chunk.
    ! Multiplying each by coef2 will turn it into the tendency due to that process.
    h_eddy_mse(:,:,:) = coef2*h_eddy_mse(:,:,:) ! advection + diffusion + surface fluxes
    sli_eddy_mse(:,:,:) = coef2*sli_eddy_mse(:,:,:)
    qtot_eddy_mse(:,:,:) = coef2*qtot_eddy_mse(:,:,:)
    
    h_lsf_mse(:,:,:) = coef2*h_lsf_mse(:,:,:) ! large-scale forcing
    sli_lsf_mse(:,:,:) = coef2*sli_lsf_mse(:,:,:)
    qtot_lsf_mse(:,:,:) = coef2*qtot_lsf_mse(:,:,:)
    
    h_misc_mse(:,:,:) = coef2*h_misc_mse(:,:,:) ! damping + nudging + upperbound
    sli_misc_mse(:,:,:) = coef2*sli_misc_mse(:,:,:)
    qtot_misc_mse(:,:,:) = coef2*qtot_misc_mse(:,:,:)
    
    ! compute tendencies due to sedimentation and average into chunks in x and y.
    call increment_chunk_averages(qtot_sed_mse,qtot_sed,coef2)
    call increment_chunk_averages(qice_sed_mse,qice_sed,coef2)

    ! sli/Cp = T + g*z/Cp - (Lcond/Cp)*(qliq+qice) - (Lfus/Cp)*qice
    sli_sed_mse(:,:,:) = - fac_cond*qtot_sed_mse(:,:,:) - fac_fus*qice_sed_mse(:,:,:)

    ! fmse/Cp = T + g*z/Cp + (Lcond/Cp)*qv - (Lfus/Cp)*qice
    h_sed_mse(:,:,:) = - fac_fus*qice_sed_mse(:,:,:)

    if(do_chunk_mkbudget) then
      !bloss(2020-10): microphysical budgets
      mkbudget_eddy_mse(:,:,:,:) = coef2*mkbudget_eddy_mse(:,:,:,:)
      mkbudget_lsf_mse(:,:,:,:) = coef2*mkbudget_lsf_mse(:,:,:,:)
      mkbudget_misc_mse(:,:,:,:) = coef2*mkbudget_misc_mse(:,:,:,:)
      mkbudget_mphy_mse(:,:,:,:) = coef2*mkbudget_mphy_mse(:,:,:,:)
      do n = 1,n_mkbudget
        call increment_chunk_averages(mkbudget_sed_mse(:,:,:,n),mkbudget_sed(:,:,:,n),coef2)
      end do
    end if

    if(do_mkbudget_extra) then
      do n = 1,n_mkbudget_extra
        call increment_chunk_averages(mkbudget_extra_mse(:,:,:,n),mkbudget_extra(:,:,:,n),coef2)
      end do
    end if

    ! increment precipitation fluxes (in W/m2) here
    tmp2d(:,:) = lcond*dz*prec_accum + lfus*dz*prec_ice_accum
    call increment_chunk_averages2D(prec_mse,tmp2d,coef2)

    tmp2d(:,:) = lsub*dz*prec_ice_accum
    call increment_chunk_averages2D(prec_ice_mse,tmp2d,coef2)

    if(do_chunked_momentum_budgets) then

      tmp(:,:,:) = ug + u(1:nx,1:ny,1:nzm)
      call increment_chunk_averages(u_stor_mse,tmp,coef2)
      
      tmp(:,:,:) = vg + v(1:nx,1:ny,1:nzm) 
      call increment_chunk_averages(v_stor_mse,tmp,coef2)

      tmp(:,:,:) = 0.5*w(1:nx,1:ny,1:nzm) + 0.5*w(1:nx,1:ny,2:nzm+1)
      call increment_chunk_averages(w_stor_mse,tmp,coef2)

      ! eddy, lsf and misc budget terms now contain the total increment due to each process
      !  multiplied by the number of columns in each chunk.
      ! Multiplying each by coef2 will turn it into the tendency due to that process.
      u_eddy_mse(:,:,:) = coef2*u_eddy_mse(:,:,:) ! advection + diffusion + surface fluxes
      v_eddy_mse(:,:,:) = coef2*v_eddy_mse(:,:,:)
      w_eddy_mse(:,:,:) = coef2*w_eddy_mse(:,:,:)

      u_lsf_mse(:,:,:) = coef2*u_lsf_mse(:,:,:) ! large-scale forcing + coriolis forcing
      v_lsf_mse(:,:,:) = coef2*v_lsf_mse(:,:,:)
      w_lsf_mse(:,:,:) = coef2*w_lsf_mse(:,:,:)

      u_misc_mse(:,:,:) = coef2*u_misc_mse(:,:,:) ! damping + nudging + upperbound
      v_misc_mse(:,:,:) = coef2*v_misc_mse(:,:,:)
      w_misc_mse(:,:,:) = coef2*w_misc_mse(:,:,:)

      u_pgrad_mse(:,:,:) = coef2*u_pgrad_mse(:,:,:) ! pressure gradient tendency
      v_pgrad_mse(:,:,:) = coef2*v_pgrad_mse(:,:,:)
      w_pgrad_mse(:,:,:) = coef2*w_pgrad_mse(:,:,:)

      w_buoy_mse(:,:,:) = coef2*w_buoy_mse(:,:,:) ! buoyancy tendency

    end if

    ! compute vertically-integrated budget terms: FMSE
    call columnint(h_stor_mse,h_stor_2d_mse,nchunk_x,nchunk_y)
    call columnint(h_eddy_mse,h_eddy_2d_mse,nchunk_x,nchunk_y)
    call columnint(h_sed_mse,h_sed_2d_mse,nchunk_x,nchunk_y)
    call columnint(h_lsf_mse,h_lsf_2d_mse,nchunk_x,nchunk_y)
    call columnint(h_misc_mse,h_misc_2d_mse,nchunk_x,nchunk_y)
    h_stor_2d_mse(:,:) = cp*h_stor_2d_mse(:,:)
    h_eddy_2d_mse(:,:) = cp*h_eddy_2d_mse(:,:)
    h_sed_2d_mse(:,:) = cp*h_sed_2d_mse(:,:)
    h_lsf_2d_mse(:,:) = cp*h_lsf_2d_mse(:,:)
    h_misc_2d_mse(:,:) = cp*h_misc_2d_mse(:,:)

    ! compute vertically-integrated budget terms: S_LI
    call columnint(sli_stor_mse,sli_stor_2d_mse,nchunk_x,nchunk_y)
    call columnint(sli_eddy_mse,sli_eddy_2d_mse,nchunk_x,nchunk_y)
    call columnint(sli_sed_mse,sli_sed_2d_mse,nchunk_x,nchunk_y)
    call columnint(sli_lsf_mse,sli_lsf_2d_mse,nchunk_x,nchunk_y)
    call columnint(sli_misc_mse,sli_misc_2d_mse,nchunk_x,nchunk_y)
    sli_stor_2d_mse(:,:) = cp*sli_stor_2d_mse(:,:)
    sli_eddy_2d_mse(:,:) = cp*sli_eddy_2d_mse(:,:)
    sli_sed_2d_mse(:,:) = cp*sli_sed_2d_mse(:,:)
    sli_lsf_2d_mse(:,:) = cp*sli_lsf_2d_mse(:,:)
    sli_misc_2d_mse(:,:) = cp*sli_misc_2d_mse(:,:)

    ! compute vertically-integrated budget terms: Q_TOT
    call columnint(qtot_stor_mse,qtot_stor_2d_mse,nchunk_x,nchunk_y)
    call columnint(qtot_eddy_mse,qtot_eddy_2d_mse,nchunk_x,nchunk_y)
    call columnint(qtot_sed_mse,qtot_sed_2d_mse,nchunk_x,nchunk_y)
    call columnint(qtot_lsf_mse,qtot_lsf_2d_mse,nchunk_x,nchunk_y)
    call columnint(qtot_misc_mse,qtot_misc_2d_mse,nchunk_x,nchunk_y)

    call columnint(qice_sed_mse,qice_sed_2d_mse,nchunk_x,nchunk_y)

    ! compute vertically-integrated budget terms: radiative flux divergence
    call columnint(qrad_lw_mse,qrad_lw_2d_mse,nchunk_x,nchunk_y)
    call columnint(qradclr_lw_mse,qradclr_lw_2d_mse,nchunk_x,nchunk_y)
    qrad_lw_2d_mse(:,:) = cp*qrad_lw_2d_mse(:,:)
    qradclr_lw_2d_mse(:,:) = cp*qradclr_lw_2d_mse(:,:)

    call columnint(qrad_sw_mse,qrad_sw_2d_mse,nchunk_x,nchunk_y)
    call columnint(qradclr_sw_mse,qradclr_sw_2d_mse,nchunk_x,nchunk_y)
    qrad_sw_2d_mse(:,:) = cp*qrad_sw_2d_mse(:,:)
    qradclr_sw_2d_mse(:,:) = cp*qradclr_sw_2d_mse(:,:)

    ! write all the data to disk
    call writeMSE()

    ! reset all fields to zero and initialize storage terms
    call resetMSE()

  end if  ! if (icycle.eq.ncycle.and.mod(nstep,nsaveMSE).eq.0), ie, time2write


  call t_stopf('mse_diagnose')

end subroutine diagnoseMSE


!-------------------------------------------------------------------
subroutine docollectMSE()
  ! decide if it is time to collect MSE statistics
  implicit none

  if( (nstep.ge.(nsaveMSEstart-navgMSE+1)).and.(nstep.le.nsaveMSEend).and. &
    ((mod(nstep,nsaveMSE).ge.(nsaveMSE-navgMSE+1)).or. &
                                  (mod(nstep,nsaveMSE).eq.0)) ) then
     doMSE=.true.
  else
     doMSE=.false.
  end if

  if(masterproc) print*,'Starting to collect MSE fluxes...'

end subroutine docollectMSE

!-------------------------------------------------------------------
subroutine writeMSE()

  implicit none
  character *240 filename
  character *80 long_name
  character *8 name
  character *10 timechar
  character *4 rankchar
  character(LEN=5) :: sepchar !bloss: output_sep
  character *6 filetype
  character *10 units
  character *10 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
  integer i,j,k,n
  integer nfields1,nfields2D,nfields3D
  real :: tmp2D(nchunk_x,nchunk_y), tmp2Drho(nchunk_x,nchunk_y), tmp3D(nchunk_x,nchunk_y,nzm)
  real(KIND=4) :: tmp(nchunk_x,nchunk_y,nzm)
  logical :: NotOpened2D = .true., NotOpened3D = .true.

  !bloss: Hard-wire for binary output
  logical, parameter :: MSE3Dbin = .true.

  real, external :: qsatw

  call t_startf('mse_write')

  nfields3D=51 ! many more coming with budget terms... !bloss: 51 = PRECIP
  nfields2D=44 ! many more coming with budget terms...

  if(do_chunk_mkbudget) then
    nfields3D = nfields3D + n_mkbudget*7 + nmicro_fields ! budget terms + eddy fluxes + micro_fields
    nfields2D = nfields2D + n_mkbudget*2 ! surface flux and precipitation

    ! add one more 3D field if each budget that combines different micro_fields
    do n = 1,n_mkbudget
      if(SUM(flag_mkbudget(:,n)).gt.1) then
        nfields3D = nfields3D + 1
      end if
    end do

    if(do_mkbudget_extra) then
      nfields3D = nfields3D + n_mkbudget_extra
    end if
  end if

  if(do_chunked_momentum_budgets) then
    nfields3D = nfields3D + 16
  end if

  !--------------------
  ! dump 3D fields in bin3D file (bin2D for two-dimensional simulation)
  nfields1=0

  if(masterproc.or.output_sep) then

    sepchar=""
    if(output_sep) then
      write(rankchar,'(i4)') rank
      sepchar="_"//rankchar(5-lenstr(rankchar):4)
    end if

    write(rankchar,'(i4.4)') nsubdomains
    write(timechar,'(i10.10)') nstep

    if(RUN3D) then
      filetype = '.bin3D'
      filename='./OUT_3D/'//trim(case)//'_MSE_'//trim(caseid)//'_'// &
           TRIM(rankchar)//'_'//TRIM(timechar)//filetype//sepchar
      open(46,file=filename,status='unknown',form='unformatted')
    else
      filetype = '.bin2D'
      filename='./OUT_3D/'//trim(case)//'_MSE_'//trim(caseid)//'_'// &
           TRIM(rankchar)//filetype//sepchar
      if(nrestart.eq.0.and.NotOpened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                            form='unformatted', position='append')
      end if
      NotOpened3D=.false.
    end if  ! end if(run3D)

    if(masterproc) then
      ! write fields in 3Dbin format
      !  Since we have averaged into chunks, we use nchunk_x and nchunk_y
      write(46) nchunk_x,nchunk_y,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields3D
      do k=1,nzm
        write(46) real(z(k),KIND=4)
      end do
      do k=1,nzm
        write(46) real(pres_mse(k),KIND=4)
      end do
      write(46) real( dx*float(nx_chunk), KIND=4)
      write(46) real( dy*float(ny_chunk), KIND=4)
      write(46) real( float(nstep)*dt/(3600.*24.)+day0, KIND=4)
    end if

  end if         ! if(masterproc)

  ! winds
  tmp(:,:,:) = u_mse(:,:,:); name = 'U'; long_name = 'X Wind Component'; units = 'm/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = v_mse(:,:,:); name = 'V'; long_name = 'Y Wind Component'; units = 'm/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = w_mse(:,:,:); name = 'W'; long_name = 'Z Wind Component'; units = 'm/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  
  ! temperatures
  tmp(:,:,:) = tabs_mse(:,:,:); name = 'TABS'; long_name = 'Absolute temperature'; units = 'K'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = tvprime_mse(:,:,:); name = 'TVPRIME'; long_name = 'Density temperature anomaly'; units = 'K'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = h_mse(:,:,:); name = 'FMSE'; long_name = 'Frozen moist static energy over cp'; units = 'K'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = sli_mse(:,:,:); name = 'SLI'; long_name = 'Liquid-ice static energy over cp'; units = 'K'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! waters
  tmp(:,:,:) = qtot_mse(:,:,:); name = 'QTOT'; long_name = 'Total water mass mixing ratio (vapor+cloud+precip)'; units = 'kg/kg'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qv_mse(:,:,:); name = 'QV'; long_name = 'Water vapor mass mixing ratio'; units = 'kg/kg'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qcl_mse(:,:,:); name = 'QCL'; long_name = 'Cloud liquid mass mixing ratio'; units = 'kg/kg'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qci_mse(:,:,:); name = 'QCI'; long_name = 'Cloud ice mass mixing ratio'; units = 'kg/kg'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qpl_mse(:,:,:); name = 'QPL'; long_name = 'Precipitating liquid mass mixing ratio'; units = 'kg/kg'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qpi_mse(:,:,:); name = 'QPI'; long_name = 'Precipitating ice mass mixing ratio'; units = 'kg/kg'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  
  if(do_chunk_mkbudget) then
    !output chunk-average versions of all microphysical fields
    do n = 1,nmicro_fields
      tmp(:,:,:) = mkoutputscale(n)*micro_field_mse(:,:,:,n); 
      name = TRIM(mkname(n)); long_name = TRIM(mklongname(n)); units = TRIM(mkunits(n)); nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    end do

    do n = 1,n_mkbudget
      if(SUM(flag_mkbudget(:,n)).gt.1) then
        ! Only output budget variables if they are not already output above.
        tmp(:,:,:) = mkbudget_vars_mse(:,:,:,n); 
        name = TRIM(mkbudget_name(n)); long_name = TRIM(mkbudget_longname(n)); 
        units = TRIM(mkbudget_units(n)); nfields1=nfields1+1
        call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      end if
    end do

    if(do_mkbudget_extra) then
      do n = 1,n_mkbudget_extra
        tmp(:,:,:) = mkbudget_extra_mse(:,:,:,n); 
        name = TRIM(mkbudget_extra_name(n)); long_name = TRIM(mkbudget_extra_longname(n)); 
        units = TRIM(mkbudget_extra_units(n)); nfields1=nfields1+1
        call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      end do
    end if

  end if

  ! extras
  tmp(:,:,:) = rho_mse(:,:,:); name = 'RHO'; long_name = 'Density (horizontally-uniform for anelastic SAM)'; 
       units = 'kg/m3'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = pp_mse(:,:,:); name = 'PP'; long_name = 'Perturbation pressure'; units = 'Pa'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! cloud fraction: level-by-level, cumulative up, cumulative down
  tmp(:,:,:) = cld_mse(:,:,:); name = 'CLOUD'; long_name = 'Level-by-level Cloud Fraction'; units = ' '; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = cldcumup_mse(:,:,:); name = 'CLDCUMUP'; long_name = 'Cumulative Cloud Fraction, Upwards From Surface'; units = ' '; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = cldcumdn_mse(:,:,:); name = 'CLDCUMDN'; long_name = 'Cumulative Cloud Fraction, Downwards from Top of Model'; units = ' '; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)


  ! fluxes
  ! winds
  tmp(:,:,:) = u_flux_mse(:,:,:); name = 'RHOWU'; long_name = 'Vertical flux of x momentum'; units = 'kg/m/s2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = v_flux_mse(:,:,:); name = 'RHOWV'; long_name = 'Vertical flux of y momentum'; units = 'kg/m/s2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = w_flux_mse(:,:,:); name = 'RHOWW'; long_name = 'Vertical flux of z momentum'; units = 'kg/m/s2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  
  ! temperatures
  tmp(:,:,:) = buoy_flux_mse(:,:,:); name = 'BUOYFLUX'; long_name = 'Buoyancy flux (Resolved)'; units = 'm2/s3'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = h_flux_mse(:,:,:); name = 'FMSEFLUX'; long_name = 'Frozen moist static energy flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = sli_flux_mse(:,:,:); name = 'SLIFLUX'; long_name = 'Liquid-ice static energy flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! waters
  tmp(:,:,:) = qtot_flux_mse(:,:,:); name = 'QTOTFLUX'; long_name = 'Total water flux (vapor+cloud+precip; Resolved)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qv_flux_mse(:,:,:); name = 'QVFLUX'; long_name = 'Water vapor flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qcl_flux_mse(:,:,:); name = 'QCLFLUX'; long_name = 'Cloud liquid flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qci_flux_mse(:,:,:); name = 'QCIFLUX'; long_name = 'Cloud ice flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qpl_flux_mse(:,:,:); name = 'QPLFLUX'; long_name = 'Precipitating liquid flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qpi_flux_mse(:,:,:); name = 'QPIFLUX'; long_name = 'Precipitating ice flux (Resolved)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  
  ! radiative heating: partitioned into lw/sw, clear-sky/full-sky
  tmp(:,:,:) = qrad_lw_mse(:,:,:); name = 'QRADLW'; long_name = 'Longwave Radiative Heating Rate (Full-sky)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qradclr_lw_mse(:,:,:); name = 'QRADLWCL'; long_name = 'Clearsky longwave Radiative Heating Rate'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qrad_sw_mse(:,:,:); name = 'QRADSW'; long_name = 'Shortwave Radiative Heating Rate (Full-sky)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qradclr_sw_mse(:,:,:); name = 'QRADSWCL'; long_name = 'Clearsky shortwave Radiative Heating Rate'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! FMSE budget terms
  tmp(:,:,:) = h_stor_mse(:,:,:); name = 'FMSESTOR'; long_name = 'Storage of frozen moist static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = h_eddy_mse(:,:,:); name = 'FMSEEDDY'; long_name = 'Advective+diffusive tendency of frozen moist static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = h_lsf_mse(:,:,:); name = 'FMSELSF'; long_name = 'Large-scale forcing tendency of frozen moist static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = h_misc_mse(:,:,:); name = 'FMSEMISC'; long_name = 'Nudging+damping+upperbound tendency of frozen moist static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = h_sed_mse(:,:,:); name = 'FMSESED'; long_name = 'Sedimentation tendency of frozen moist static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! SLI budget terms
  tmp(:,:,:) = sli_stor_mse(:,:,:); name = 'SLISTOR'; long_name = 'Storage of liquid-ice static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = sli_eddy_mse(:,:,:); name = 'SLIEDDY'; long_name = 'Advective+diffusive tendency of liquid-ice static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = sli_lsf_mse(:,:,:); name = 'SLILSF'; long_name = 'Large-scale forcing tendency of liquid-ice static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = sli_misc_mse(:,:,:); name = 'SLIMISC'; long_name = 'Nudging+damping+upperbound tendency of liquid-ice static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = sli_sed_mse(:,:,:); name = 'SLISED'; long_name = 'Sedimentation tendency of liquid-ice static energy (over cp)'; units = 'K/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! QTOT budget terms
  tmp(:,:,:) = qtot_stor_mse(:,:,:); name = 'QTOTSTOR'; long_name = 'Storage of total water (vapor+cloud+precip)'; units = 'kg/kg/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qtot_eddy_mse(:,:,:); name = 'QTOTEDDY'; long_name = 'Advective+diffusive tendency of total water (vapor+cloud+precip)'; units = 'kg/kg/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qtot_lsf_mse(:,:,:); name = 'QTOTLSF'; long_name = 'Large-scale forcing tendency of total water (vapor+cloud+precip)'; units = 'kg/kg/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qtot_misc_mse(:,:,:); name = 'QTOTMISC'; long_name = 'Nudging+damping+upperbound tendency of total water (vapor+cloud+precip)'; units = 'kg/kg/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,:) = qtot_sed_mse(:,:,:); name = 'QTOTSED'; long_name = 'Sedimentation tendency of total water (vapor+cloud+precip)'; units = 'kg/kg/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,:) = qice_sed_mse(:,:,:); name = 'QICESED'; long_name = 'Sedimentation tendency of ice condensate (cloud+precip)'; units = 'kg/kg/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  !bloss(2020-10): Integrate qtot_sed_mse downwards from top of model to yield precipitation rate
  tmp(:,:,:) = 0.
  tmp(:,:,nzm) = - rho(nzm)*dz*adz(nzm)*qtot_sed_mse(:,:,nzm)
  do k = nzm-1,1,-1
    tmp(:,:,k) = tmp(:,:,k+1) - rho(k)*dz*adz(k)*qtot_sed_mse(:,:,k)
  end do
  tmp(:,:,:) = 86400.*tmp(:,:,:) ! convert from kg/m2/s to mm/day
  name = 'PRECIP'; long_name = 'Precipitation rate'; units = 'mm/day'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  !bloss(2020-10): Microphysical budgets for individual species or
  !    the sum of different microphysical species
  if(do_chunk_mkbudget) then
    do n = 1,n_mkbudget
      tmp(:,:,:) = mkbudget_stor_mse(:,:,:,n); name = TRIM(mkbudget_name(n)) // 'STOR'; 
      long_name = 'Storage of ' // TRIM(mkbudget_longname(n)); units = TRIM(mkbudget_units(n)) // '/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      tmp(:,:,:) = mkbudget_eddy_mse(:,:,:,n); name = TRIM(mkbudget_name(n)) // 'EDDY'; 
      long_name = 'Advective+diffusive tendency of ' // TRIM(mkbudget_longname(n)); units = TRIM(mkbudget_units(n)) // '/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      tmp(:,:,:) = mkbudget_lsf_mse(:,:,:,n); name = TRIM(mkbudget_name(n)) // 'LSF'; 
      long_name = 'Large-scale forcing tendency of ' // TRIM(mkbudget_longname(n)); units = TRIM(mkbudget_units(n)) // '/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      tmp(:,:,:) = mkbudget_misc_mse(:,:,:,n); name = TRIM(mkbudget_name(n)) // 'MISC'; 
      long_name = 'Nudging+damping+upperbound tendency of ' // TRIM(mkbudget_longname(n)); units = TRIM(mkbudget_units(n)) // '/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      tmp(:,:,:) = mkbudget_sed_mse(:,:,:,n); name = TRIM(mkbudget_name(n)) // 'SED'; 
      long_name = 'Sedimentation tendency of ' // TRIM(mkbudget_longname(n)); units = TRIM(mkbudget_units(n)) // '/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
      tmp(:,:,:) = mkbudget_mphy_mse(:,:,:,n) - mkbudget_sed_mse(:,:,:,n); name = TRIM(mkbudget_name(n)) // 'MICR'; 
      long_name = 'Microphysical tendency of ' // TRIM(mkbudget_longname(n)); units = TRIM(mkbudget_units(n)) // '/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

      !bloss: Integrate the EDDY tendency to get a flux of each budget term at the interface heights
        !bloss(2020-10): Integrate mkbudget_eddy_mse downwards from top of model to yield precipitation rate
      tmp(:,:,:) = 0.
      tmp(:,:,nzm) = - rho(nzm)*dz*adz(nzm)*mkbudget_eddy_mse(:,:,nzm,n)
      do k = nzm-1,1,-1
        tmp(:,:,k) = tmp(:,:,k+1) - rho(k)*dz*adz(k)*mkbudget_eddy_mse(:,:,k,n)
      end do
      name = TRIM(mkbudget_name(n))//'FLUX'; long_name = 'Resolved + subgrid flux of '//TRIM(mkbudget_longname(n))//' at w-levels';
      units = TRIM(mkbudget_units(n))//'*m/s'; nfields1=nfields1+1
      call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

    end do

  end if

  if(do_chunked_momentum_budgets) then

    ! u budget terms
    tmp(:,:,:) = u_stor_mse(:,:,:); name = 'USTOR'; long_name = 'Storage of zonal velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = u_eddy_mse(:,:,:); name = 'UEDDY'; long_name = 'Advective+diffusive tendency of zonal velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = u_lsf_mse(:,:,:); name = 'ULSF'; long_name = 'Large-scale forcing tendency of zonal velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = u_misc_mse(:,:,:); name = 'UMISC'; long_name = 'Nudging+damping+upperbound tendency of zonal velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = u_pgrad_mse(:,:,:); name = 'UPGRAD'; long_name = 'Pressure gradient tendency of zonal velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

    ! v budget terms
    tmp(:,:,:) = v_stor_mse(:,:,:); name = 'VSTOR'; long_name = 'Storage of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = v_eddy_mse(:,:,:); name = 'VEDDY'; long_name = 'Advective+diffusive tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = v_lsf_mse(:,:,:); name = 'VLSF'; long_name = 'Large-scale forcing tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = v_misc_mse(:,:,:); name = 'VMISC'; long_name = 'Nudging+damping+upperbound tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = v_pgrad_mse(:,:,:); name = 'VPGRAD'; long_name = 'Pressure gradient tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

    ! w budget terms
    tmp(:,:,:) = w_stor_mse(:,:,:); name = 'WSTOR'; long_name = 'Storage of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = w_eddy_mse(:,:,:); name = 'WEDDY'; long_name = 'Advective+diffusive tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = w_lsf_mse(:,:,:); name = 'WLSF'; long_name = 'Large-scale forcing tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = w_misc_mse(:,:,:); name = 'WMISC'; long_name = 'Nudging+damping+upperbound tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = w_pgrad_mse(:,:,:); name = 'WPGRAD'; long_name = 'Pressure gradient tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    tmp(:,:,:) = w_buoy_mse(:,:,:); name = 'WBUOY'; long_name = 'Buoyancy tendency of meridional velocity'; units = 'm/s2'; nfields1=nfields1+1
    call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  end if

  ! finished outputing 3D fields
  call task_barrier()

  if(nfields3D.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: expected nfields = ', nfields3D, ', actual nfields = ', nfields1
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting MSE data. file:'//filename
    else
       print*, 'Appending MSE data. file:'//filename
    end if
  endif

  !--------------------
  ! dump separate 2Dbin file for surface/top-of-atmosphere/top-of-model fluxes and 
  !   vertically-integrated quantities
  !--------------------
  ! dump 2D fields in 2Dbin file
  nfields1=0

  if(masterproc.or.output_sep) then

    sepchar=""
    if(output_sep) then
      write(rankchar,'(i4)') rank
      sepchar="_"//rankchar(5-lenstr(rankchar):4)
    end if

    write(rankchar,'(i4.4)') nsubdomains
    write(timechar,'(i10.10)') nstep

    filetype = '.2Dbin'
    filename='./OUT_2D/'//trim(case)//'_MSE_'//trim(caseid)//'_'// &
         TRIM(rankchar)//filetype//sepchar
    if(nrestart.eq.0.and.NotOpened2D) then
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      open(46,file=filename,status='unknown', &
           form='unformatted', position='append')
    end if
    NotOpened2D=.false.

    ! write fields in 2Dbin format
    !  Since we have averaged into chunks, we use nchunk_x and nchunk_y
    write(46) nstep
    write(46) nchunk_x,nchunk_y,nzm,nsubdomains,nsubdomains_x,nsubdomains_y, nfields2D
    write(46) real( dx*float(nx_chunk), KIND=4)
    write(46) real( dy*float(ny_chunk), KIND=4)
    write(46) real( float(nstep)*dt/(3600.*24.)+day0, KIND=4)

  end if         ! if(masterproc)

  ! re-initialize single precision 3D array to zero.
  !   Note: only first level will be used to output 2D fields
  tmp(:,:,:) = 0.

  ! surface fluxes
  tmp(:,:,1)=prec_mse; name='PREC'; long_name='Surface Precip. Rate'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=prec_ice_mse; name='PRECICE'; long_name='Surface Ice Precip. Rate'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=shf_mse; name='SHF'; long_name='Sensible Heat Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=lhf_mse; name='LHF'; long_name='Latent Heat Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=sst_mse; name='SST'; long_name='Sea Surface Temperature'; units='K'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  if(do_chunk_mkbudget) then
    ! First output surface flux
    do n = 1,n_mkbudget
      tmp(:,:,1)=mkbudget_surfflux_mse(:,:,n); 
      name=TRIM(mkbudget_name(n)) // 'SFLX'; 
      long_name='Surface flux of '//TRIM(mkbudget_longname(n)); 
      units=TRIM(mkbudget_units(n)) // '*m/s'; 
      nfields1 = nfields1 + 1
      call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    end do

    ! next, output precipitation by integrating sedimentation tendency down from the model top.
    do n = 1,n_mkbudget
      tmp(:,:,1)=0.
      do k = nzm,1,-1
        tmp(:,:,1) = tmp(:,:,1) - rho(k)*dz*adz(k)*mkbudget_sed_mse(:,:,k,n)
      end do

      name=TRIM(mkbudget_name(n)) // 'PREC'; 
      long_name='Surface precipitation rate for '//TRIM(mkbudget_longname(n)); 
      units=TRIM(mkbudget_units(n)) // '*kg/m2/s'; 
      nfields1 = nfields1 + 1
      call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
    end do

  end if

  !--------------------
  ! MSE budget stuff

  !------------------
  ! vertically-integrated FMSE, SLI, WVP, etc.

  ! save column mass
  tmp3D(:,:,:) = 1. !bloss(2016-11-23): Note that factor of rho is included in columnint routine.
  call columnint(tmp3D,tmp2Drho,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2Drho(:,:); name='MASS'; long_name='Column mass'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  !bloss(2016-11-23): Remove factor of rho from argument of columnint since factor of rho is included inside columnint
  call columnint(h_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:)/tmp2Drho(:,:); name='FMSE'; long_name='Density-weighted frozen moist static energy (over Cp)'; units='K'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(sli_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:)/tmp2Drho(:,:); name='SLI'; long_name='Density-weighted liquid-ice static energy (over Cp)'; units='K'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(qtot_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='TWP'; long_name='Total water path (vapor+cloud+precip)'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(qv_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='WVP'; long_name='Water vapor path'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(qcl_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='CLIQWP'; long_name='Cloud liquid water path'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(qci_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='CICEWP'; long_name='Cloud ice water path'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(qpl_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='PLIQWP'; long_name='Precipitating liquid water path'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  call columnint(qpi_mse,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='PICEWP'; long_name='Precipitating ice water path'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  do k = 1,nzm
    do j = 1,nchunk_y
      do i = 1,nchunk_x
        tmp3D(i,j,k) = qsatw(tabs_mse(i,j,k),pres_mse(k))
      end do
    end do
  end do
  call columnint(tmp3D,tmp2D,nchunk_x,nchunk_y)
  tmp(:,:,1)=tmp2D(:,:); name='SWVP'; long_name='Saturation water vapor path (\int rho q_sat(p,T) dz)'; units='kg/m2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! Cloud fraction (at least one level has qcl+qci>1e-5 kg/kg)
  tmp(:,:,1)=cldcumup_mse(:,:,nzm); name='CLDSHD'; long_name='Shaded Cloud Fraction (a cloudy column has at least one grid level w/qcl+qci>1e-5 kg/kg)'; units=' '; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  



  !------------------
  ! radiative fluxes
  tmp(:,:,1)=lwns_mse; name='LWNS'; long_name='Net Surface Upward LW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=lwnt_mse; name='LWNT'; long_name='Net TOA Upward LW Flux (OLR)'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=lwnsc_mse; name='LWNSC'; long_name='Clear Sky Surface Upward LW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=lwntc_mse; name='LWNTC'; long_name='Clear Sky TOA Upward LW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=swns_mse; name='SWNS'; long_name='Net Surface Downward SW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  tmp(:,:,1)=swnt_mse; name='SWNT'; long_name='Net TOA Downward SW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1)=swnsc_mse; name='SWNSC'; long_name='Clear Sky Surface Downward SW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
                                                                                
  tmp(:,:,1)=swntc_mse; name='SWNTC'; long_name='Clear Sky TOA Downward SW Flux'; units='W/m^2'; nfields1 = nfields1 + 1
  call compress3D(tmp,nchunk_x,nchunk_y,1,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! radiative flux divergence: partitioned into lw/sw, clear-sky/full-sky
  tmp(:,:,1) = qrad_lw_2d_mse(:,:); name = 'DRADLW'; long_name = 'Longwave Radiative Flux Divergence (Full-sky)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qradclr_lw_2d_mse(:,:); name = 'DRADLWCL'; long_name = 'Clearsky longwave Radiative Flux Divergence'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qrad_sw_2d_mse(:,:); name = 'DRADSW'; long_name = 'Shortwave Radiative Flux Divergence (Full-sky)'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qradclr_sw_2d_mse(:,:); name = 'DRADSWCL'; long_name = 'Clearsky shortwave Radiative Flux Divergence'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! vertically-integrated FMSE budget terms
  tmp(:,:,1) = h_stor_2d_mse(:,:); name = 'FMSESTOR'; long_name = 'Storage of column-integrated frozen moist static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = h_eddy_2d_mse(:,:); name = 'FMSEEDDY'; long_name = 'Advective+diffusive tendency of column-integrated frozen moist static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = h_lsf_2d_mse(:,:); name = 'FMSELSF'; long_name = 'Large-scale forcing tendency of column-integrated frozen moist static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = h_misc_2d_mse(:,:); name = 'FMSEMISC'; long_name = 'Nudging+damping+upperbound of column-integrated frozen moist static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = h_sed_2d_mse(:,:); name = 'FMSESED'; long_name = 'Sedimentation tendency of column-integrated frozen moist static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1) = sli_stor_2d_mse(:,:); name = 'SLISTOR'; long_name = 'Storage of column-integrated liquid-ice static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = sli_eddy_2d_mse(:,:); name = 'SLIEDDY'; long_name = 'Advective+diffusive tendency of column-integrated liquid-ice static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = sli_lsf_2d_mse(:,:); name = 'SLILSF'; long_name = 'Large-scale forcing tendency of column-integrated liquid-ice static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = sli_misc_2d_mse(:,:); name = 'SLIMISC'; long_name = 'Nudging+damping+upperbound of column-integrated liquid-ice static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = sli_sed_2d_mse(:,:); name = 'SLISED'; long_name = 'Sedimentation tendency of column-integrated liquid-ice static energy'; units = 'W/m2'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1) = qtot_stor_2d_mse(:,:); name = 'QTOTSTOR'; long_name = 'Storage of column-integrated total water (vapor+cloud+precip)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qtot_eddy_2d_mse(:,:); name = 'QTOTEDDY'; long_name = 'Advective+diffusive tendency of column-integrated total water (vapor+cloud+precip)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qtot_lsf_2d_mse(:,:); name = 'QTOTLSF'; long_name = 'Large-scale forcing tendency of column-integrated total water (vapor+cloud+precip)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qtot_misc_2d_mse(:,:); name = 'QTOTMISC'; long_name = 'Nudging+damping+upperbound of column-integrated total water (vapor+cloud+precip)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)
  tmp(:,:,1) = qtot_sed_2d_mse(:,:); name = 'QTOTSED'; long_name = 'Sedimentation tendency of column-integrated total water (vapor+cloud+precip)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  tmp(:,:,1) = qice_sed_2d_mse(:,:); name = 'QICESED'; long_name = 'Sedimentation tendency of column-integrated ice condensate (cloud+precip)'; units = 'kg/m2/s'; nfields1=nfields1+1
  call compress3D(tmp,nchunk_x,nchunk_y,nzm,name,long_name,units,MSE3Dbin,dompi,rank,nsubdomains)

  ! finished outputing 2D fields
  call task_barrier()

  if(nfields2D.ne.nfields1) then
    if(masterproc) print*,'write_fields2D error: expected nfields = ', nfields2D, ', actual nfields = ', nfields1
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting MSE data. file:'//filename
    else
       print*, 'Appending MSE data. file:'//filename
    end if
  endif

  ! reset MSE fields for a new averaging period
  call resetMSE()

  call t_stopf('mse_write')

end subroutine writeMSE


!-------------------------------------------------------------------
subroutine columnint(f,value,ncolx,ncoly)
  !  mass weighted column integral of scalar field
  !        returns int_0^TOA (rho*f dz)


implicit none

real, intent(in)  :: f(1:ncolx,1:ncoly,1:nzm)
real, intent(out) :: value(1:ncolx,1:ncoly)
integer, intent(in) :: ncolx, ncoly

integer k, i, j

!-----------------------------
value(:,:) = 0.
do k=1,nzm
  value(:,:) = value(:,:) + f(:,:,k)*(rho(k)*dz*adz(k))
end do

end subroutine columnint


!-------------------------------------------------------------------
subroutine initializeMSE()
  ! compute number of chunks on each processor, allocate all the MSE variables
  !    and set their initial values to 0

  implicit none

  integer :: ierr

  ! error checking
  if(MOD(nchunk_x_gl,nsubdomains_x).ne.0) then
    write(*,*) ' ERROR in initializeMSE within mse.f90: nchunk_x_gl must be an even multiple of nsubdomains_x'
    call task_abort()
  end if

  if(MOD(nchunk_y_gl,nsubdomains_y).ne.0) then
    write(*,*) ' ERROR in initializeMSE within mse.f90: nchunk_y_gl must be an even multiple of nsubdomains_y'
    call task_abort()
  end if

  if(MOD(nx_gl,nchunk_x_gl).ne.0) then
    write(*,*) ' ERROR in initializeMSE within mse.f90: nx_gl must be an even multiple of nchunk_x_gl'
    call task_abort()
  end if

  if(MOD(ny_gl,nchunk_y_gl).ne.0) then
    write(*,*) ' ERROR in initializeMSE within mse.f90: ny_gl must be an even multiple of nchunk_y_gl'
    call task_abort()
  end if

  if(do_chunked_momentum_budgets.AND.(.NOT.do_chunked_energy_budgets)) then
    write(*,*) ' ERROR in initializeMSE within mse.f90: cannot output momentum budgets without energy budgets'
    call task_abort()
  end if

  if(MOD(nsaveMSE,nstat).ne.0) then
    write(*,*) ' ERROR in initializeMSE within mse.f90: nsaveMSE should be an even multiple of nstat.'
    call task_abort()
  end if

  if(do_chunk_mkbudget) then
    if(n_mkbudget.lt.1) then
      write(*,*) 'ERROR in initializeMSE within mse.f90: If do_chunk_mkbudget==.true., n_mkbudget should be positive'
      call task_abort()
    end if
    if(n_mkbudget.gt.nmicro_fields) then
      write(*,*) 'ERROR in initializeMSE within mse.f90: If do_chunk_mkbudget==.true., n_mkbudget should be smaller than nmicro_fields'
      call task_abort()
    end if
  end if

  ! compute number of chunks on each processor
  nchunk_x = nchunk_x_gl/nsubdomains_x
  nchunk_y = nchunk_y_gl/nsubdomains_y

  ! compute number of grid points in each chunk
  nx_chunk = nx/nchunk_x
  ny_chunk = ny/nchunk_y
  
  if(do_chunked_energy_budgets) then

    allocate(  sli_before(nx,ny,nzm), fmse_before(nx,ny,nzm), qtot_before(nx,ny,nzm), &
         u_before(nx,ny,nzm), v_before(nx,ny,nzm), w_before(nx,ny,nzm), &
         u_mse(nchunk_x,nchunk_y,nzm), v_mse(nchunk_x,nchunk_y,nzm), w_mse(nchunk_x,nchunk_y,nzm), &
         tabs_mse(nchunk_x,nchunk_y,nzm), tvprime_mse(nchunk_x,nchunk_y,nzm), &
         h_mse(nchunk_x,nchunk_y,nzm), sli_mse(nchunk_x,nchunk_y,nzm), &
         qtot_mse(nchunk_x,nchunk_y,nzm), qv_mse(nchunk_x,nchunk_y,nzm), qcl_mse(nchunk_x,nchunk_y,nzm), &
         qci_mse(nchunk_x,nchunk_y,nzm), qpl_mse(nchunk_x,nchunk_y,nzm), qpi_mse(nchunk_x,nchunk_y,nzm), &
         rho_mse(nchunk_x,nchunk_y,nzm), pp_mse(nchunk_x,nchunk_y,nzm), &
         cld_mse(nchunk_x,nchunk_y,nzm), cldcumup_mse(nchunk_x,nchunk_y,nzm), cldcumdn_mse(nchunk_x,nchunk_y,nzm), &
         STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot allocate mean fields'
      call task_abort()
    end if

    !  fluxes to output -- computed from anomalies to domain-mean quantities
    allocate( u_flux_mse(nchunk_x,nchunk_y,nzm), v_flux_mse(nchunk_x,nchunk_y,nzm), w_flux_mse(nchunk_x,nchunk_y,nzm), &
         buoy_flux_mse(nchunk_x,nchunk_y,nzm), h_flux_mse(nchunk_x,nchunk_y,nzm), sli_flux_mse(nchunk_x,nchunk_y,nzm) , &
         qtot_flux_mse(nchunk_x,nchunk_y,nzm), qv_flux_mse(nchunk_x,nchunk_y,nzm), qcl_flux_mse(nchunk_x,nchunk_y,nzm), &
         qci_flux_mse(nchunk_x,nchunk_y,nzm), qpl_flux_mse(nchunk_x,nchunk_y,nzm), qpi_flux_mse(nchunk_x,nchunk_y,nzm), &
         STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot allocate fluxes'
      call task_abort()
    end if
    
    !  vertically-resolved FMSE, SLI and QTOT budget terms for output in K/s (kg/kg/s for QTOT)
    allocate( h_stor_mse(nchunk_x,nchunk_y,nzm), sli_stor_mse(nchunk_x,nchunk_y,nzm), qtot_stor_mse(nchunk_x,nchunk_y,nzm), &
         h_eddy_mse(nchunk_x,nchunk_y,nzm), sli_eddy_mse(nchunk_x,nchunk_y,nzm), qtot_eddy_mse(nchunk_x,nchunk_y,nzm), &
         h_sed_mse(nchunk_x,nchunk_y,nzm), sli_sed_mse(nchunk_x,nchunk_y,nzm), &
         qtot_sed_mse(nchunk_x,nchunk_y,nzm), qice_sed_mse(nchunk_x,nchunk_y,nzm), &
         h_lsf_mse(nchunk_x,nchunk_y,nzm), sli_lsf_mse(nchunk_x,nchunk_y,nzm), qtot_lsf_mse(nchunk_x,nchunk_y,nzm), &
         h_misc_mse(nchunk_x,nchunk_y,nzm), sli_misc_mse(nchunk_x,nchunk_y,nzm), qtot_misc_mse(nchunk_x,nchunk_y,nzm) , &
         qradclr_lw_mse(nchunk_x,nchunk_y,nzm), qrad_lw_mse(nchunk_x,nchunk_y,nzm), &
         qradclr_sw_mse(nchunk_x,nchunk_y,nzm), qrad_sw_mse(nchunk_x,nchunk_y,nzm), STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot allocate arrays for 3D budgets'
      call task_abort()
    end if
    
    !  2D fields to output for budgets
    !  All 2D output fields have units of W/m^2 except for sst=K
    allocate( prec_mse(nchunk_x,nchunk_y), prec_ice_mse(nchunk_x,nchunk_y), shf_mse(nchunk_x,nchunk_y), lhf_mse(nchunk_x,nchunk_y), sst_mse(nchunk_x,nchunk_y), &
         h_stor_2d_mse(nchunk_x,nchunk_y), sli_stor_2d_mse(nchunk_x,nchunk_y), qtot_stor_2d_mse(nchunk_x,nchunk_y), &
         h_eddy_2d_mse(nchunk_x,nchunk_y), sli_eddy_2d_mse(nchunk_x,nchunk_y), qtot_eddy_2d_mse(nchunk_x,nchunk_y), &
         h_sed_2d_mse(nchunk_x,nchunk_y), sli_sed_2d_mse(nchunk_x,nchunk_y), &
         qtot_sed_2d_mse(nchunk_x,nchunk_y), qice_sed_2d_mse(nchunk_x,nchunk_y), &
         h_lsf_2d_mse(nchunk_x,nchunk_y), sli_lsf_2d_mse(nchunk_x,nchunk_y), qtot_lsf_2d_mse(nchunk_x,nchunk_y), &
         h_misc_2d_mse(nchunk_x,nchunk_y), sli_misc_2d_mse(nchunk_x,nchunk_y), qtot_misc_2d_mse(nchunk_x,nchunk_y), &
         qrad_lw_2d_mse(nchunk_x,nchunk_y), qradclr_lw_2d_mse(nchunk_x,nchunk_y), &
         qrad_sw_2d_mse(nchunk_x,nchunk_y), qradclr_sw_2d_mse(nchunk_x,nchunk_y), &
         lwns_mse(nchunk_x,nchunk_y), lwnt_mse(nchunk_x,nchunk_y), lwnsc_mse(nchunk_x,nchunk_y), lwntc_mse(nchunk_x,nchunk_y), &
         swns_mse(nchunk_x,nchunk_y), swnt_mse(nchunk_x,nchunk_y), swnsc_mse(nchunk_x,nchunk_y), swntc_mse(nchunk_x,nchunk_y), STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot allocate arrays for vertically-integrated budget terms'
      call task_abort()
    end if

    if(do_chunked_momentum_budgets) then
      !  vertically-resolved momentum budget terms for output in m/s2
      allocate( u_stor_mse(nchunk_x,nchunk_y,nzm), v_stor_mse(nchunk_x,nchunk_y,nzm), w_stor_mse(nchunk_x,nchunk_y,nzm), &
           u_eddy_mse(nchunk_x,nchunk_y,nzm), v_eddy_mse(nchunk_x,nchunk_y,nzm), w_eddy_mse(nchunk_x,nchunk_y,nzm), &
           u_pgrad_mse(nchunk_x,nchunk_y,nzm), v_pgrad_mse(nchunk_x,nchunk_y,nzm), w_pgrad_mse(nchunk_x,nchunk_y,nzm), &
           u_lsf_mse(nchunk_x,nchunk_y,nzm), v_lsf_mse(nchunk_x,nchunk_y,nzm), w_lsf_mse(nchunk_x,nchunk_y,nzm), &
           u_misc_mse(nchunk_x,nchunk_y,nzm), v_misc_mse(nchunk_x,nchunk_y,nzm), w_misc_mse(nchunk_x,nchunk_y,nzm), &
           w_buoy_mse(nchunk_x,nchunk_y,nzm), STAT=ierr)
      if(ierr.ne.0) then
        write(*,*) ' ERROR in initializeMSE within mse.f90: cannot allocate fluxes'
        call task_abort()
      end if
    end if

    ! 1D field to output
    allocate(pres_mse(nzm), z_mse(nzm), zi_mse(nzm), layermass_mse(nzm), STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot allocate 1D quantities'
      call task_abort()
    end if

    isAllocatedMSE = .true.

    if(do_chunked_energy_budgets.AND.(.NOT.isAllocatedIndividualQrad)) then
      ! allocate individual arrays with lw/sw, full sky/clear sky heating rates
      !  these will be averaged horizontally into chunks and output in mse.f90.
      allocate(qrad_lw(nx,ny,nzm), qradclr_lw(nx,ny,nzm), qrad_sw(nx,ny,nzm), qradclr_sw(nx,ny,nzm), STAT=ierr)
      if(ierr.ne.0) then
        write(*,*) 'Cannot allocate individual qrad arrays in initialize_radiation'
        call task_abort()
      end if
      isAllocatedIndividualQrad = .true.

      ! just in case the run does not use interactive radiation, zero out fluxes and heating rates
      swntxy(:,:) = 0.
      swntcxy(:,:) = 0.
      swnsxy(:,:) = 0.
      swnscxy(:,:) = 0.
      lwntxy(:,:) = 0.
      lwntcxy(:,:) = 0.
      lwnsxy(:,:) = 0.
      lwnscxy(:,:) = 0.
      qrad_lw(:,:,:) = 0.
      qradclr_lw(:,:,:) = 0.
      qrad_sw(:,:,:) = 0.
      qradclr_sw(:,:,:) = 0.
      qrad(:,:,:) = 0.
    end if

    if(do_chunk_mkbudget.and.(n_mkbudget.gt.0)) then
      allocate(  micro_field_mse(nchunk_x,nchunk_y,nzm,nmicro_fields), &
           mkbudget_extra_mse(nchunk_x,nchunk_y,nzm,n_mkbudget_extra), &
           mkbudget_vars_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_before(nx,ny,nzm,n_mkbudget), &
           mkbudget_stor_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_eddy_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_sed_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_lsf_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_misc_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_mphy_mse(nchunk_x,nchunk_y,nzm,n_mkbudget), &
           mkbudget_surfflux_mse(nchunk_x,nchunk_y,n_mkbudget), STAT=ierr) 
      if(ierr.ne.0) STOP 'when allocating mkbudget arrays in mse.f90'
    end if

  end if

  ! initialize all fields to zero
  call resetMSE()

end subroutine initializeMSE

!--------------------------------------------------------------------------
subroutine resetMSE()
  ! set all budget-related fields to zero.
  ! initialize storage terms with negative of current value
  !   for each quantity, so that they can be differenced at end
  !   to give storage tendency.

  implicit none

  real :: tmp(nx,ny,nzm), coef2
  integer :: n, m

  u_mse=0.
  v_mse=0.
  w_mse=0.
  tabs_mse=0.
  tvprime_mse=0.
  h_mse=0.
  sli_mse=0.
  qtot_mse=0.
  qv_mse=0.
  qcl_mse=0.
  qci_mse=0.
  qpl_mse=0.
  qpi_mse=0.
  rho_mse=0.
  pp_mse=0.
  cld_mse=0.
  cldcumup_mse=0.
  cldcumdn_mse=0.

  !  fluxes to output -- computed from anomalies to domain-mean quantities
  u_flux_mse=0.
  v_flux_mse=0.
  w_flux_mse=0.
  buoy_flux_mse=0.
  h_flux_mse=0.
  sli_flux_mse = 0.
  qtot_flux_mse=0.
  qv_flux_mse=0.
  qcl_flux_mse=0.
  qci_flux_mse=0.
  qpl_flux_mse=0.
  qpi_flux_mse=0.

  !  vertically-resolved FMSE, SLI and QTOT budget terms for output in K/s (kg/kg/s for QTOT)
  h_stor_mse=0.
  sli_stor_mse=0.
  qtot_stor_mse=0.
  h_eddy_mse=0.
  sli_eddy_mse=0.
  qtot_eddy_mse=0.
  h_sed_mse=0.
  sli_sed_mse=0.
  qtot_sed_mse=0.
  qice_sed_mse=0.
  h_lsf_mse=0.
  sli_lsf_mse=0.
  qtot_lsf_mse=0.
  h_misc_mse=0.
  sli_misc_mse=0.
  qtot_misc_mse = 0.
  qradclr_lw_mse=0.
  qrad_lw_mse=0.
  qradclr_sw_mse=0.
  qrad_sw_mse=0.
  
  if(do_chunk_mkbudget) then
    ! microphysical budgets
    micro_field_mse = 0.
    mkbudget_vars_mse = 0.
    mkbudget_extra_mse=0.
    mkbudget_surfflux_mse = 0.
    mkbudget_stor_mse=0.
    mkbudget_eddy_mse=0.
    mkbudget_sed_mse=0.
    mkbudget_lsf_mse=0.
    mkbudget_misc_mse=0.
    mkbudget_mphy_mse=0.
  end if

  !  2D fields to output for budgets
  prec_mse=0.
  prec_ice_mse=0.
  shf_mse=0.
  lhf_mse=0.
  sst_mse=0.
  h_stor_2d_mse=0.
  sli_stor_2d_mse=0.
  qtot_stor_2d_mse=0.
  h_eddy_2d_mse=0.
  sli_eddy_2d_mse=0.
  qtot_eddy_2d_mse=0.
  h_sed_2d_mse=0.
  sli_sed_2d_mse=0.
  qtot_sed_2d_mse=0.
  qice_sed_2d_mse=0.
  h_lsf_2d_mse=0.
  sli_lsf_2d_mse=0.
  qtot_lsf_2d_mse=0.
  h_misc_2d_mse=0.
  sli_misc_2d_mse=0.
  qtot_misc_2d_mse=0.
  lwns_mse=0.
  lwnt_mse=0.
  lwnsc_mse=0.
  lwntc_mse=0.
  swns_mse=0.
  swnt_mse=0.
  swnsc_mse=0.
  swntc_mse=0.
  qrad_lw_2d_mse=0.
  qradclr_lw_2d_mse=0.
  qrad_sw_2d_mse=0.
  qradclr_sw_2d_mse=0.

  if(do_chunked_momentum_budgets) then
    u_stor_mse=0.
    v_stor_mse=0.
    w_stor_mse=0.
    u_eddy_mse=0.
    v_eddy_mse=0.
    w_eddy_mse=0.
    u_pgrad_mse=0.
    v_pgrad_mse=0.
    w_pgrad_mse=0.
    u_lsf_mse=0.
    v_lsf_mse=0.
    w_lsf_mse=0.
    u_misc_mse=0.
    v_misc_mse=0.
    w_misc_mse=0.
    w_buoy_mse=0.
  end if

  ! initialize time-invariant quantities to their values
  pres_mse(:) = pres(:)
  z_mse(1:nzm) = z(1:nzm)
  zi_mse = zi(1:nzm)
  layermass_mse(:) = rho(:)*dz*adz(:)

  ! initialize storage terms
  ! NOTE NEGATIVE SIGN: WE ARE SUBTRACTING THIS OFF THE FINAL VALUE
  coef2 = - 1./float(nx_chunk*ny_chunk) /( float(nsaveMSE)*dt ) ! normalize sums by total time

  tmp(:,:,:) = t(1:nx,1:ny,1:nzm)
  call increment_chunk_averages(sli_stor_mse,tmp,coef2)

  tmp(:,:,:) = 0.
  do n = 1,nmicro_fields
    if(flag_wmass(n).gt.0) then
      tmp(:,:,:) = tmp(:,:,:) + micro_field(1:nx,1:ny,1:nzm,n)
    end if
  end do
  call increment_chunk_averages(qtot_stor_mse,tmp,coef2)

  h_stor_mse = sli_stor_mse + fac_cond*qtot_stor_mse

  if(do_chunk_mkbudget) then
    do n = 1,n_mkbudget
      tmp(:,:,:) = 0.
      do m = 1,nmicro_fields
        if(flag_mkbudget(m,n).gt.0) then
          tmp(:,:,:) = tmp(:,:,:) + micro_field(1:nx,1:ny,1:nzm,m)
        end if
      end do
      call increment_chunk_averages(mkbudget_stor_mse(:,:,:,n),tmp,coef2)
    end do
  end if

  if(do_chunked_momentum_budgets) then

    tmp(:,:,:) = ug + u(1:nx,1:ny,1:nzm)
    call increment_chunk_averages(u_stor_mse,tmp,coef2)

    tmp(:,:,:) = vg + v(1:nx,1:ny,1:nzm) 
    call increment_chunk_averages(v_stor_mse,tmp,coef2)

    tmp(:,:,:) = 0.5*w(1:nx,1:ny,1:nzm) + 0.5*w(1:nx,1:ny,2:nzm+1)
    call increment_chunk_averages(w_stor_mse,tmp,coef2)
  end if

end subroutine resetMSE

!----------------------------------------------------------------------
subroutine increment_chunk_averages(chunk_avg,field,factor)
  implicit none

  real, intent(inout) :: chunk_avg(nchunk_x,nchunk_y,nzm)
  real, intent(in) :: field(nx,ny,nzm)
  real, intent(in) :: factor

  integer :: i, ii, j, jj, k, idx_x, idx_y

  do k = 1,nzm

    idx_y = 0
    do jj = 1,nchunk_y
      do j = 1,ny_chunk
        idx_y = idx_y + 1

        idx_x = 0
        do ii = 1,nchunk_x
          do i = 1,nx_chunk
            idx_x = idx_x + 1

            chunk_avg(ii,jj,k) = chunk_avg(ii,jj,k) + field(idx_x,idx_y,k)*factor
          end do
        end do
      end do
    end do
  end do

end subroutine increment_chunk_averages

!----------------------------------------------------------------------
subroutine increment_chunk_averages2D(chunk_avg,field,factor)
  implicit none

  real, intent(inout) :: chunk_avg(nchunk_x,nchunk_y)
  real, intent(in) :: field(nx,ny)
  real, intent(in) :: factor

  integer :: i, ii, j, jj, idx_x, idx_y

  idx_y = 0
  do jj = 1,nchunk_y
    do j = 1,ny_chunk
      idx_y = idx_y + 1

      idx_x = 0
      do ii = 1,nchunk_x
        do i = 1,nx_chunk
          idx_x = idx_x + 1

          chunk_avg(ii,jj) = chunk_avg(ii,jj) + field(idx_x,idx_y)*factor
        end do
      end do
    end do
  end do
  
end subroutine increment_chunk_averages2D

!----------------------------------------------------------------------
subroutine deallocateMSE()
  ! deallocates all the MSE variables
  implicit none

  integer :: ierr

  if(.NOT.isAllocatedMSE) return

  if(do_chunked_energy_budgets) then

    deallocate( sli_before, fmse_before, qtot_before, &
         u_before, v_before, w_before, &
         u_mse, v_mse, w_mse, &
         tabs_mse, tvprime_mse, h_mse, sli_mse, &
         qtot_mse, qv_mse, qcl_mse, &
         qci_mse, qpl_mse, qpi_mse, &
         rho_mse, pp_mse, &
         cld_mse, cldcumup_mse, cldcumdn_mse, &
         STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot deallocate mean fields'
      call task_abort()
    end if

    !  fluxes to output -- computed from anomalies to domain-mean quantities
    deallocate( u_flux_mse, v_flux_mse, w_flux_mse, &
         buoy_flux_mse, h_flux_mse, sli_flux_mse , &
         qtot_flux_mse, qv_flux_mse, qcl_flux_mse, &
         qci_flux_mse, qpl_flux_mse, qpi_flux_mse, &
         STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot deallocate fluxes'
      call task_abort()
    end if

    !  vertically-resolved FMSE, SLI and QTOT budget terms for output in K/s (kg/kg/s for QTOT)
    deallocate( h_stor_mse, sli_stor_mse, qtot_stor_mse, &
         h_eddy_mse, sli_eddy_mse, qtot_eddy_mse, &
         h_sed_mse, sli_sed_mse, qtot_sed_mse, qice_sed_mse, &
         h_lsf_mse, sli_lsf_mse, qtot_lsf_mse, &
         h_misc_mse, sli_misc_mse, qtot_misc_mse , &
         qradclr_lw_mse, qrad_lw_mse, &
         qradclr_sw_mse, qrad_sw_mse, STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot deallocate arrays for 3D budgets'
      call task_abort()
    end if

    !  2D fields to output for budgets
    !  All 2D output fields have units of W/m^2 except for sst=K
    deallocate( prec_mse, prec_ice_mse, shf_mse, lhf_mse, sst_mse, &
         h_stor_2d_mse, sli_stor_2d_mse, qtot_stor_2d_mse, &
         h_eddy_2d_mse, sli_eddy_2d_mse, qtot_eddy_2d_mse, &
         h_sed_2d_mse, sli_sed_2d_mse, qtot_sed_2d_mse, qice_sed_2d_mse, &
         h_lsf_2d_mse, sli_lsf_2d_mse, qtot_lsf_2d_mse, &
         h_misc_2d_mse, sli_misc_2d_mse, qtot_misc_2d_mse, &
         qrad_lw_2d_mse, qradclr_lw_2d_mse, &
         qrad_sw_2d_mse, qradclr_sw_2d_mse, &
         lwns_mse, lwnt_mse, lwnsc_mse, lwntc_mse, &
         swns_mse, swnt_mse, swnsc_mse, swntc_mse, STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot deallocate arrays for vertically-integrated budget terms'
      call task_abort()
    end if

    if(do_chunked_momentum_budgets) then
      !  vertically-resolved momentum budget terms for output in m/s2
      deallocate( u_stor_mse, v_stor_mse, w_stor_mse, &
           u_eddy_mse, v_eddy_mse, w_eddy_mse, &
           u_pgrad_mse, v_pgrad_mse, w_pgrad_mse, &
           u_lsf_mse, v_lsf_mse, w_lsf_mse, &
           u_misc_mse, v_misc_mse, w_misc_mse, &
           w_buoy_mse, STAT=ierr)
      if(ierr.ne.0) then
        write(*,*) ' ERROR in initializeMSE within mse.f90: cannot deallocate fluxes'
        call task_abort()
      end if
    end if

    ! 1D field to output
    deallocate(pres_mse, z_mse, zi_mse, layermass_mse, STAT=ierr)
    if(ierr.ne.0) then
      write(*,*) ' ERROR in initializeMSE within mse.f90: cannot deallocate 1D quantities'
      call task_abort()
    end if

    if(do_chunk_mkbudget.and.(n_mkbudget.gt.0)) then
      deallocate(  micro_field_mse, mkbudget_vars_mse, mkbudget_before,  &
           mkbudget_stor_mse, mkbudget_eddy_mse, mkbudget_sed_mse, &
           mkbudget_lsf_mse,  mkbudget_misc_mse, mkbudget_mphy_mse, &
           mkbudget_surfflux_mse, mkbudget_extra_mse, STAT=ierr) 
      if(ierr.ne.0) then 
        write(*,*) 'ERROR when deallocating mkbudget arrays in mse.f90'
        call task_abort()
      end if
    end if

    isAllocatedMSE = .false.

  end if

end subroutine deallocateMSE


end module mse

