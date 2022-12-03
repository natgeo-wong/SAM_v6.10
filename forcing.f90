
subroutine forcing

use vars
use params
use microphysics, only: micro_field, index_water_vapor, total_water, is_water_vapor, nmicro_fields, mk0, mklsadv, mk_ref
use simple_ocean, only: sst_evolve

implicit none


integer i,j,k,n,nn,m,iz,iday0,iday
real coef, radtend, dayy
real tt(nzm,2),qq(nzm,2),uu(nzm,2),vv(nzm,2),ww(nzm,2),tp(nzm,2),pp(nzm,2)
real tpm(nzm)
real ratio1, ratio2, ratio_t1, ratio_t2
logical zgrid, pgrid

real Derbyshire_qv(nzm), Derbyshire_theta(nzm)
real hadvs(nzm), hadvq(nzm), tmp_hadvq(nzm)
real tmp_qsat, wgt_High

integer :: ktrop
real, external :: erfff

call t_startf ('forcing')


! if doseasons=.false. do perpetual forcing

 if(doseasons) then
   dayy = day
 else
   iday0 = day0
   iday = day
   dayy = day-iday
   dayy = iday0 + dayy
 end if


! ---------------------------------------------------------------
! Large-scale sounding:

  nn=1
  do i=1,nsnd-1
   if(day.gt.daysnd(i)) then
     nn=i
   endif
  end do

  do n=1,2

    m = nn+n-1
    if(zsnd(2,m).gt.zsnd(1,m)) then
      zgrid=.true.
    else
      zgrid=.false.
    endif

    if(psnd(2,m).lt.psnd(1,m)) then
      pgrid=.true.
    else
      pgrid=.false.
    end if

    if((.not.zgrid).and.(.not.pgrid)) then
      if(masterproc) print*,'error in grid in snd'
      stop
    end if

    do iz = 1,nzm
      if(zgrid) then
        do i = 2,nzsnd
          if(z(iz).le.zsnd(i,m)) then
            coef = (z(iz)-zsnd(i-1,m))/(zsnd(i,m)-zsnd(i-1,m))
            tt(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            if(pgrid) then
              pp(iz,n)=psnd(i-1,m)+(psnd(i,m)-psnd(i-1,m))*coef
              tt(iz,n)=tt(iz,n)/((1000./pp(iz,n))**(rgas/cp))
            else
              tt(iz,n)=tt(iz,n)/prespoti(iz)
            endif
            tp(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            qq(iz,n)=qsnd(i-1,m)+(qsnd(i,m)-qsnd(i-1,m))*coef
            uu(iz,n)=usnd(i-1,m)+(usnd(i,m)-usnd(i-1,m))*coef
            vv(iz,n)=vsnd(i-1,m)+(vsnd(i,m)-vsnd(i-1,m))*coef
           goto 11
          endif
        end do
      else
        do i = 2,nzsnd
          if(pres(iz).ge.psnd(i,m)) then
            coef = (pres(iz)-psnd(i-1,m))/(psnd(i,m)-psnd(i-1,m))
            tt(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            tt(iz,n)=tt(iz,n)/prespoti(iz)
            tp(iz,n)=tsnd(i-1,m)+(tsnd(i,m)-tsnd(i-1,m))*coef
            qq(iz,n)=qsnd(i-1,m)+(qsnd(i,m)-qsnd(i-1,m))*coef
            uu(iz,n)=usnd(i-1,m)+(usnd(i,m)-usnd(i-1,m))*coef
            vv(iz,n)=vsnd(i-1,m)+(vsnd(i,m)-vsnd(i-1,m))*coef
            pp(iz,n)=psnd(i-1,m)+(psnd(i,m)-psnd(i-1,m))*coef
           goto 11
          endif
        end do
      end if

      call atmosphere(z(iz-1)/1000.,ratio1,ratio2,ratio_t1)
      call atmosphere(z(iz)/1000.,ratio1,ratio2,ratio_t2)

      tt(iz,n)=ratio_t2/ratio_t1*tt(iz-1,n)
!      qq(iz,n)=max(0.,2.*qq(iz-1,n)-qq(iz-2,n))
      qq(iz,n) = qq(iz-1,n)*exp(-(z(iz)-z(iz-1))/3000.)
      uu(iz,n)=uu(iz-1,n)
      vv(iz,n)=vv(iz-1,n)

   11 continue

    end do ! iz

  end do ! n

  coef=(day-daysnd(nn))/(daysnd(nn+1)-daysnd(nn))
  do k=1,nzm
    tg0(k)=tt(k,1)+(tt(k,2)-tt(k,1))*coef
    tp0(k)=tp(k,1)+(tp(k,2)-tp(k,1))*coef
    qg0(k)=qq(k,1)+(qq(k,2)-qq(k,1))*coef
    qg0(k)=qg0(k)*1.e-3
! Note that ug0 and vg0 maybe reset if dolargescale is true)
    ug0(k)=uu(k,1)+(uu(k,2)-uu(k,1))*coef - ug
    vg0(k)=vv(k,1)+(vv(k,2)-vv(k,1))*coef - vg
    pg0(k)=pp(k,1)+(pp(k,2)-pp(k,1))*coef - vg
   end do

   !bloss: Save sounding winds in case we want to nudge to them
   usounding0(:) = ug0(:)
   vsounding0(:) = vg0(:)

! ---------------------------------------------------------------
! Initialize tendencies:


do k=1,nzm
 ttend(k)=0.
 qtend(k)=0.
end do


! Large-Scale Advection Forcing:


if(dolargescale.and.time.gt.timelargescale) then

  nn=1
  do i=1,nlsf-1
    if(day.gt.dayls(i)) then
      nn=i
    endif
  end do

  do n=1,2

    m = nn+n-1

    if(zls(2,m).gt.zls(1,m)) then
      zgrid=.true.
    else if(pls(2,m).lt.pls(1,m)) then
      zgrid=.false.
    else
      if(masterproc) print*,'error in grid in lsf'
      stop
    end if
    do iz = 1,nzm
      if(zgrid) then
        do i = 2,nzlsf
        if(z(iz).le.zls(i,m)) then
          coef = (z(iz)-zls(i-1,m))/(zls(i,m)-zls(i-1,m))
          tt(iz,n)=dtls(i-1,m)+(dtls(i,m)-dtls(i-1,m))*coef
          qq(iz,n)=dqls(i-1,m)+(dqls(i,m)-dqls(i-1,m))*coef
          uu(iz,n)=ugls(i-1,m)+(ugls(i,m)-ugls(i-1,m))*coef
          vv(iz,n)=vgls(i-1,m)+(vgls(i,m)-vgls(i-1,m))*coef
          ww(iz,n)=wgls(i-1,m)+(wgls(i,m)-wgls(i-1,m))*coef
          goto 12
        endif
        end do
      else
        do i = 2,nzlsf
        if(pres(iz).ge.pls(i,m)) then
          coef = (pres(iz)-pls(i-1,m))/(pls(i,m)-pls(i-1,m))
          tt(iz,n)=dtls(i-1,m)+(dtls(i,m)-dtls(i-1,m))*coef
          qq(iz,n)=dqls(i-1,m)+(dqls(i,m)-dqls(i-1,m))*coef
          uu(iz,n)=ugls(i-1,m)+(ugls(i,m)-ugls(i-1,m))*coef
          vv(iz,n)=vgls(i-1,m)+(vgls(i,m)-vgls(i-1,m))*coef
          ww(iz,n)=wgls(i-1,m)+(wgls(i,m)-wgls(i-1,m))*coef
          goto 12
        endif
        end do
      end if
      tt(iz,n)=0.
      qq(iz,n)=0.
      uu(iz,n)=uu(iz-1,n)
      vv(iz,n)=vv(iz-1,n)
      ww(iz,n)=0.
  12 continue

    end do

  end do ! n

  coef=(day-dayls(nn))/(dayls(nn+1)-dayls(nn))
  dosubsidence = .false.
    do k=1,nzm
    ttend(k)=tt(k,1)+(tt(k,2)-tt(k,1))*coef
    qtend(k)=qq(k,1)+(qq(k,2)-qq(k,1))*coef
    ug0(k)=uu(k,1)+(uu(k,2)-uu(k,1))*coef - ug
    vg0(k)=vv(k,1)+(vv(k,2)-vv(k,1))*coef - vg
    wsub(k)=ww(k,1)+(ww(k,2)-ww(k,1))*coef
    dosubsidence = dosubsidence .or. wsub(k).ne.0.
    do j=1,ny
      do i=1,nx
      t(i,j,k)=t(i,j,k)+ttend(k) * dtn
      micro_field(i,j,k,index_water_vapor) = &
                  max(0.,micro_field(i,j,k,index_water_vapor) + qtend(k) * dtn)
      end do
    end do
  end do

  pres0 = pres0ls(nn)+(pres0ls(nn+1)-pres0ls(nn))*coef

  if(wgls_holds_omega) then
    ! convert omega (sitting in wsub) into large-scale vertical velocity.
    ! Note that omega was read in from SCAM IOP netcdf input file.
    do k = 1,nzm
      wsub(k) = -wsub(k)/rho(k)/ggr
    end do
  end if

  if((dofixdivg).AND.(do_linear_subsidence)) then
    !build subsidence velocity with fixed divergence.
    do k = 1,nzm
      wsub(k) = -divg_ls*z(k)
    end do
    dosubsidence = .true.
  end if

   !-------------------------------------------------------------------------------
   !UW additions
   ! save reference copy of large-scale vertical velocity before it's
   ! modified by WTG or scaling techniques
  wsub_ref(1:nzm) = wsub(1:nzm)

  if(dodgw) then

    if(am_tscale.gt.0) then
      twtgmax = (nstop * dt - timelargescale) * wtgscale_time
      twtg = time-timelargescale
      if(twtg.gt.twtgmax) then
        am_wtg_time = am_wtg
      else
        am_wtg_time = am_wtg * twtgmax / twtg
      endif
    else
      am_wtg_time = am_wtg
    endif

    if (dowtg_raymondzeng_QJRMS2005) call wtg_james2009(nzm, &
        100.*pres, tg0, qg0, tabs0, qv0, qn0+qp0, &
        fcor, lambda_wtg, am_wtg_time, am_wtg_exp, o_wtg, ktrop)
    if (dowtg_decompdgw) then
      call wtg_james2009(nzm, &
        100.*pres, tg0, qg0, tabs0, qv0, qn0+qp0, &
        fcor, lambda_wtg, am_wtg_time, am_wtg_exp, owtgr, ktrop)
      call wtg_decompdgw(masterproc, &
        nzm, nz, z, 100.*pg0, tg0, qg0, tabs0, qv0, qn0+qp0, &
        lambda_wtg, am_wtg_time, wtgscale_vertmodenum, wtgscale_vertmode, &
        o_wtg, wwtgc, ktrop)
    end if

    ! convert from omega in Pa/s to wsub in m/s
    w_wtg(1:nzm) = -o_wtg(1:nzm)/rho(1:nzm)/ggr
    if (dowtg_decompdgw) wwtgr(1:nzm) = -owtgr(1:nzm)/rho(1:nzm)/ggr

    ! add to reference large-scale vertical velocity.
    wsub(1:nzm) = wsub(1:nzm) + w_wtg(1:nzm)

    dosubsidence = .true.

    if(dowtg_tnudge.and.ktrop.lt.nzm-5) then
      ! nudge temperature to observed sounding,
      !   ramping up the nudging starting a couple of
      !   kilometers below the tropopause.
      do k = 1,nzm
        tnudge(k) = -(t0(k)+t00-tg0(k)-gamaz(k)) * itau_wtg_tnudge/86400. &
            *0.5*(1. + erfff( (z(k) - z(ktrop) + taulz_wtg_tnudge)/taulz_wtg_tnudge ) )
        do j=1,ny
          do i=1,nx
            t(i,j,k)=t(i,j,k)+tnudge(k) * dtn
          end do
        end do
      end do
    end if

    if(dowtg_qnudge) then
      ! nudge moisture profile using tau_wtg_qnudge timescale
      !   convert inverse nudging timescale from 1/d to 1/s
      ! use error function profile, so that inverse nudging
      !   timescale increases smoothly from zero at the surface
      !   to a uniform value in the free troposphere.
      do k = 1,nzm
        qnudge(k) = (qg0(k) - q0(k))*itau_wtg_qnudge/86400. &
            *0.5*(1. + erfff( (z(k) - tauz0_wtg_qnudge)/taulz_wtg_qnudge ) )
        do j=1,ny
          do i=1,nx
            micro_field(i,j,k,index_water_vapor) = &
                micro_field(i,j,k,index_water_vapor) + qnudge(k) * dtn
          end do
        end do
      end do
    end if
  end if

  if (dotgr) then

    if(tau_tscale.gt.0) then
      twtgmax = (nstop * dt - timelargescale) * wtgscale_time
      twtg = time-timelargescale
      if(twtg.gt.twtgmax) then
        tau_wtg_time = tau_wtg
      else
        tau_wtg_time = tau_wtg * twtg / twtgmax
      endif
    else
      tau_wtg_time = tau_wtg
    endif

    do k = 1,nzm
      tpm(k) = tabs0(k) * prespot(k)
    end do

    if (dowtg_raymondzeng_QJRMS2005) call wtg_qjrms2005(masterproc, nzm, nz, z, &
                        tp0, tpm, tabs0, tau_wtg_time, dowtgLBL, boundstatic, &
                        dthetadz_min, w_wtg, wwtgr)
    if (dowtg_daleuetal_JAMES2015)   call wtg_james2015(masterproc, nzm, nz, z, &
                        tp0, tpm, tabs0, tau_wtg_time, boundstatic, dthetadz_min, w_wtg)
    if (dowtg_decomptgr)            call wtg_decomptgr(masterproc, nzm, nz, z, &
                        tp0, tpm, tabs0, tau_wtg_time, &
                        wtgscale_vertmode1, wtgscale_vertmode2, &
                        dowtgLBL, boundstatic, dthetadz_min, w_wtg, wwtgr, wwtgc)

    wsub(1:nzm) = wsub(1:nzm) + w_wtg(1:nzm)

    dosubsidence = .true.

  end if

  mklsadv(:,:) = 0. ! move initialization here

  if(doDerbyshire) then

    ! initialize tendencies
    hadvq(1:nzm) = 0.
    hadvs(1:nzm) = 0.

    do k = 1,nzm
      Derbyshire_theta(k) = Derbyshire_theta0 + z(k)*Derbyshire_LapseRate

      tmp_qsat = MIN( qsatw(tabs0(k),pres(k)), qsati(tabs0(k),pres(k)))
      if(z(k).lt.Derbyshire_z2-250.) then
        ! below Derbyshire_z2 ~ 2km
        Derbyshire_qv(k) = Derbyshire_RelH_Low*tmp_qsat
      elseif(z(k).gt.Derbyshire_z2+250.) then
        ! above Derbyshire_z2 ~ 2km, below Derbyshire_z3 ~ 15km
        Derbyshire_qv(k) = Derbyshire_RelH_High*tmp_qsat
      else
        ! around Derbyshire_z2 ~ 2km, blend RH over 500m
        !   with half a sine wave
        wgt_High = 0.5*( 1. + sin(pi*(z(k)-Derbyshire_z2)/500.) )
        Derbyshire_qv(k) = tmp_qsat &
            *( Derbyshire_RelH_Low*(1.-wgt_High) &
            + Derbyshire_RelH_High*wgt_High )
      end if
    end do

    do k = 1,nzm
      if((z(k).gt.Derbyshire_z1).AND.(z(k).lt.Derbyshire_z3)) then
        ! compute nudging tendency for humidity and temperature
        hadvq(k) = (Derbyshire_qv(k) - qv0(k)) / Derbyshire_tau

        ! compute theta tendency and divide by prespot to convert back
        !   to a DSE tendency
        hadvs(k) = (Derbyshire_theta(k) - tabs0(k)*prespot(k)) / prespot(k) / Derbyshire_tau
      end if
    end do

    do k = 1,nzm
      qtend(k) = qtend(k) + hadvq(k)
      ttend(k) = ttend(k) + hadvs(k)

      do j = 1,ny
        do i = 1,nx
          t(i,j,k)=t(i,j,k)+hadvs(k) * dtn
          micro_field(i,j,k,index_water_vapor) = &
              MAX(0., micro_field(i,j,k,index_water_vapor) + hadvq(k) * dtn)
        end do
      end do
    end do

    ! now handle tendency for non-standard forms of water vapor (e.g. heavy isotopes)
    do n = 1,nmicro_fields
      if( is_water_vapor(n) .AND. (n.ne.index_water_vapor) ) then
!!$         tmp_hadvq(:) = 0.
!!$           if(hadvq(k).le.0.) then
!!$             ! if tendency is drying, remove moisture with mean isotopic composition of domain
!!$             tmp_hadvq(k) = hadvq(k)*mk0(k,n)/mk0(k,index_water_vapor)
!!$           else
!!$             ! if tendency is moistening, add moisture with reference isotopic composition
!!$             tmp_hadvq(k) = hadvq(k)*mk_ref(k,n)/mk_ref(k,index_water_vapor)
!!$           end if
        tmp_hadvq(:) = hadvq(:)*mk0(1:nzm,n)/mk0(1:nzm,index_water_vapor)
        mklsadv(1:nzm,n) = mklsadv(1:nzm,n) + tmp_hadvq(1:nzm)*float(nx*ny)
        do k = 1,nzm
          do j = 1,ny
            do i = 1,nx
              micro_field(i,j,k,n)= MAX(0., micro_field(i,j,k,n) + dtn*tmp_hadvq(k) )
            end do
          end do
        end do
      end if
    end do

  end if

  if(dosubsidence) call subsidence()

end if

!---------------------------------------------------------------------
! Prescribed Radiation Forcing:


if(doradforcing.and.time.gt.timelargescale) then

  nn=1
  do i=1,nrfc-1
   if(day.gt.dayrfc(i)) then
     nn=i
   endif
  end do

  do n=1,2

    m = nn+n-1

    if(prfc(2,m).gt.prfc(1,m)) then
     zgrid=.true.
    else if(prfc(2,m).lt.prfc(1,m)) then
     zgrid=.false.
    else
     if(masterproc) print*,'error in grid in rad'
     stop
    end if
    do iz = 1,nzm
     if(zgrid) then
      do i = 2,nzrfc
       if(z(iz).le.prfc(i,m)) then
        tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                           *(z(iz)-prfc(i-1,m))
        goto 13
       endif
      end do
     else
      do i = 2,nzrfc
       if(pres(iz).ge.prfc(i,m)) then
        tt(iz,n)=dtrfc(i-1,m)+(dtrfc(i,m)-dtrfc(i-1,m))/(prfc(i,m)-prfc(i-1,m)) &
                                                           *(pres(iz)-prfc(i-1,m))
        goto 13
       endif
      end do
     end if
     tt(iz,n)=0.
  13 continue
    end do

  end do ! n

  coef=(day-dayrfc(nn))/(dayrfc(nn+1)-dayrfc(nn))
  do k=1,nzm
    radtend=tt(k,1)+(tt(k,2)-tt(k,1))*coef
    radqrlw(k)=radtend*float(nx*ny)
    radqrsw(k)=0.
    do j=1,ny
     do i=1,nx
       t(i,j,k)=t(i,j,k)+radtend*dtn
     end do
    end do
  end do

endif


!----------------------------------------------------------------------------
! Surface flux forcing:

if(dosfcforcing.and.time.gt.timelargescale) then

  nn=1
  do i=1,nsfc-1
    if(day.gt.daysfc(i)) then
      nn=i
    endif
  end do

  coef=(day-daysfc(nn))/(daysfc(nn+1)-daysfc(nn))
  tabs_s=sstsfc(nn)+(sstsfc(nn+1)-sstsfc(nn))*coef
  fluxt0=(shsfc(nn)+(shsfc(nn+1)-shsfc(nn))*coef)/(rhow(1)*cp)
  fluxq0=(lhsfc(nn)+(lhsfc(nn+1)-lhsfc(nn))*coef)/(rhow(1)*lcond)
  tau0=tausfc(nn)+(tausfc(nn+1)-tausfc(nn))*coef

  do j=1,ny
    do i=1,nx
      sstxy(i,j) = tabs_s - t00
    end do
  end do

  if(dostatis) then
    sstobs = tabs_s  ! sst is not averaged over the sampling period
    lhobs = lhobs + fluxq0 * rhow(1)*lcond
    shobs = shobs + fluxt0 * rhow(1)*cp
  end if

endif

!----------------------------------------------------------------------------
! Temperature Tendency Forcing:
!	Simple Radiative Tendencies taken from Pauluis & Garner [2006]

if(doradtendency.and.time.gt.timelargescale) then

  do k = 1,nzm
    do j=1,ny
      do i=1,nx
        if (tabs(i,j,k)>207.5) then
            t(i,j,k) = t(i,j,k) - dtn * troptend / 86400
        else
            t(i,j,k) = t(i,j,k) + dtn * (200 - tabs(i,j,k)) / (5*86400)
        end if
      end do
    end do
  end do

endif

!-------------------------------------------------------------------------------

if(.not.dosfcforcing.and.dodynamicocean) then

  call sst_evolve()

endif

!-------------------------------------------------------------------------------
call t_stopf ('forcing')


end subroutine forcing
