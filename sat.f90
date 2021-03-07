! Saturation vapor pressure and mixing ratio. 
! Based on Flatau et.al, (JAM, 1992:1507)
!

real function esatw(t)
implicit none
real t	! temperature (K)
real a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
!	6.105851, 0.4440316, 0.1430341e-1, &
!        0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
!        0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
     	6.11239921, 0.443987641, 0.142986287e-1, &
       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
real dt
dt = max(-80.,t-273.16)
esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
end
        
        
        
real function qsatw(t,p)
implicit none
real t	! temperature (K)
real p	! pressure    (mb)
real esat,esatw
esat = esatw(t)
qsatw = 0.622 * esat/max(esat,p-esat)
end
        
        
real function dtesatw(t)
implicit none
real t	! temperature (K)
real a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
          0.443956472, 0.285976452e-1, 0.794747212e-3, &
          0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
         -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
real dt
dt = max(-80.,t-273.16)
dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
end
        
        
real function dtqsatw(t,p)
implicit none
real t	! temperature (K)
real p	! pressure    (mb)
real dtesatw
dtqsatw=0.622*dtesatw(t)/p
end
   
   
real function esati(t)
implicit none
real t	! temperature (K)
real a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	
real dt
real esatw
if(t.gt.273.15) then
  esati = esatw(t)
else if(t.gt.185.) then
  dt = t-273.16
  esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
else   ! use some additional interpolation below 184K
  dt = max(-100.,t-273.16)
  esati = 0.00763685 + dt*(0.000151069+dt*7.48215e-07)
end if
end
        
        
        
real function qsati(t,p)
implicit none
real t	! temperature (K)
real p	! pressure    (mb)
real esat,esati
esat=esati(t)
qsati=0.622 * esat/max(esat,p-esat)
end
        
        
real function dtesati(t)
implicit none
real t	! temperature (K)
real a0,a1,a2,a3,a4,a5,a6,a7,a8 
data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
	0.503223089, 0.377174432e-1,0.126710138e-2, &
	0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
	0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
real dt
real dtesatw
if(t.gt.273.15) then
  dtesati = dtesatw(t)
else if(t.gt.185.) then
  dt = t-273.16
  dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
else  ! use additional interpolation below 185K
  dt = max(-100.,t-273.16)
  dtesati = 0.0013186 + dt*(2.60269e-05+dt*1.28676e-07)
end if
end
        
        
real function dtqsati(t,p)
implicit none
real t	! temperature (K)
real p	! pressure    (mb)
real dtesati
dtqsati=0.622*dtesati(t)/p
end
      

real function get_inversion_height(nzm,z_in,p_in,tli_in,tabs_in,q_in)
implicit none
! inversion height is estimated as the height of the model
!   interface where d(RH)/dz * d(t_li)/dz takes on its
!   minimum value

integer, intent(in) :: nzm
real, intent(in) ::  z_in(nzm)	! altitude (m)
real, intent(in) ::  p_in(nzm)	! pressure    (mb)
real, intent(in) ::  tli_in(nzm)	! liquid-ice static energy/Cp (K)
real, intent(in) ::  tabs_in(nzm)	! absolute temperature (K)
real, intent(in) ::  q_in(nzm)	! total water mass mixing ratio (kg/kg)

real :: rh(nzm) ! relative humidity
real :: drh_dtli(nzm-1) ! d(rh)/dz * d(t_li)/dz
real :: drh_dtli_min ! min value of d(rh)/dz * d(t_li)/dz
integer :: kinv(1)

integer :: k
real, external :: qsatw

do k = 1,nzm
  rh(k) = q_in(k) / qsatw( tabs_in(k), p_in(k) )
end do

do k = 1,nzm-1
  drh_dtli(k) = (rh(k+1)-rh(k)) / (z_in(k+1)-z_in(k)) &
       * (tli_in(k+1) - tli_in(k)) / (z_in(k+1)-z_in(k))
end do
  
kinv = MINLOC(drh_dtli)
get_inversion_height = 0.5*(z_in(kinv(1)) + z_in(kinv(1)+1))

end function get_inversion_height
      
!-----------------------------------------
subroutine saturation_adjustment_liquid(N_in,Tl_in,qt_in,p_in,T_out,qc_out)
  use params, only: fac_cond
  !bloss/qt: Utility routine based on cloud.f90 in 
  !  MICRO_SAM1MOM that was written by Marat Khairoutdinov.
  !  This routine performs a saturation adjustment for
  !  cloud liquid water only using a Newton method.
  !  While 20 iterations are allowed, most often this
  !  routine should exit in five iterations or less.
  !  Only a single calculation of the saturation vapor
  !  pressure is required in subsaturated air.

  implicit none

  integer, intent(in) :: N_in
  real, intent(in), dimension(N_in) :: Tl_in ! liquid water temperature, K
  real, intent(in), dimension(N_in) :: qt_in ! total water mass mixing ratio, kg/kg
  real, intent(in), dimension(N_in) :: p_in ! pressure, Pa

  real, intent(out), dimension(N_in) :: T_out  ! absolute temperature, K
  real, intent(out), dimension(N_in) :: qc_out  ! cloud liquid mass mixing ratio, kg/kg
  ! NOTE: If you need qv_out, compute it as qv_out = qt_in - qc_out

  real, external :: qsatw, dtqsatw

  real tabs1, dtabs, thresh, qsat1, fff, dfff
  integer k, niter

  integer, parameter :: maxiter = 20

  !bloss/qt: quick saturation adjustment to compute cloud liquid
  !water content.
  ! NOTE: In saturation adjustment, we seek a temperature T
  !  such that Tl_in = T_out - (L/Cp)*( qt_in - qsatw(T_out,p) )
  !  This is implemented as finding the root of the function
  !    f(T_out) = T_out - (L/Cp)*( qt_in - qsatw(T_out,p) - Tl_in
  !  using Newton's method.
  do k = 1,N_in
    ! first guess is unsaturated air.
    tabs1 = Tl_in(k)
    qsat1 = qsatw( tabs1, 0.01*p_in(k) ) ! convert p_in to hPa
    
    if (qt_in(k).le.qsat1) then
      ! unsaturated, put all moisture in water vapor and account for latent heat release
      T_out(k) = tabs1
      qc_out(k) = 0.

    else
      ! if saturated, do saturation adjustment (modeled after Marat's cloud.f90).

      ! take one step with Newton's method using the above calculation of qsat
      dtabs = + fac_cond*MAX(0.,qt_in(k) - qsat1) &
           / ( 1. + fac_cond*dtqsatw( tabs1, 0.01*p_in(k) ) ) !convert p_in to hPa
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(0.01, 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        qsat1 = qsatw(tabs1, 0.01*p_in(k) ) ! saturation mixing ratio, convert pressure to hPa

        fff = Tl_in(k) - tabs1 + fac_cond*MAX(0.,qt_in(k) - qsat1)
        dfff = 1. + fac_cond*dtqsatw( tabs1, 0.01*p_in(k) ) ! convert p_in to hPa
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      T_out(k) = tabs1 ! absolute temperature for output in K
      qc_out(k) = MAX( 0.,tabs1 - Tl_in(k) )/fac_cond ! cloud liquid mass mixing ratio in kg/kg
      
      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,N_in

end subroutine saturation_adjustment_liquid

