
subroutine damping()

!  "Spange"-layer damping at the domain top region

use vars
use microphysics, only: micro_field, index_water_vapor
use params, only: pi
implicit none

real tau_min	! minimum damping time-scale (at the top)
real tau_max    ! maxim damping time-scale (base of damping layer)
real damp_depth ! damping depth as a fraction of the domain height
parameter(tau_min=60., tau_max=1800., damp_depth=0.3)
real tau(nzm)   
integer i, j, k, n_damp

!bloss: Option for smoother onset to damping
real :: xi, tau_wb, tau_wc
integer :: kb, kc

call t_startf ('damping')

if(tau_min.lt.2*dt) then
   print*,'Error: in damping() tau_min is too small!'
   call task_abort()
end if

if(doSmoothDamping) then
  tau(:) = 0.
  do k = 1,nzm
    if(z(k).gt.zbot_SmoothDamping) then
      xi = (z(k)-zbot_SmoothDamping)/(zi(nz)-zbot_SmoothDamping)
      tau(k) = (1./tau_SmoothDamping)*sin(0.5*pi*xi)**2
      if(nstep.eq.1.AND.icycle.eq.1.AND.masterproc) then
        write(*,*) 'Smooth Damping, k, z(k), xi, 1/tau = ',k,z(k),xi, tau(k)
      end if
    end if
  end do

  n_damp=nzm-1
  do k = 1,nzm
    if(tau(k).gt.0.) then
      n_damp = nzm-k
      EXIT
    end if
  end do
  if(nstep.eq.1.AND.icycle.eq.1.AND.masterproc) then
    write(*,*) 'Smooth Damping, n_damp = ', n_damp
  end if

else
  ! standard nudging setup
  do k=nzm,1,-1
    if(z(nzm)-z(k).lt.damp_depth*z(nzm)) then 
      n_damp=nzm-k+1
    endif
  end do

  do k=nzm,nzm-n_damp,-1
    tau(k) = tau_min *(tau_max/tau_min)**((z(nzm)-z(k))/(z(nzm)-z(nzm-n_damp)))
    tau(k)=1./tau(k)
  end do
end if


do k = nzm, nzm-n_damp, -1
  kb = MAX(1,k-1)
  tau_wb = 0.5*( tau(kb) + tau(k) ) 
  do j=1,ny
    do i=1,nx
      dudt(i,j,k,na)= dudt(i,j,k,na)-(u(i,j,k)-u0(k)) * tau(k)
      dvdt(i,j,k,na)= dvdt(i,j,k,na)-(v(i,j,k)-v0(k)) * tau(k)
      dwdt(i,j,k,na)= dwdt(i,j,k,na)-w(i,j,k) * tau_wb
      !      t(i,j,k)= t(i,j,k)-dtn*(t(i,j,k)-t0(k)) * tau(k)
      !      micro_field(i,j,k,index_water_vapor)= micro_field(i,j,k,index_water_vapor)- &
      !                                    dtn*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)-q0(k)) * tau(k)
    end do! i 
  end do! j
end do ! k

if(dostatis) then
  tkedamp(:) = 0.
  do k = nzm, nzm-n_damp, -1
    kb = MAX(1,k-1)
    kc = MIN(nzm,k+1)
    tau_wb = 0.5*( tau(kb) + tau(k) ) 
    tau_wc = 0.5*( tau(k) + tau(kc) )
   do j=1,ny
      do i=1,nx
        tkedamp(k) = 0.5*( - tau(k)*(u(i,j,k)-u0(k))**2 &
                           - tau(k)*(v(i,j,k)-v0(k))**2 &
                           - 0.5*tau_wb*w(i,j,k)**2 &
                           - 0.5*tau_wc*w(i,j,k+1)**2 )
      end do! i 
    end do! j
  end do ! k
end if


call t_stopf('damping')
end subroutine damping
