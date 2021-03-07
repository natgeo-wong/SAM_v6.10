subroutine InterpolateFromForcings(Nt_in,Nz_in,day_in,z_in,pres_in,field_in, &
     Nz_out,day_out,z_out,pres_out,field_out,ExtrapolateAboveSounding)
  use grid, only: masterproc
  implicit none

  integer, intent(in) :: Nt_in, Nz_in, Nz_out  
  real, intent(in) :: day_in(Nt_in), z_in(Nz_in,Nt_in), pres_in(Nz_in,Nt_in), field_in(Nz_in,Nt_in)
  real, intent(in) :: day_out, z_out(Nz_out), pres_out(Nz_out)
  logical, intent(in) :: ExtrapolateAboveSounding
  real, intent(out) :: field_out(Nz_out)

  integer :: i,iz,k,m,n,nn
  logical :: zgrid
  real :: coef
  real :: field_tmp(Nz_out,2)
  real :: missing_value = -9999.

  ! find index of day_in(:) so that day_in(nn) < day_out <= day_in(nn+1)
  nn=1
  do i=1,Nt_in-1
    if(day_out.gt.day_in(i)) then
      nn=i
    endif
  end do

  ! extract field_in(:,nn:nn+1)
  do n=1,2
  
    m = nn+n-1
    ! figure out whether to use height or pressure as vertical coordinate
    if(z_in(2,m).gt.z_in(1,m)) then
      zgrid=.true.
    else if(pres_in(2,m).lt.pres_in(1,m)) then
      zgrid=.false.
    else
      if(masterproc) print*,'error in grid in snd'
      call task_abort()
    end if
    do iz = 1,Nz_out
      if(zgrid) then
        do i = 2,Nz_in
          if(z_out(iz).le.z_in(i,m)) then
           coef = (z_out(iz)-z_in(i-1,m))/(z_in(i,m)-z_in(i-1,m)) 
           field_tmp(iz,n)=field_in(i-1,m)+(field_in(i,m)-field_in(i-1,m))*coef
           goto 11
          endif
        end do
      else
        do i = 2,Nz_in
          if(pres_out(iz).ge.pres_in(i,m)) then
           coef = (pres_out(iz)-pres_in(i-1,m))/(pres_in(i,m)-pres_in(i-1,m))
           field_tmp(iz,n)=field_in(i-1,m)+(field_in(i,m)-field_in(i-1,m))*coef
           goto 11
          endif
        end do
      end if

      if(ExtrapolateAboveSounding) then
        ! constant value above input sounding
        field_tmp(iz,n)=field_tmp(iz-1,n)
      else
        field_tmp(iz,n)=missing_value
      end if

   11 continue

    end do ! iz 
  end do ! n

  coef=(day_out-day_in(nn))/(day_in(nn+1)-day_in(nn))
  do k=1,Nz_out
    field_out(k)=field_tmp(k,1)+(field_tmp(k,2)-field_tmp(k,1))*coef
  end do

end subroutine InterpolateFromForcings
