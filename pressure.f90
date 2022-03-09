subroutine pressure

! call a pressure solver

use grid
use params, only: dompiensemble
implicit none

call t_startf ('pressure')

! MPI Ensemble run:
! use a modified version of pressure_orig, as done by Song Qiyu (2022)
! pressure_orig for Kuang ensemble moved to pressure_kuang by Nathanael Wong (2022)
if(dompiensemble) then
  call pressure_mpiensemble
else
  if(RUN3D) then
    if(mod(nx_gl,nsubdomains).ne.0.or.mod(ny_gl,nsubdomains).ne.0) then
      call pressure_orig
    else
      call pressure_big
    end if
  else
    call pressure_orig
  end if
end if

call t_stopf ('pressure')

end
