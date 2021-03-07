module aerosol_utils
  !bloss(2020-11): Added a function to compute aerosol mass fraction.

  implicit none

  private
  public :: DryAerosolMassFraction, PartitionAerosolMass

  integer, parameter :: rb = selected_real_kind(12) ! 8 byte real

contains
  
  subroutine PartitionAerosolMass( N, nacc, qacc, SigmaG, ncl, nad, qad, qaw)
    integer, intent(in) :: N
    real, dimension(N), intent(in) :: nacc ! total aerosol number
    real, dimension(N), intent(in) :: qacc ! total aerosol mass
    real, intent(in) :: SigmaG ! geometric standard deviation of lognormal aerosol size distn.
    real, dimension(N), intent(inout) :: ncl  ! cloud droplet number
    real, dimension(N), intent(out) :: nad  ! dry aerosol number
    real, dimension(N), intent(out) :: qad  ! dry aerosol mass
    real, dimension(N), intent(out) :: qaw  ! wet (in-cloud-droplet) aerosol mass

    integer :: k

    !bloss(2020-11): Partition total (dry+wet) aerosol mass into
    !   dry aerosol mass and wet (in-cloud-droplet) aerosol mass
    !   based on the fraction of total aerosol number associated
    !   with cloud droplets.  Here, we assume that each cloud droplet
    !   contains a single aerosol particle, and that the total aerosol
    !   size distribution is divided, with the larger aerosols in
    !   cloud droplets and the smaller ones being dry.
    
    call t_startf('AerosolPartition')

    ! Use total aerosol to populate dry and wet aerosol categories
    do k = 1,N
      ! limit NC to be no larger than the total aerosol (NAcc=NAd+NC)
      ncl(k) = MAX(0., MIN(ncl(k),nacc(k)) )

      ! compute dry aerosol number as the NAcc - NC
      nad(k) = MAX(0., nacc(k) - ncl(k) )
      qad(k) = qacc(k) ! assume no cloud at first pass
      qaw(k) = 0.

      if(nad(k).lt.nacc(k) ) then
        qad(k) = qacc(k) &
             *DryAerosolMassFraction( nad(k)/nacc(k) , &
             SigmaG , 1 )
        qaw(k) = MAX(0., qacc(k) - qad(k) )
      end if
    end do

    call t_stopf('AerosolPartition')

  end subroutine PartitionAerosolMass

  real function DryAerosolMassFraction( DryAerosolNumberFraction, SigmaG, itag )

    real, intent(in) :: DryAerosolNumberFraction, SigmaG
    integer, intent(in) :: itag ! records which routine called this function for error message

    real(kind=rb) :: tmp1, xx
    character(LEN=80) :: message

    real(kind=rb), external :: derfi ! Netlib function

    ! This function assumes that we have an aerosol size distribution
    ! that is lognormal with a given value of its geometric standard
    ! deviation, SigmaG.  Within a cloud, part of the underlying
    ! aerosol population will be activated and serve as the nuclei of
    ! cloud droplets.  We assume that the aerosol population is
    ! divided by size, with the larger aerosols serving as cloud
    ! nuclei and the smaller aerosols being dry, interstitial
    ! aerosols.  With this assumption, we can use partial moments of
    ! the lognormal size distribution to relate the number fraction
    ! of dry aerosols (Ndry/(Ndry+NC) input here) to the mass
    ! fraction (Mdry/(Mdry+Mwet), output here).

    ! The partial moments of the lognormal give:
    !   Ndry/(Ndry+Nwet) = Phi(xx)
    ! where xx = log(D'/Dg)/log(SigmaG) in which D' is the diameter
    ! dividing dry and wet aerosol, Dg is a characteristic diameter
    ! of the lognormal size distribution.  The function Phi can be
    ! computed using the error function:
    !   Phi(xx) = 1/2 + 1/2*erf( xx/sqrt(2) )
    ! Since we have Ndry/(Ndry+Nwet) as an input, we can use the
    ! inverse error function (taken from netlib in derfi.f) to
    ! compute xx below.  Then, the mass fraction is computed as:
    !   Mdry/(Mdry+Mwet) = Phi(-3*log(SigmaG) + xx)


    if((DryAerosolNumberFraction.lt.0.) &
         .OR.(DryAerosolNumberFraction.gt.1.)) then

      select case (itag)
      case (1)
        message = TRIM('ParitionAerosolMass')
      case (2)
        message = TRIM('micro_init')
      case (3)
        message = TRIM('First Activate')
      case (4)
        message = TRIM('Second Activate')
      case (5)
        message = TRIM('Third Activate')
      end select

      write(*,*) '******************************************************'
      write(*,*) 'Warning: Out of bounds DryAerosolNumberFraction = ', DryAerosolNumberFraction
      write(*,*) 'Input from ', TRIM(message)
      write(*,*) '******************************************************'
      if((DryAerosolNumberFraction.lt.-0.01) &
           .OR.(DryAerosolNumberFraction.gt.1.01)) then
        write(*,*) 'Error greater than 1%.  Model stopping...'
        call task_abort()
      end if
    end if

!bloss    call t_startf('AerosolPartition')

    tmp1 = 2.*( DryAerosolNumberFraction - 0.5 )

    tmp1 = MAX(-1.+EPSILON(1.), MIN( 1.-EPSILON(1.), tmp1 ) )

    xx = sqrt(2.) * derfi( tmp1 )

    DryAerosolMassFraction = 0.5 + 0.5* erf( (-3.*log(SigmaG) + xx ) / sqrt(2.) )

!bloss    call t_stopf('AerosolPartition')

  end function DryAerosolMassFraction

end module aerosol_utils
