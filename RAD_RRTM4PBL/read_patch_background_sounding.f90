SUBROUTINE read_patch_background_sounding( SoundingFileName, ptop_model, masterproc )
	
	USE netcdf
	
	USE rad, ONLY: &
		npatch_start, npatch_end, & ! indices of start/end of patched sounding
		nzsnd, & ! number of levels
		ntsnd, & ! number of time record
		psnd, &  ! pressure sounding, mb
		tsnd, &  ! time-interpolated temperature
		qsnd, &  ! time-interpolatedwater vapor
		o3snd, & ! time-interpolated ozone mass mixing ratio
		tsndng, &  ! temperature sounding, K
		qsndng, &  ! water vapor sounding, kg/kg
		o3sndng, & ! ozone mass mixing ratio, kg/kg
		have_o3mmr ! TRUE if iopfile has ozone profile
	
        USE rad_driver, ONLY: constructRestartFileName

        USE grid, only: dompi, nrestart, rank, restart_sep, nsubdomains, &
             case, case_restart, caseid, caseid_restart

	IMPLICIT NONE
	
	! Inputs
	CHARACTER(LEN=250), INTENT(IN) :: SoundingFileName
	REAL, INTENT(IN) :: ptop_model ! pressure at top of model in mb
	LOGICAL :: masterproc ! true if MPI rank==0
	
	! Local variables
	INTEGER :: nlev
	INTEGER :: status(14)
	INTEGER :: k, ierr, nn, i, count, itmp(2)
	INTEGER :: ncid, dimIDp, dimIDt, varID, ntime
	REAL, ALLOCATABLE, DIMENSION(:) :: psnd_in, rtmp
	REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: tsnd_in, qsnd_in, o3snd_in
	CHARACTER(LEN=nf90_max_name) :: tmpName
        character(LEN=256) :: filename
        integer :: idummy1, idummy2

	!++++++++++++++++++++++++ Soundings when input file is not specified +++++++++++++++++++++++++++
	! The soundings are created and tested with WRF by Cavallo et al. (2010, MWR).
	! The soundings are weighted mean pressure and temperature profiles from midlatitude summer
	! (MLS), midlatitude winter (MLW), sub-Arctic winter (SAW),sub-Arctic summer (SAS), and tropical
	! (TROP) standard atmospheres.
	INTEGER, PARAMETER :: nz60 = 60
	REAL, DIMENSION(nz60) :: &
		pprof = (/1000.00, 855.47, 731.82, 626.05, 535.57, 458.16, 391.94, 335.29, 286.83, 245.38, &
		           209.91, 179.57, 153.62, 131.41, 112.42,  96.17,  82.27,  70.38,  60.21,  51.51, &
		            44.06,  37.69,  32.25,  27.59,  23.60,  20.19,  17.27,  14.77,  12.64,  10.81, &
		             9.25,   7.91,   6.77,   5.79,   4.95,   4.24,   3.63,   3.10,   2.65,   2.27, &
		             1.94,   1.66,   1.42,   1.22,   1.04,   0.89,   0.76,   0.65,   0.56,   0.48, &
		             0.41,   0.35,   0.30,   0.26,   0.22,   0.19,   0.16,   0.14,   0.12,   0.10 /)
	REAL, DIMENSION(nz60) :: &
		tprof = (/ 286.96, 281.07, 275.16, 268.11, 260.56, 253.02, 245.62, 238.41, 231.57, 225.91, &
		           221.72, 217.79, 215.06, 212.74, 210.25, 210.16, 210.69, 212.14, 213.74, 215.37, &
		           216.82, 217.94, 219.03, 220.18, 221.37, 222.64, 224.16, 225.88, 227.63, 229.51, &
		           231.50, 233.73, 236.18, 238.78, 241.60, 244.44, 247.35, 250.33, 253.32, 256.30, &
		           259.22, 262.12, 264.80, 266.50, 267.59, 268.44, 268.69, 267.76, 266.13, 263.96, &
		           261.54, 258.93, 256.15, 253.23, 249.89, 246.67, 243.48, 240.25, 236.66, 233.86 /)
	REAL :: get_qv ! function to get qv from snd
	REAL :: get_tabs ! function to get tabs from snd
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	! Check if input file is specified or not
	IF ( TRIM(SoundingFileName) == 'NoInputFile' ) THEN
		nlev = nz60
		nzsnd = nlev + 1 ! one exlayer will be added just above the domain top
		ntsnd = 1
		ALLOCATE( psnd(nzsnd), tsnd(nzsnd), qsnd(nzsnd), o3snd(nzsnd), STAT=ierr )
		IF (ierr /= 0) THEN
			WRITE(*,*) 'ERROR in read_patch_background_sounding.f90:'
			WRITE(*,*) '***** Could not allocate arrays to store soundings'
			CALL rad_error()
		ENDIF
		psnd(:) = 0.0
		tsnd(:)  = 0.0
		qsnd(:)  = 0.0
		o3snd(:) = 0.0
		psnd(1:nlev)  = pprof(1:nlev) ! mb
		tsnd(1:nlev)  = tprof(1:nlev) ! K
		qsnd(1:nlev)  = get_qv( ptop_model ) ! kg/kg, constant profile above domain top
		o3snd(1:nlev) = 0.0 ! kg/kg
		have_o3mmr = .FALSE.
	ELSE
   ! iopfile exists, so read profiles from SCAM IOP data file
   !   do this on master processor and then broadcast results to all other nodes.
   if(masterproc) then

		status(:) = nf90_NoErr
		status(1) = nf90_open( TRIM(SoundingFileName), nf90_nowrite, ncid)

		! Get number of pressure levels
		status(2) = nf90_inq_dimid( ncid, "lev", dimIDp)
		status(3) = nf90_inquire_dimension( ncid, dimIDp, tmpName, nlev )

		! Get number of time levels
		status(4) = nf90_inq_dimid( ncid, "time", dimIDt )
		IF ( status(4) /= nf90_NoErr ) THEN
			status(4) = nf90_inq_dimid( ncid, "tsec", dimIDt )
		ENDIF
		status(5) = nf90_inquire_dimension( ncid, dimIDt, tmpName, ntime )

		ALLOCATE( psnd_in(nlev), tsnd_in(1,1,nlev,ntime), qsnd_in(1,1,nlev,ntime), &
			o3snd_in(1,1,nlev,ntime), rtmp(nlev + 3*nlev*ntime), STAT=ierr)
		IF (ierr /= 0) THEN
			WRITE(*,*) 'ERROR in read_patch_background_sounding.f90:'
			WRITE(*,*) '***** Could not allocate arrays to read in sounding'
			CALL rad_error()
		ENDIF
		
		! Get pressure levels (in Pascal)
		status(6) = nf90_inq_varid( ncid, "lev", varID )
		status(7) = nf90_get_var( ncid, varID, psnd_in )
		
		! Get temperature and moisture (K and kg/kg, respectively)
		status(8)  = nf90_inq_varid( ncid, "T", varID )
		status(9)  = nf90_get_var( ncid, varID, tsnd_in )
		status(10) = nf90_inq_varid( ncid, "q", varID )
		status(11) = nf90_get_var( ncid, varID, qsnd_in )

		status(12) = nf90_inq_varid( ncid, "o3mmr", varID )
		IF (status(12) == 0) THEN
			status(13) = nf90_get_var( ncid, varID, o3snd_in )
		ELSE
			o3snd_in(:,:,:,:) = -999.
		ENDIF

                status(14) = nf90_close( ncid )

                ! error checking on moisture sounding
                IF (MINVAL(qsnd_in).lt.0.) THEN
                  write(*,*) '============================================'
                  WRITE(*,*) 'ERROR in read_patch_background_sounding.f90:'
                  WRITE(*,*) '***** Negative water in background sounding.'
                  write(*,*) '  Stopping model .....'
                  call task_abort()
                END IF

                ! load temporary array with all data read on masterproc
                ! Reverse order from top-down to bottom-up for SAM/RRTMG.
                count = 0
                do k = 1,nlev
                  rtmp(count+1) = psnd_in(nlev+1-k)/100. ! convert from Pa to hPa
                  count = count + 1
                end do

                do nn = 1,ntime
                  do k = 1,nlev
                    rtmp(count+1) = tsnd_in(1,1,nlev+1-k,nn)
                    rtmp(count+2) = qsnd_in(1,1,nlev+1-k,nn)
                    rtmp(count+3) = o3snd_in(1,1,nlev+1-k,nn)
                    count = count + 3
                  end do
                end do

                DEALLOCATE( psnd_in, tsnd_in, qsnd_in, o3snd_in, STAT=ierr )
                IF (ierr /= 0) THEN
                  WRITE(*,*) 'ERROR in read_patch_background_sounding.f90:'
                  WRITE(*,*) '***** Could not allocate arrays to read in sounding'
                  CALL rad_error()
                ENDIF

              end if
              
              if(dompi) then
                ! if not the masterproc, wait here until data is read
                call task_barrier()
                if(masterproc) then
                  itmp(1) = nlev
                  itmp(2) = ntime
                end if

                ! share the values of nlev and ntime from masterproc
                call task_bcast_integer(0,itmp,2)
                nlev = itmp(1)
                ntime = itmp(2)
              end if

              ! allocate arrays that will hold soundings.
              nzsnd = nlev
              ntsnd = ntime ! has to be >= 2
              ALLOCATE(psnd(nzsnd), tsnd(nzsnd), qsnd(nzsnd), o3snd(nzsnd), &
                   tsndng(nzsnd,ntsnd), qsndng(nzsnd,ntsnd), o3sndng(nzsnd,ntsnd), STAT=ierr)
              IF (ierr /= 0) THEN
                WRITE(*,*) 'ERROR in read_patch_background_sounding.f90:'
                WRITE(*,*) '***** Could not allocate arrays to save sounding'
                CALL rad_error()
              ENDIF

              ! initialize soundings to zero
              psnd(:) = 0.0
              tsnd(:) = 0.0
              qsnd(:) = 0.0
              o3snd(:) = 0.0

              ! allocate rtmp array on other nodes and broadcast from master proc
              if(dompi) then
                if(.NOT.masterproc) then
                  ALLOCATE(rtmp(nlev+3*ntime*nlev), STAT=ierr)
                  IF (ierr /= 0) THEN
                    WRITE(*,*) 'ERROR in read_patch_background_sounding.f90:'
                    WRITE(*,*) '***** Could not allocate arrays to broadcast sounding'
                    CALL rad_error()
                  ENDIF
                end if

                call task_bcast_real(0,rtmp,nlev+3*ntime*nlev)
              end if

              ! fill sounding arrays with data read from netcdf
              count = 0
              do k = 1,nlev
                psnd(k) = rtmp(count+1)
                count = count + 1
              end do

              do nn = 1,ntime
                do k = 1,nlev
                  tsndng(k,nn) = rtmp(count+1)
                  qsndng(k,nn) = rtmp(count+2)
                  o3sndng(k,nn) = rtmp(count+3)
                  count = count + 3
                end do
              end do

              have_o3mmr = .FALSE.
              if(MAXVAL(o3sndng(:,:)).gt.0.) then
                have_o3mmr = .TRUE.
              end if

	ENDIF ! iopfile exists or not
	
        IF(nrestart.eq.0) THEN
          ! Find whether we need this sounding to be patched on top of model.
          npatch_start = nlev + 1
          npatch_end = nlev + 1
          DO k = 1, nlev
            IF (ptop_model.LT.100.) THEN
              ! start patch at p = 0.9*ptop_model
              IF ( psnd(k) < 0.9*ptop_model ) THEN
                ! Start patch ~10hPa above model top.
                npatch_start = k
                EXIT
              END IF
            ELSE
              ! Start patch ~10hPa above model top.
              IF ( psnd(k) < ptop_model - 10.0 ) THEN
                npatch_start = k
                EXIT
              ENDIF
            ENDIF
          ENDDO

        ELSE 
          ! nrestart == 1 or 2 bloss(2020Aug): Add patch indices to
          !  the restart file output and read them on restart.  This
          !  should prevent the size of the patch from changing on
          !  restart due to rising/falling model top pressure during
          !  the course of a simulation
          IF (restart_sep) THEN
            IF (nrestart /= 2) THEN
              filename = TRIM(constructRestartFileName(case, caseId, rank, 'bin'))
            ELSE
              filename = TRIM(constructRestartFileName(case_restart, caseId_restart, rank, 'bin'))
            ENDIF
          ELSE
            IF (nrestart /= 2) THEN
              filename = TRIM(constructRestartFileName(case, caseId, nSubdomains, 'bin'))
            ELSE
              filename = TRIM(constructRestartFileName(case_restart, caseId_restart, nSubdomains, 'bin'))
            ENDIF
          END IF

          ! read patch information from restart file
          open(56, file = trim(filename), status='unknown',form='unformatted')
          READ(56) idummy1, idummy2, npatch_start, npatch_end
          CLOSE(56)

        END IF

	IF ( npatch_start <= nlev ) npatch_end = nlev
	
	! Add one extra layer just above the domain top when the specified patched profiles are
	! used in order to avoide having artificial cooling near the domain top due to unrealistic
	! dtabs/dz between tabs at the domain top and the first patched level.
	IF ( TRIM(SoundingFileName) == 'NoInputFile' .AND. npatch_start <= nlev ) THEN
		! Shift patched soundings by 1 from the top to npatch_start
		DO k = nlev+1, npatch_start+1, -1
			psnd(k) = psnd(k-1)
			tsnd(k) = tsnd(k-1)
			qsnd(k) = qsnd(k-1)
		ENDDO
                IF (ptop_model.LT.100.) THEN
                  ! Compute tabs at 0.9*ptop_model (mb)
                  psnd(npatch_start) = 0.9*ptop_model
                ELSE
                  ! Compute tabs at ptop_model-10 (mb)
                  psnd(npatch_start) = ptop_model - 10.0
                END IF
		tsnd(npatch_start) = get_tabs( psnd(npatch_start) )
		! qsnd is not necessary since qsnd is constant
		! Adjust npatch_end
		npatch_end = nlev + 1
	ENDIF
	
        IF (  TRIM(SoundingFileName) /= 'NoInputFile' ) THEN
          tsnd(:) = tsndng(:,1)
          qsnd(:) = qsndng(:,1)
          o3snd(:) = o3sndng(:,1)
        ENDIF

	IF (masterproc) THEN
		WRITE(*,*)
		WRITE(*,*) 'Number of time background sounding:', ntsnd
		IF (have_o3mmr) THEN
			WRITE(*,*) 'Background sounding, p (mb), T (K), q (kg/kg), o3 (kg/kg)'
		ELSE
			WRITE(*,*) 'Background sounding, p (mb), T (K), q (kg/kg)'
		ENDIF
		DO k = 1, nzsnd
			IF (k == npatch_start) WRITE(*,*) '**** Start Patch Here *****'
			IF (have_o3mmr) THEN
				WRITE(*,998) psnd(k), tsnd(k), qsnd(k), o3snd(k)
			ELSE
				WRITE(*,998) psnd(k), tsnd(k), qsnd(k)
			ENDIF
			998 FORMAT(f8.3,f8.3,2e12.4)
			IF (k == npatch_end) WRITE(*,*) '**** End Patch Here *****'
		ENDDO
		IF (npatch_start > nlev) WRITE(*,*) '**** No patching required -- model top is deeper than sounding ****'
		WRITE(*,*)
	ENDIF
	
        if(dompi) call task_barrier()

END SUBROUTINE read_patch_background_sounding

!========== Functions & subroutine for the specified patched sounding ==========
! The following procedures assume:
! 1. snd file contains values at the domain top.
! 2. constant dz is used in snd file if soundings are defined with z.
! 3. if nzsnd is odd number, zsnd(1,:) = 0.0
! 4. if nzsnd is even number, zsnd(1,:) > 0.0

REAL FUNCTION get_qv( pres_in )
! Return qv at pres_in with linear interpolation (or extrapolation)
	USE vars, ONLY: qsnd
	IMPLICIT NONE
	REAL, INTENT(IN) :: pres_in
	REAL, DIMENSION(2) :: pres_out
	INTEGER, DIMENSION(2) :: index
	REAL :: coef
	! Get two pressure values and indeices closes to pres_in
	CALL pressure_index( pres_in, pres_out, index )
	! Linear interpolation: pres_out(1) > pers_out(2)
	coef = ( pres_out(1) - pres_in ) / ( pres_out(1) - pres_out(2) )
	get_qv = ( 1.0 - coef ) * qsnd(index(1),1) + coef * qsnd(index(2),1)
	! Return value has to be in kg/kg
	get_qv = get_qv * 1.0E-3 ! qsnd (g/kg) ==> get_qv (kg/kg)
	RETURN
END FUNCTION get_qv

REAL FUNCTION get_tabs( pres_in )
! Return tabs at pres_in with linear interpolation (or extrapolation)
	USE vars, ONLY: tsnd
	USE params, ONLY: cp, rgas
	IMPLICIT NONE
	REAL, INTENT(IN) :: pres_in
	REAL, DIMENSION(2) :: pres_out
	INTEGER, DIMENSION(2) :: index
	REAL :: tabs1, tabs2, coef
	! Get two pressure values and indeices closes to pres_in
	CALL pressure_index( pres_in, pres_out, index )
	! Absolute temperature at index
	tabs1 = tsnd(index(1),1) * ( pres_out(1) * 1.0E-3 ) ** ( rgas/cp )
	tabs2 = tsnd(index(2),1) * ( pres_out(2) * 1.0E-3 ) ** ( rgas/cp )
	! Linear interpolation: pres_out(1) > pers_out(2)
	coef = ( pres_out(1) - pres_in ) / ( pres_out(1) - pres_out(2) )
	get_tabs = ( 1.0 - coef ) * tabs1 + coef * tabs2
END FUNCTION get_tabs

SUBROUTINE pressure_index( pres_in, pres_out, index )
! Return two pressure levels (value and index) closest to pres_in
	USE grid, ONLY: masterproc, pres0 ! reference surface pressure, mb
	USE vars, ONLY: nzsnd, zsnd, psnd, tsnd, qsnd, nzsnd
	USE params, ONLY: cp, rgas, epsv, ggr
	IMPLICIT NONE
	REAL, INTENT(IN) :: pres_in
	REAL, DIMENSION(2), INTENT(OUT) :: pres_out
	INTEGER, DIMENSION(2), INTENT(OUT) :: index
	REAL :: dzsnd, tvsnd
	REAL, ALLOCATABLE, DIMENSION(:) :: zisnd, prsnd, pisnd, pmsnd
	INTEGER :: k, m(1)
	! Sounding defined with z or p?
	IF ( zsnd(2,1) > zsnd(1,1) ) THEN
		! Soundings are defined with z
		! Allocate arrays
		ALLOCATE( zisnd(nzsnd+1), prsnd(nzsnd+1), pisnd(nzsnd+1), pmsnd(nzsnd) )
		! Interface levels
		dzsnd = zsnd(2,1) - zsnd(1,1)
		zisnd(1) = 0.0
		IF ( MOD(nzsnd,2) /= 0 ) THEN
			zisnd(2) = 0.5 * dzsnd
			DO k = 3, nzsnd+1
				zisnd(k) = zisnd(k-1) + dzsnd
			ENDDO
		ELSE
			DO k = 2, nzsnd+1
				zisnd(k) = zisnd(k-1) + dzsnd
			ENDDO
		ENDIF
		! Compute pressure levels
		prsnd(1) = ( pres0 * 1.0E-3 ) ** ( rgas/cp )
		pisnd(1) = pres0
	!	IF ( masterproc ) PRINT*, 'k, zsnd(k), pmsnd(k)'
		DO k = 1, nzsnd
			tvsnd = tsnd(k,1) * ( 1.0 + epsv * qsnd(k,1) * 1.0E-3 )
			prsnd(k+1) = prsnd(k) - ggr / cp / tvsnd * ( zisnd(k+1) - zisnd(k) )
			pisnd(k+1) = 1.0E3 * prsnd(k+1) ** ( cp/rgas )
			pmsnd(k) = EXP( LOG( pisnd(k) ) + LOG( pisnd(k+1) / pisnd(k) ) &
					  * ( zsnd(k,1) - zisnd(k) ) / ( zisnd(k+1) - zisnd(k) ) )
		!	IF ( masterproc ) PRINT*, k, zsnd(k,1), pmsnd(k)
		ENDDO
		! Deallocate arrays except pmsnd
		DEALLOCATE( zisnd, prsnd, pisnd )
	ELSE IF ( psnd(2,1) < psnd(1,1) ) THEN
		! Soundings are defined with p
		ALLOCATE( pmsnd(nzsnd) )
		! Store psnd into pmsnd
		pmsnd(:) = psnd(:,1)
	ENDIF
	! Pressure level indices closest to pres_in: index(1) < index(2)
	m = MINLOC( ABS( pmsnd(:) - pres_in ) )
	index(1) = m(1) - 1
	index(2) = m(1)
	! Return pressure values: pres_out(1) > pres_out(2)
	pres_out(1) = pmsnd( index(1) )
	pres_out(2) = pmsnd( index(2) )
	! Deallocate pmsnd
	DEALLOCATE( pmsnd )
!	IF ( masterproc ) THEN
!		PRINT*, 'index(1), pres_out(1)', index(1), pres_out(1)
!		PRINT*, 'index(2), pres_out(2)', index(2), pres_out(2)
!	ENDIF
END SUBROUTINE pressure_index
