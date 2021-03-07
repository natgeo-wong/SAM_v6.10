	subroutine write_rad()
	
	use rad
        use radae, only: abstot_3d, absnxt_3d, emstot_3d
	implicit none
	character *4 rankchar
	character *256 filename, filename_save
	integer irank, index
        integer lenstr
        external lenstr

        if(masterproc) print*,'Writting radiation restart file...'

        !bloss(2020-08-11): Build filename up here to make code easier to read
        if(restart_sep) then
          index = rank
        else
          index = nsubdomains
        end if
        filename= constructRestartFileName( case, caseid, index, 'bin' )

        !bloss(2020Aug): copy the last restart file to *.old just in case
        !  something happens later.
        filename_save= constructRestartFileName( case, caseid, index, 'old' )
        
        !bloss: Now back to the old routine
        if(restart_sep) then

          !bloss(2020Aug) remove old restart file (if it exists)
          call system('rm -f ' // TRIM(filename_save) )

          !bloss(2020Aug) next, move last restart file to *.old
          call system('mv '//TRIM(filename)//' '//TRIM(filename_save) )

          open(56,file=TRIM(filename),status='unknown',form='unformatted')
          write(56) nsubdomains
	  write(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad, &
                    cld_rad, rel_rad, rei_rad, &
	      qrad,absnxt_3d,abstot_3d,emstot_3d 
          close(56)

        else

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               if(masterproc) then

                 !bloss(2020Aug) remove old restart file (if it exists)
                 call system('rm -f ' // TRIM(filename_save) )

                 !bloss(2020Aug) next, move last restart file to *.old
                 call system('mv '//TRIM(filename)//' '//TRIM(filename_save) )

                  open(56,file=TRIM(filename),status='unknown',form='unformatted')
                  write(56) nsubdomains

               else

                  open(56,file=TRIM(filename),status='unknown',form='unformatted', position='append')

               end if

	       write(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad, &
                         cld_rad, rel_rad, rei_rad, &
	      	   qrad,absnxt_3d,abstot_3d,emstot_3d 
               close(56)
           end if
        end do

        end if ! restart_sep

	if(masterproc) then
           print *,'Saved radiation restart file. nstep=',nstep
	endif

        call task_barrier()

        return
        end
 
 
 
 
     
	subroutine read_rad()
	
	use rad
        use radae, only: abstot_3d, absnxt_3d, emstot_3d
	implicit none
	character *4 rankchar
	character *256 filename
	integer irank,ii, index
        integer lenstr
        external lenstr
	
        if(masterproc) print*,'Reading radiation restart file...'

        !bloss(2020-08-11): Build filename up here to make code easier to read
        if(restart_sep) then
          index = rank
        else
          index = nsubdomains
        end if
        if(nrestart.ne.2) then
          filename= constructRestartFileName( case, caseid, index, 'bin' )
        else
          filename= constructRestartFileName( case_restart, caseid_restart, index, 'bin' )
        end if
        
        if(restart_sep) then

          open(56,file=TRIM(filename),status='unknown',form='unformatted')
          read (56)
	  read(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad, &
                         cld_rad, rel_rad, rei_rad, &
	     qrad,absnxt_3d,abstot_3d,emstot_3d 
          close(56)

        else

          do irank=0,nsubdomains-1

             call task_barrier()

             if(irank.eq.rank) then

               open(56,file=TRIM(filename),status='unknown',form='unformatted')
               read (56)

               do ii=0,irank-1 ! skip records
                 read(56)
               end do

	       read(56) initrad,nradsteps,tabs_rad,qc_rad,qi_rad,qv_rad, &
                         cld_rad, rel_rad, rei_rad, &
	  	 qrad,absnxt_3d,abstot_3d,emstot_3d 
               close(56)
             end if

          end do

        end if ! restart_sep

        if(rank.eq.nsubdomains-1) then
             print *,'Case:',caseid
             print *,'Restart radiation at step:',nstep
             print *,'Time:',nstep*dt
        endif

        call task_barrier()


        return
        end

