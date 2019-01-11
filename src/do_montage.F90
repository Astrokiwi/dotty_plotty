    subroutine do_montage(intime)
        use plot_prams, only: montage_files,anim_list_file,run_name,istart,istop,nrot,anim_type, n_pram
        implicit none
        
        character*128 :: monty_file
        character*2048 :: command
        integer :: intime

        character*5 :: width_string

        print *,"Montaging images"
        
        write(width_string,"(I5)")  n_pram
        
        write(monty_file,"('montages/montage',A,I5.5,'.png')") trim(run_name),intime
        if ( anim_type==1 ) then
            command = "montage "//trim(montage_files)//" -geometry +0+0 -tile "//trim(width_string)//"x1 "//trim(monty_file)
        else
            command = "montage "//trim(montage_files)//" -geometry +0+0 "//trim(monty_file)
        endif
        print *,trim(command)
        call system(command)
        
        if ( istop>istart .or. nrot>1  ) then
            open(unit=42,file=anim_list_file,access='append')
            write(42,*) trim(monty_file)
            close(42)
        endif
        
        montage_files = ""
    end subroutine do_montage

