    subroutine do_strips
        use plot_prams
        implicit none
        
        character*128 :: strip_file
        character*2048 :: command
        character*2 :: width

        strip_file = trim(anim_file_base)//"_ensemble.png"

        write(width,"(I2.2)") n_runs

        print *,"Combining strips in ",strip_file
        command = "montage -tile "//width//"x1 "//trim(anim_files)//" -geometry +0+0 "//trim(strip_file)
        call system(command)
        
    end subroutine do_strips