    subroutine do_anim
        use plot_prams
        implicit none
        
        character*128 :: anim_file,mpeg_file
        character*2048 :: command
        
!        if ( n_pram==1 ) then        
!            write(animy_file,"('anims/anim',A,A,'.gif')") trim(run_name),trim(all_prams(1)%plot_name)
!        else
!            write(animy_file,"('anims/anim',A,'.gif')") trim(run_name)
!        endif

        if ( anim_type==0 ) then        
            anim_file = trim(anim_file_base)//trim(run_name)//".gif"

            print *,"Producing animation in ",anim_file
            !command = "convert @"//trim(anim_list_file)//" -loop 0 "//trim(anim_file)
            command = "convert -delay 10 @"//trim(anim_list_file)//" -loop 0 "//trim(anim_file)
            ! delay 5 for production
            call system(command)
            
            print *,"Converting to mpeg (if ffmpeg installed) in",mpeg_file
            mpeg_file = trim(anim_file_base)//trim(run_name)//".mp4"
            command="ffmpeg -y -r 24 -i "//anim_file//" -c:v mpeg4 -q:v 1 "//mpeg_file
            call system(command)
            print *,"Tidying up and deleting ",anim_file
            command = "rm "//anim_file
            call system(command)
            
        else if ( anim_type==1 ) then
            anim_file = trim(anim_file_base)//trim(run_name)//".png"

            print *,"Producing montage strip in ",anim_file
            command = "montage -tile 1x50 @"//trim(anim_list_file)//" -geometry +0+0 "//trim(anim_file)
            call system(command)
            anim_files = trim(anim_files)//" "//trim(anim_file)//" "
        else if ( anim_type==-1 ) then
            close(99)
            anim_file = trim(anim_file_base)//trim(run_name)//".dat"
            print *,"Data dumped to",anim_file
        else
            print *,"lolwut"
            stop
        endif
        
    end subroutine do_anim