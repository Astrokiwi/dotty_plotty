program dotty_plotty
    use plot_prams, only: istart,istop,istep,n_pram,run_name,n_runs,runlist,nrot,rot_only,anim_files,anim_type
    use flatkernel, only: kernel_init
    use gcdp_cool, only: read_cool
    use data_vals, only: time
    implicit none
    integer :: intime, rot_step
    integer :: i_pram,i_run
    logical :: read_its
    real(kind=8) :: oldtime
    
    call kernel_init
    
    call read_prams

    !call set_res

    call init_rgb_table
    call read_cool
    
    if ( istop==-1 .and. istart==-1 ) then ! Do all available files
        read_its = .true.
    else
        read_its = .false.
    endif

    if ( anim_type==1 ) then
        anim_files = ""
    endif
    
    do i_run=1,n_runs
        run_name = trim(runlist(i_run))
        print *,i_run," Plotting:",trim(run_name)

        if ( read_its ) then
            call read_stop
        endif

        if ( nrot/=1 .and. nrot/=(istop-istart+1)/istep .and. (istop-istart+1)/istep/=1 ) then
            print *,"Incorrect number of rotation steps"
            print *,"Possibilities are:"
            print *,"  1. nrot=1, i.e. no rotation"
            print *,"  2. nrot=number of steps, i.e. rotation is simultaneous with time"
            print *,"  3. 1 step, and nrot is arbitrary - i.e. stay still and spin around"
            stop
        endif
        if ( (istop-istart+1)/istep==1 .and. nrot>1 ) then
            rot_only = .true.
        else
            rot_only = .false.
        endif
        
        call init_anim
        
        if ( rot_only ) then ! rotate while staying still in time
            call read_data(istart)
            do rot_step = 0,nrot-1
                call set_rot(rot_step)
                do i_pram=1,n_pram
                    call plot_data(i_pram,rot_step)
                end do
                if ( n_pram>1 .and. anim_type>=0 ) then
                    call do_montage(rot_step)
                endif
            end do
            
            call do_anim
        else ! Either dont rotate, or rotate while moving through time
            rot_step = 0
            oldtime = -1
            do intime=istart,istop,istep
                call read_data(intime)
                
                if ( time>oldtime ) then
                    oldtime=time
                    if ( nrot/=1 ) then ! Rotating simultaneously with animating
                        call set_rot(rot_step)
                        rot_step = rot_step + 1
                    endif
                    
                    do i_pram=1,n_pram
                        call plot_data(i_pram,intime)
                    end do
                    if ( n_pram>1 .and. anim_type>=0 ) then
                        call do_montage(intime)
                    endif
                endif
            end do
    
            if ( istop-istart >= 1 .and. anim_type>=0 ) then
                call do_anim
            endif
        endif
        
    end do
    
    if ( anim_type==1 ) then
        call do_strips
    endif

    print *,"Exiting successfully"

end program dotty_plotty

subroutine set_rot(rot_step)
    use plot_prams, only: phi_min, phi_max, theta_min, theta_max, g_phi, g_theta,nrot
    implicit none
    integer :: rot_step

    if ( phi_min==phi_max .or. nrot<=1 ) then
        g_phi = phi_min
    else
        g_phi = (phi_max-phi_min)/(nrot-1)*rot_step + phi_min
    endif
    if ( theta_min==theta_max .or. nrot<=1 ) then
        g_theta = theta_min
    else
        g_theta = (theta_max-theta_min)/(nrot-1)*rot_step + theta_min
    endif
    
    print *,"Global phi=",g_phi,"Global theta=",g_theta

end subroutine set_rot
