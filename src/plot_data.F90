subroutine plot_data(i_pram,intime)
!    use data_vals
    use plot_prams
    implicit none
    
    integer, intent(in) :: i_pram,intime
    character*256 :: outfile
    
    
    
    type (plot_pram_t),pointer :: this_plot ! All the data I need for this particular plot
    
    
    this_plot=>all_prams(i_pram)

#ifndef nopng
    if ( anim_type/=-1 ) then
        call init_plot(i_pram,intime)
    endif
#endif
    
    call dot_plot(i_pram)
    
#ifndef nopng
    if ( anim_type/=-1 ) then
        call finish_plot(i_pram,intime)
    endif
#endif    
end subroutine plot_data

#ifndef nopng
subroutine init_plot(i_pram,intime)
    use plot_prams
	implicit none
	character*256 :: outfile
	integer :: ier,j
    integer :: pgbeg
    integer, intent(in) :: i_pram,intime
    type (plot_pram_t),pointer :: this_plot ! All the data I need for this particular plot

    this_plot=>all_prams(i_pram)

    write(outfile,"('pics/',A,A,I5.5,'_',I2.2,A,'.png')") trim(run_name),trim(this_plot%plot_name),intime,i_pram,trim(pid_str)

    if ( istop-istart >= 1 .or. nrot>1 ) then
        if ( n_pram==1 ) then
            open(unit=42,file=anim_list_file,access='append')
            write(42,*) trim(outfile)
            close(42)
        else
            write(montage_files,"(A,' ',A)") trim(montage_files),outfile
        endif
    else
        if ( n_pram>1 ) then
            write(montage_files,"(A,' ',A)") trim(montage_files),outfile
        endif
    endif
    write(outfile,"(A,'/PNG')") trim(outfile)
    
    
	ier = pgbeg(0,outfile,1,1)
    IF (ier/=1) STOP "Can't open png"

    if ( plot_inches>0. ) then
        call pgpap(plot_inches,1.)
    endif
    
    ! Set the colour spectrum
    do j = 0,ncols+2
        call pgscr(j,rgbtable(this_plot%palette,1,j+1),rgbtable(this_plot%palette,2,j+1),rgbtable(this_plot%palette,3,j+1))
    end do
    call pgscr(255,rgbtable(this_plot%palette,1,ncols+2),rgbtable(this_plot%palette,2,ncols+2),&
            rgbtable(this_plot%palette,3,ncols+2))

    
    ! Plot the boundary
    if ( .not.clean_plot ) then
        call pgsci(255)
    endif

!    call system("export PGPLOT_ENVOPT=I")
!    call system("setenv PGPLOT_ENVOPT I")

    if ( .not.clean_plot ) then
        CALL PGENV(this_plot%rbounds(1,1),this_plot%rbounds(1,2),this_plot%rbounds(2,1),this_plot%rbounds(2,2),1,0)
    else
        CALL PGENV(this_plot%rbounds(1,1),this_plot%rbounds(1,2),this_plot%rbounds(2,1),this_plot%rbounds(2,2),1,-2)
    endif
    
end subroutine init_plot
#endif

subroutine dot_plot(i_pram)
	use plot_prams
	use data_vals
	use flatkernel
	implicit none
	integer,intent(in) :: i_pram
    type (plot_pram_t),pointer :: this_plot ! All the data I need for this particular plot

	integer :: ii
	integer :: isca,ivec,jvec
	
	integer :: ip
	
	real(kind=8) :: dot_value,dot_value2
	real(kind=8) :: dnmax,dnmin
	real(kind=8) :: plot_x,plot_y,plot_z
    real(kind=8) :: costh,sinth,cosph,sinph
	logical :: maxmin_set
	
	integer :: icol
	
	integer :: icol0
	
	
    real(kind=4) :: x1,y1,x2,y2 ! Bounds of drawn region, in pixels
    
    integer,dimension(2) :: ir
    integer :: ix,iy
    integer,dimension(2) :: lb ! how many pixels in x,y direction?
    integer,allocatable,dimension(:,:) :: pix
    
    real(kind=8),allocatable,dimension(:,:,:) :: valgrids
    
    real(kind=8),allocatable,dimension(:,:,:) :: valgrids2

    real(kind=8) :: cell_area, cell_width

    real(kind=8) :: maxh
    real(kind=8) :: phys_sm
    integer :: maxsm,thismaxsm
    integer :: mnix,mniy,mxix,mxiy,six,siy
    real(kind=8) :: smdist,smx,smy

    real(kind=8), dimension(2) :: rcent
    real(kind=8) :: rad,radmax
    integer :: irad
    integer, parameter :: nmean=100
    integer, dimension(nmean) :: mean_n
    real(kind=8),save,dimension(nmean) :: mean_vals

    real(kind=8) :: interp_low,interp_high,interp_f,interp_mean


    real(kind=8) :: norman
    
    real(kind=8) :: zdist2
    
    real(kind=8), parameter :: PI = 3.14159265359
    
    this_plot=>all_prams(i_pram)
    
    if ( this_plot%lx>0 ) then
        ! Arbitrary grid of cells - i.e. if you need a greater column density
        lb(1) = this_plot%lx
        lb(2) = this_plot%ly
    else
#ifndef nopng
        if ( anim_type/=-1 ) then
            ! One cell per pixel of output image
            call pgqvp(3,x1,x2,y1,y2)
            lb(1) = (x2-x1)
            lb(2) = (y2-y1)
        else
#endif
            print *,"Must specify pixels for output data - use grid not dot"
            stop
#ifndef nopng
        endif
#endif
    endif
    
    print *,"Grid size:",lb(1),lb(2)

    if ( anim_type==-1 ) then
        write(99,"(A6,2I8,6E15.6)") this_plot%plot_var,lb,this_plot%rbounds,this_plot%dbounds
    endif
    
    allocate(pix(lb(1),lb(2)))

    if ( this_plot%iover_type==isum_plot .or. &
        this_plot%iover_type==isurf_plot .or. &
        this_plot%iover_type==ismooth_plot ) then
        if ( this_plot%massnorm .and. this_plot%iover_type==ismooth_plot) then
            allocate(valgrids(2,lb(1),lb(2)))
        else
            allocate(valgrids(1,lb(1),lb(2)))
        endif
        valgrids = 0.
    else
        if ( anim_type==-1 ) then
            allocate(valgrids(1,lb(1),lb(2)))
        endif
    endif

	
	ivec = -1 ! Which vector variable? (-1 -> scalar)
	jvec = -1 ! x/y/z component (-1 -> scalar)
	isca = -1 ! Which scalar variable? (-1 -> vector)
	
	select case(this_plot%plot_var(4:4))
		case('x')
			jvec=1
		case('y')
			jvec=2
		case('z')
			jvec=3
	end select
	
	if ( jvec>=0 ) then
		do ii=1,nsca
			if ( vec_codes(ii)==this_plot%plot_var(1:3) ) then
				ivec = ii
				exit
			endif
		end do
	
		if ( ivec == -1 ) then
			print *,"Can not find vector variable ",this_plot%plot_var(1:3)
			print *,"All vectors:",vec_codes(1:nsca)
			stop
		endif
	else if ( this_plot%plot_var/="type" ) then
		do ii=1,nsca
			if ( sca_codes(ii)==this_plot%plot_var ) then
				isca = ii
				exit
			endif
		end do
	
		if ( isca == -1 ) then
			print *,"Can not find scalar variable ",this_plot%plot_var
			print *,"All scalars:",sca_codes(1:nsca)
			stop
		endif
	endif
	! If this_plot%plot_var=="type", then leave ivec,jvec,isca == -1
	
    this_plot%isca = isca
    
	! Set value max/min if not already set
	if ( this_plot%dbounds(1)==this_plot%dbounds(2) ) then
	    this_plot%dbounds(2) = -huge(this_plot%dbounds(2))
	    this_plot%dbounds(1) = huge(this_plot%dbounds(1))
	    do ip=1,np
            if ( itype(ip)==this_plot%iptype .or. this_plot%iptype==itany) then
                if ( isca>=0 ) then
                    dot_value = (sca_data(ip,isca))
                else if ( ivec>=0 ) then
                    dot_value = (vec_data(jvec,ip,ivec))
                else
                    dot_value = ifeedtype(ip)
                endif
                if ( this_plot%massnorm ) then
                    dot_value2 = sca_data(ip,1)
                    dot_value = dot_value * sca_data(ip,1)
                endif
                if ( this_plot%doabs ) then
                    dot_value = abs(dot_value)
                endif
                if ( isca>=0 .and. (isca/=24 .and. isca/=29 .and. isca/=30) .and. &
                    this_plot%iover_type/=isum_plot .and. &
                    this_plot%iover_type/=isurf_plot.and. &
                    this_plot%iover_type/=ismooth_plot.and. &
                    this_plot%iover_type/=ifar_plot.and. &
                        .not.this_plot%dodsub  ) then
                    dot_value = log10(dot_value)
                endif
                this_plot%dbounds(2) = max(this_plot%dbounds(2),dot_value)
                this_plot%dbounds(1) = min(this_plot%dbounds(1),dot_value)
		    endif
	    end do
	endif
	
	if ( this_plot%iover_type==ifar_plot ) then
		! What is zero?
		icol0 = (((0.)-this_plot%dbounds(1))*(ncols-2))/(this_plot%dbounds(2)-this_plot%dbounds(1))+2
		!pix = icol0
        pix = ncols+2
	else if ( this_plot%iover_type==imin_plot ) then
	    pix = ncols+2
	else
		pix = 0
	endif

    if ( anim_type==-1 ) then
        if ( this_plot%iover_type/=isum_plot .and. &
        this_plot%iover_type/=isurf_plot .and. &
        this_plot%iover_type/=ismooth_plot ) then
            if ( this_plot%iover_type==imin_plot ) then
                valgrids(1,:,:) = huge(valgrids(1,0,0))
            else
                valgrids(1,:,:) = -huge(valgrids(1,0,0))
            endif
        endif
    endif
	
	if ( this_plot%iover_type==ismooth_plot .or.&
         this_plot%iover_type==isurf_plot ) then
        cell_area = ((this_plot%rbounds(1,2)-this_plot%rbounds(1,1))/lb(1))&
            * ((this_plot%rbounds(2,2)-this_plot%rbounds(2,1))/lb(2))
    endif
	if ( this_plot%iover_type==ismooth_plot ) then
        cell_width = (this_plot%rbounds(1,2)-this_plot%rbounds(1,1))/lb(1) ! ONLY WORKS FOR SQUARE PLOTS
        maxh = 0.
        do ip=1,np
            if ( itype(ip)==this_plot%iptype .or. this_plot%iptype==itany ) then
                if ( sca_data(ip,4)>maxh ) then
                    maxh = sca_data(ip,4)
                endif
            endif
        end do
        maxsm = maxh/(((this_plot%rbounds(1,2)-this_plot%rbounds(1,1))/lb(1)))+1
        print *,"Max smoothing:",maxh," kpc ",maxsm," px"
!        if ( maxsm>256 ) then
!            maxsm = 256
        if ( maxsm>2048 ) then
            maxsm = 2048
            print *,"Forcing ",maxsm," px"
        endif
	endif

    ! Get mean if necessary
    if ( this_plot%iplot_type==iradmean_plot .or. &
         this_plot%iplot_type==izmean_plot ) then
        rcent = (this_plot%rbounds(1:2,2)+this_plot%rbounds(1:2,1))/2.
        if ( this_plot%iplot_type==iradmean_plot ) then
            radmax = sqrt(sum(((this_plot%rbounds(1:2,2)-this_plot%rbounds(1:2,1))/2)**2.))
        else
            radmax = (this_plot%rbounds(2,2)-this_plot%rbounds(2,1))/2.
        endif

        mean_vals = 0.
        mean_n = 0

        do ip=1,np

            if ( itype(ip)==this_plot%iptype .or. this_plot%iptype==itany ) then
            
                ! Transform r to view from different angles

                costh = cos(this_plot%theta+g_theta)
                sinth = sin(this_plot%theta+g_theta)
                cosph = cos(this_plot%phi+g_phi)
                sinph = sin(this_plot%phi+g_phi)

                !is this more correct?
                plot_x = r_p(1,ip)*costh       +   r_p(2,ip)*sinth
                plot_y =-r_p(1,ip)*sinth       +   r_p(2,ip)*costh

                plot_z =-plot_y*sinph          +   r_p(3,ip)*cosph
                plot_y = plot_y*cosph          +   r_p(3,ip)*sinph ! ????

                if ( plot_z>this_plot%zbounds(1) .and. plot_z<this_plot%zbounds(2) ) then
                    if ( this_plot%isca>=0 ) then
                        dot_value = (sca_data(ip,this_plot%isca))
                    else
                        print *,"can't do vector meanness yet"
                        stop
                    endif


                    if ( this_plot%iplot_type==iradmean_plot ) then
                        rad = sqrt(((plot_x-rcent(1))**2+(plot_y-rcent(2))**2))/radmax
                    else
                        rad = abs(plot_y-rcent(2))/radmax
                    endif
                    irad = rad*nmean+1
                    if ( irad>=1 .and. irad<=nmean ) then
                        mean_n(irad) = mean_n(irad)+1
                        !if ( this_plot%dodsub ) then
                        !    mean_vals(irad) = mean_vals(irad)+dot_value-this_plot%dsub
                        !else
                            mean_vals(irad) = mean_vals(irad)+dot_value
                        !endif
                    endif
                endif
            endif
        end do

        if ( this_plot%iover_type==isurf_plot ) then
            do irad=1,nmean
                mean_vals(irad) = mean_vals(irad)/( ((irad)**2-(irad-1)**2)*PI*(radmax/nmean)**2)
                
                mean_vals(irad) = mean_vals(irad)*cell_area ! weightings to subtract from added values
            end do
        else
            do irad=1,nmean
                if ( mean_n(irad)>0 ) then
                    mean_vals(irad) = mean_vals(irad)/mean_n(irad)
                endif
            end do
        endif
    endif
    
	! TODO: Move main loop (building up "pix") into different routine
	maxmin_set = .false.
	do ip=1,np
!	    if (mod(ip,10000)==0 ) then
!	        print *,ip
!	    endif
	
		! Do a "cut" by particle type (e.g. density is only for gas particles)
		if ( itype(ip)==this_plot%iptype .or. this_plot%iptype==itany ) then
		
			! Transform r to view from different angles
			!plot_x = r_p(1,ip)*cos(this_plot%theta+g_theta) + r_p(2,ip)*sin(this_plot%theta+g_theta)
			!plot_y = -r_p(1,ip)*cos(this_plot%phi+g_phi)*sin(this_plot%theta+g_theta) + r_p(2,ip)*cos(this_plot%phi+g_phi)*cos(this_plot%theta+g_theta)+r_p(3,ip)*sin(this_plot%phi+g_phi)
            costh = cos(this_plot%theta+g_theta)
            sinth = sin(this_plot%theta+g_theta)
            cosph = cos(this_plot%phi+g_phi)
            sinph = sin(this_plot%phi+g_phi)
            !plot_x = r_p(1,ip)*costh    +   r_p(2,ip)*sinth*cosph   +   r_p(3,ip)*sinth*sinph
            !plot_y =-r_p(1,ip)*sinth    +   r_p(2,ip)*costh*cosph   +  -r_p(3,ip)*costh*sinph
            !plot_z =                        r_p(2,ip)*sinph         +  r_p(3,ip)*cosph

            !is this more correct?
            plot_x = r_p(1,ip)*costh       +   r_p(2,ip)*sinth
            plot_y =-r_p(1,ip)*sinth       +   r_p(2,ip)*costh

            plot_z =-plot_y*sinph          +   r_p(3,ip)*cosph
            plot_y = plot_y*cosph          +   r_p(3,ip)*sinph

            if ( plot_z>this_plot%zbounds(1) .and. plot_z<this_plot%zbounds(2) ) then


                if ( isca>=0 ) then
                    dot_value = (sca_data(ip,isca))
                else if ( ivec>=0 ) then
                    dot_value = (vec_data(jvec,ip,ivec))
                else
                    dot_value = ifeedtype(ip)
                endif

                if ( this_plot%massnorm ) then
                    dot_value2 = sca_data(ip,1)
                    dot_value = sca_data(ip,isca)
                endif
    		

                if ( this_plot%doabs ) then
                    dot_value = abs(dot_value)
                endif
                if ( this_plot%iplot_type==iradmean_plot .or. &
                     this_plot%iplot_type==izmean_plot ) then
                    if ( this_plot%iplot_type==iradmean_plot ) then
                        rad = sqrt(((plot_x-rcent(1))**2+(plot_y-rcent(2))**2))/radmax
                    else
                        rad = abs(plot_y-rcent(2))/radmax
                    endif

                    irad = rad*nmean+1
                    if ( irad>=1 .and. irad<=nmean ) then
                        if ( mean_n(irad)>0 ) then
                            interp_low = mean_vals(irad)
                            if ( irad==nmean ) then
                                interp_high = mean_vals(irad)
                            else
                                interp_high = mean_vals(irad+1)
                            endif
                            interp_f = rad*nmean+1-irad
                            interp_mean = interp_low*(1.-interp_f) + interp_high*interp_f

                            if ( this_plot%dodsub ) then
                                !dot_value = (dot_value-this_plot%dsub)-interp_mean
                                dot_value = dot_value-interp_mean
                            else
                                dot_value = dot_value/interp_mean
                            endif
                        else
                            print *,"I should not be here, it's mean"
                            stop
                        endif
                    endif
                endif

                ! Abitrarily do the log
                ! NOT FINISHED
                if ( isca>=0 .and. (isca/=24 .and. isca/=29 .and. isca/=30) .and. &
                    this_plot%iover_type/=isum_plot .and. &
                    this_plot%iover_type/=isurf_plot .and. &
                    this_plot%iover_type/=ismooth_plot.and. &
                        this_plot%iover_type/=ifar_plot.and. &
                        .not.this_plot%dodsub ) then
                    dot_value = log10(dot_value)
                endif
    			
    			!call draw_dot(plot_x,plot_y,dot_value,i_pram)

                ! Find which colour this pixel should be
    			icol = (((dot_value)-this_plot%dbounds(1))*(ncols-2))/(this_plot%dbounds(2)-this_plot%dbounds(1))+2
    			
    			if ( icol>ncols ) icol=ncols

                if ( icol<2 .and. (this_plot%palette==4 .or. this_plot%palette==5) ) then
                    icol = 2
                endif

!			call pgsci(icol)
    			if ( icol>=2 .or. this_plot%iover_type==isum_plot .or. &
                    this_plot%iover_type==isurf_plot .or. &
                    this_plot%iover_type==ismooth_plot) then
    				!if ( (dot_value>=this_plot%dbounds(1) .and. dot_value<=this_plot%dbounds(2))) then
    					! add to array
    					ir(1) = ((plot_x-this_plot%rbounds(1,1))/(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*lb(1))+1
    					ir(2) = ((plot_y-this_plot%rbounds(2,1))/(this_plot%rbounds(2,2)-this_plot%rbounds(2,1))*lb(2))+1

    					if ( all(ir>0) .and. all(ir<lb) .or. &
                            this_plot%iover_type==ismooth_plot) then
    						select case(this_plot%iover_type)
                                ! Do these assignments in terms of *colour of pixel*
    							case(imax_plot)
                                    if ( anim_type==-1 ) then
                                        valgrids(1,ir(1),ir(2))=max(dot_value,valgrids(1,ir(1),ir(2)))
                                    else
                                        pix(ir(1),ir(2))=max(icol,pix(ir(1),ir(2)))
                                    endif
    							case(imin_plot)
                                    if ( anim_type==-1 ) then
                                        valgrids(1,ir(1),ir(2))=min(dot_value,valgrids(1,ir(1),ir(2)))
                                    else
                                        pix(ir(1),ir(2))=min(icol,pix(ir(1),ir(2)))
                                    endif
    							case(ifar_plot)
                                    if ( anim_type==-1 ) then
                                        print *,"Data dump doesn't work with far plot yet"
                                        stop
                                    endif

                                    if (pix(ir(1),ir(2))==ncols+2) then
                                        pix(ir(1),ir(2))=icol
                                    else
        								if ( abs(pix(ir(1),ir(2))-icol0)<abs(icol-icol0) ) then
        									pix(ir(1),ir(2))=icol
        								endif
                                    endif
                                ! Do these assignments in terms of *actual physical value* - convert into colour later
                                case(isum_plot,isurf_plot)
                                    valgrids(1,ir(1),ir(2)) = valgrids(1,ir(1),ir(2)) + dot_value
                                case(ismooth_plot)
                                    !do ii=1,2
                                    thismaxsm = min(maxsm,int(sca_data(ip,4)*lb(2)/(this_plot%rbounds(2,2)-this_plot%rbounds(2,1))))
                                    phys_sm = thismaxsm*cell_width

                                    if ( this_plot%smoo3d ) then
                                        zdist2 = ((plot_z-this_plot%smooslice)*lb(2) &
                                          /(this_plot%rbounds(2,2)-this_plot%rbounds(2,1)))**2
                                    else
                                        zdist2 = 0.d0
                                    endif

                                    !thismaxsm = max(2,thismaxsm)
                                    if ( thismaxsm<1 ) then
                                        if ( zdist2<=1. ) then ! otherwise it's too far away to see at all
                                            ! Everything is within one pixel - just add the value here
                                            ! but only if we are within the grid
                                            if ( all(ir>0) .and. all(ir<lb) ) then
                                             if ( this_plot%massnorm ) then
                                                 if ( this_plot%smoo3d ) then
                                                    valgrids(2,ir(1),ir(2)) = valgrids(2,ir(1),ir(2)) + dot_value2/sca_data(ip,4)**3
                                        valgrids(1,ir(1),ir(2)) = valgrids(1,ir(1),ir(2)) + dot_value * dot_value2/sca_data(ip,4)**3
                                                 else
                                                     valgrids(2,ir(1),ir(2)) = valgrids(2,ir(1),ir(2)) + dot_value2
                                                     valgrids(1,ir(1),ir(2)) = valgrids(1,ir(1),ir(2)) + dot_value * dot_value2
                                                 endif
                                             else
                                                 if ( this_plot%smoo3d ) then
                                                     valgrids(1,ir(1),ir(2)) = valgrids(1,ir(1),ir(2)) + dot_value/sca_data(ip,4)**3
                                                 else
                                                     if ( this_plot%smoo3d ) then
                                                     valgrids(1,ir(1),ir(2)) = valgrids(1,ir(1),ir(2)) + dot_value/sca_data(ip,4)**3
                                                     else
                                                         valgrids(1,ir(1),ir(2)) = valgrids(1,ir(1),ir(2)) + dot_value
                                                     endif
                                                 endif
                                                 if ( valgrids(1,ir(1),ir(2))>1.e20 ) then
                                                     print *,"BIG INNER VALUE",ir(1),ir(2),dot_value,sca_data(ip,4)
                                                 endif 
                                             endif
                                            endif
                                        endif
                                    else
                                        ! Do smoothing
                                        mnix = max(1,ir(1)-thismaxsm)
                                        mniy = max(1,ir(2)-thismaxsm)
                                        mxix = min(lb(1),ir(1)+thismaxsm)
                                        mxiy = min(lb(2),ir(2)+thismaxsm)

                                        !mnix = ir(1)-thismaxsm
                                        !mniy = ir(2)-thismaxsm
                                        !mxix = ir(1)+thismaxsm
                                        !mxiy = ir(2)+thismaxsm
                                        !norman = 0.
                                        do six=mnix,mxix
                                            do siy=mniy,mxiy
                                                !if ( six/=ir(1) .and. siy/=ir(2) ) then
                                                    !smx = (six-1.)/lb(1)*(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))+this_plot%rbounds(1,1)
                                                    !smy = (siy-1.)/lb(2)*(this_plot%rbounds(2,2)-this_plot%rbounds(2,1))+this_plot%rbounds(2,1)
                                                    !smdist = sqrt((smx-plot_x)**2+(smy-plot_y)**2)
                                                    !smdist = smdist/sca_data(ip,4) ! divide by smoothing length
                                                
                                                    smdist = sqrt(real((six-ir(1))**2+(siy-ir(2))**2)+zdist2)/thismaxsm
                                        ! moved to the left because of gfortran limits
                if ( this_plot%massnorm ) then
                    ! We are doing a mass-weighted smoothed sum

                    if ( this_plot%smoo3d ) then
                    ! Add up masses
                        valgrids(2,six,siy) = valgrids(2,six,siy) + dot_value2 * kern(smdist)/kern_norm/phys_sm**3
                        ! Add up (scalar value) with mass weighting
                        valgrids(1,six,siy) = valgrids(1,six,siy) + dot_value * dot_value2 * kern(smdist)/kern_norm/phys_sm**3
                    else
                        ! Add up masses
                        valgrids(2,six,siy) = valgrids(2,six,siy) + dot_value2 * fkern(smdist)*cell_area/phys_sm**2
                        ! Add up (scalar value) with mass weighting
                        valgrids(1,six,siy) = valgrids(1,six,siy) + dot_value * dot_value2 * fkern(smdist)*cell_area/phys_sm**2
                    endif
                else 
                    if ( this_plot%smoo3d ) then
                        valgrids(1,six,siy) = valgrids(1,six,siy) + dot_value * kern(smdist)/kern_norm/phys_sm**3
                    else
                        valgrids(1,six,siy) = valgrids(1,six,siy) + dot_value * fkern(smdist)*cell_area/phys_sm**2
                    endif
                endif
                !norman = norman + dot_value * fkern(smdist) * cell_area/phys_sm**2
                if ( valgrids(1,six,siy)>1.e20 ) then
                    print *,"BIG VALUE",six,siy,dot_value,kern(smdist),smdist,kern_norm,phys_sm
                endif 

                                                !endif
                                            end do
                                        end do

                                        !print *,ip,sca_data(ip,4),thismaxsm,phys_sm,norman
                                        !sca_data(ip,4) = sca_data(ip,4)/2.
                                        
                                    endif
                                    !end do
                                    !stop
    						end select
    					endif
    				!endif
    			endif


                ! Find the max/min
                if ( .not.maxmin_set ) then
                    maxmin_set = .true.
                    dnmax = (dot_value)
                    dnmin = (dot_value)
                endif

            endif ! zcut

			if ( (dot_value)>dnmax ) dnmax = (dot_value)
			if ( (dot_value)<dnmin ) dnmin = (dot_value)
		endif

	end do




    
    ! Convert from grid of physical values into colours
    ! (another subroutine?)
    if ( this_plot%iover_type==isum_plot .or. &
         this_plot%iover_type==isurf_plot .or. &
         this_plot%iover_type==ismooth_plot  ) then
        if ( this_plot%iover_type==isurf_plot .or. &
         this_plot%iover_type==ismooth_plot ) then
            cell_area = ((this_plot%rbounds(1,2)-this_plot%rbounds(1,1))/lb(1)) &
            * ((this_plot%rbounds(2,2)-this_plot%rbounds(2,1))/lb(2))
        endif

        maxmin_set = .false. ! We want to find the max/min of *final* summed values, not of the component densities etc
        do ix=1,lb(1) ! I'd love to loop with ir(1:2) as the loopy variables, but that's not allowed :(
            do iy=1,lb(2)
                ir(1) = ix
                ir(2) = iy
                if ( valgrids(1,ir(1),ir(2))/=0 ) then
                    if ( this_plot%massnorm ) then
                        dot_value = valgrids(1,ir(1),ir(2))/valgrids(2,ir(1),ir(2))
                    else
                        dot_value = valgrids(1,ir(1),ir(2))
                    endif
                    
                    if ( this_plot%iover_type==isurf_plot .or. &
                         this_plot%iover_type==ismooth_plot ) then
                        if ( .not. this_plot%massnorm .and. .not. this_plot%smoo3d) then
                            dot_value = dot_value/cell_area ! divide to get column density
                        endif
                    endif
                    
                    ! CHECK IF WE SHOULD LOG OR NOT
                    ! NOT FINISHED
                    if ( (isca/=24 .and. isca/=29 .and. isca/=30).and. &
                        .not.this_plot%dodsub  ) then
                        dot_value = log10(dot_value)
                        if ( dot_value/=dot_value ) then
                            print *,"NaN!",ix,iy,dot_value
                        endif
                    endif
                    
                    if ( anim_type==-1 ) then
                        !write(99,"(E15.6)") dot_value
                        write(99,*) dot_value
                    else
                        ! Get colour as normal
                        icol = (((dot_value)-this_plot%dbounds(1))*(ncols-2))/(this_plot%dbounds(2)-this_plot%dbounds(1))+2
                        if ( icol>ncols ) icol=ncols
                        if ( icol>=2 ) then
                            pix(ir(1),ir(2)) = icol
                        endif
                        ! Find the max/min
                        if ( .not.maxmin_set ) then
                            maxmin_set = .true.
                            dnmax = (dot_value)
                            dnmin = (dot_value)
                        endif
            
                        if ( (dot_value)>dnmax ) dnmax = (dot_value)
                        if ( (dot_value)<dnmin ) dnmin = (dot_value)
                    endif
                else if ( anim_type==-1 ) then
                    write(99,*) huge(dot_value)
                endif
            end do
        end do
        deallocate(valgrids)
    else if ( anim_type==-1 ) then
        do ix=1,lb(1) ! I'd love to loop with ir(1:2) as the loopy variables, but that's not allowed :(
            do iy=1,lb(2)
               ! write(99,"(E15.6)") valgrids(1,ix,iy)
                write(99,*) valgrids(1,ix,iy)
            end do
        end do
        deallocate(valgrids)
    endif
	
!	call pgpixl(pix,lb(1),lb(2),1,lb(1),1,lb(2),x1,x2,y1,y2)
    if ( anim_type/=-1 ) then
#ifndef nopng
    	call pgpixl(pix,lb(1),lb(2),1,lb(1),1,lb(2),&
    	    this_plot%rbounds(1,1),this_plot%rbounds(1,2),this_plot%rbounds(2,1),this_plot%rbounds(2,2))
#else
		print *,"ANIM TYPE IS INCORRECT FOR NO PGPLOT BUILD"
		stop
#endif
    endif
	
	print *,"max,min:",dnmax,dnmin
	
	deallocate(pix)
end subroutine dot_plot

#ifndef nopng
subroutine finish_plot(i_pram,intime)
	use plot_prams
    use data_vals, only: sca_units, sca_codes, time
	implicit none
    type (plot_pram_t),pointer :: this_plot ! All the data I need for this particular plot
    integer,intent(in) :: i_pram,intime
    
    integer :: icol
    character*256 :: text
    real(kind=4),dimension(2) :: line_x,line_y
    
    integer :: jj
    
    real(kind=4) :: cbtick_y,cbtick_val
	
	if ( .not.clean_plot ) then
    
        this_plot=>all_prams(i_pram)

        ! Print label for this particular plot at this time
        call pgsci(255)
        if ( cleanish_plot ) then
            write(text,'("Time = ",F10.3," Myr")') time*.471e3
        else
            write(text,'(A," ",A4," ",I5.5)') trim(this_plot%plot_name),this_plot%plot_var,intime
        endif
        call pgtext((this_plot%rbounds(1,1)+this_plot%rbounds(1,2))/2.,&
            this_plot%rbounds(2,2)+(this_plot%rbounds(2,2)-this_plot%rbounds(2,1))/20.,text)


        if ( .not.cleanish_plot ) then


            ! Plot axis
            do jj=2,ncols
                call pgsci(jj)
                line_y(1) = this_plot%rbounds(2,1)+((jj-2)*(this_plot%rbounds(2,2)-this_plot%rbounds(2,1)))/(ncols-1)
                line_y(2) = this_plot%rbounds(2,1)+((jj-1)*(this_plot%rbounds(2,2)-this_plot%rbounds(2,1)))/(ncols-1)
        !		line_x(1) = this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*1.05
                line_x(1) = this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*.95
                line_x(2) = this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*1.00
                call pgrect(line_x(1),line_x(2),line_y(1),line_y(2))
            end do
            call pgsci(255)
        
        !    call pgaxis("NL",this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*1.05,this_plot%rbounds(2,1),this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*1.05,this_plot%rbounds(2,2),this_plot%dbounds(1),this_plot%dbounds(2),.0,.0,4,4,.5,2,180)

            write(text,'(F8.2)') this_plot%dbounds(1)
            call pgtext(this_plot%rbounds(1,2),this_plot%rbounds(2,1),text)

            do jj=1,3
                cbtick_y = this_plot%rbounds(2,1) + ((this_plot%rbounds(2,2)-this_plot%rbounds(2,1))/(3+1))*jj
                cbtick_val = this_plot%dbounds(1) + ((this_plot%dbounds(2)-this_plot%dbounds(1))/(3+1))*jj
                write(text,'(F8.2)') cbtick_val
                call pgtext(this_plot%rbounds(1,2),cbtick_y,text)
            end do


        !     Cblabel axis, with units
            if ( this_plot%isca/=-1 ) then
                if ( this_plot%iover_type==isurf_plot .or. this_plot%iover_type==ismooth_plot ) then
                    text = "log10("//trim(sca_codes(this_plot%isca))//"surf)"
                    text = trim(text)//" ("//trim(sca_units(this_plot%isca))//"/kpc**2)"
                else
                    text = "log10("//trim(sca_codes(this_plot%isca))//")"
                    text = trim(text)//" ("//trim(sca_units(this_plot%isca))//")"
                endif
            endif
        
            call pgptext(this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*1.2,&
            (this_plot%rbounds(2,2)+this_plot%rbounds(2,1))/2.,270.,.5,text)

            write(text,'(F8.2)') this_plot%dbounds(2)
            call pgtext(this_plot%rbounds(1,2),this_plot%rbounds(2,2),text)

        endif

        text = "x (kpc)"
        call pgptext((this_plot%rbounds(1,2)+this_plot%rbounds(1,1))/2.,&
        this_plot%rbounds(2,1)+(this_plot%rbounds(2,2)-this_plot%rbounds(2,1))*(-.1),0.,.5,text)

        text = "y (kpc)"
        call pgptext(this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*(-0.1),&
        (this_plot%rbounds(2,2)+this_plot%rbounds(2,1))/2.,90.,.5,text)

    endif

	CALL PGEND

end subroutine finish_plot
#endif
