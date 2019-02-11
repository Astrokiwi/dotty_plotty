subroutine read_prams
    use plot_prams!, only: all_prams,istart,istop,istep,dirname,run_name,plot_pram_t,n_pram,idot_plot,isurf
    use data_vals
    implicit none
    character*256 pram_file
    character*64 plot_type_str,plot_over_str,plot_ptype_str
    character*1024 instr
    type (plot_pram_t),pointer :: this_plot
    
    integer :: i_pram
    
    integer :: ipos1,ipos2
    
    integer :: ii

    call getarg(1,pram_file)
    
    if ( pram_file=="" ) then
    write(6,'(''Input parameter file: '',$)')
    read(5,'(a)') pram_file
    endif
	open(unit=42,file=pram_file)
	read(42,"(A)") dirname
	read(42,"(A)") instr
    ! Separate input string into space-separated list of run names
    ipos1 = 1
    n_runs = 0
    ipos2 = 1
    do while ( ipos2/=0 )
        ipos2 = index(trim(instr(ipos1:))," ")
        n_runs = n_runs + 1
        if ( ipos2==0 ) then
            runlist(n_runs) = instr(ipos1:)        
        else
            runlist(n_runs) = instr(ipos1:(ipos1+ipos2-2))
            ipos1 = ipos2 + ipos1
        endif
    end do
!    do ii=1,n_runs
!        print *,trim(runlist(ii))
!    end do
!    stop
    read(42,*) istart,istop,istep
    read(42,*) phi_min,phi_max,theta_min,theta_max,nrot
    ! Convert to radians
    phi_min = phi_min*0.01745329251
    phi_max = phi_max*0.01745329251
    theta_min = theta_min*0.01745329251
    theta_max = theta_max*0.01745329251
    
    !read(42,*) palette_type
    read(42,*) plot_inches
	read(42,*) n_pram
	read(42,*) anim_type
    read(42,*) ii
    select case(ii)
    case(0)
        tidedisc = .false.
        xydisc = .false.
        centdisconly = .false.
    case(1)
        tidedisc = .false.
        xydisc = .true.
        centdisconly = .false.
    case(2)
        tidedisc = .true.
        xydisc = .false.
        centdisconly = .false.
    case(3)
        tidedisc = .false.
        xydisc = .false.
        centdisconly = .true.
    end select

    read(42,"(A)") select_file
    if ( select_file(1:3)=="all" ) then
        do_selection = .false.
    else
        do_selection = .true.
    endif
    
    read(42,"(A)") anim_file_base
    if ( anim_type==0 ) then
        print *,"Standard gif animation"        
    else if (anim_type==1) then
        print *,"Montage strips - not a real animation"
    else if (anim_type==-1) then
        print *,"Dumping everything to data file!"
    else
        print *,"anim_type==0 or ==1 or ==-1 !!"
        stop
    endif

#ifdef nopng
	print *,"FORCING TO DATA DUMP - PGPLOT NOT ACTIVATED IN THIS BUILD"
	anim_type = -1
#endif

	if ( n_pram>n_pram_max ) then
	    print *,n_pram_max," max plots per file"
	    stop
	endif
    do i_pram=1,n_pram
        this_plot=>all_prams(i_pram)
        read(42,"(A)") this_plot%plot_var ! Plot rho/m/z/?
        read(42,"(A)") plot_over_str ! How do we overwrite points?
        read(42,"(A)") instr ! grid (with dimensions) or as many dots as possible?
        read(42,"(A)") plot_ptype_str ! gas/dark/star
        read(42,*) this_plot%palette
        read(42,*) this_plot%phi,this_plot%theta ! in degrees
        read(42,*) this_plot%dbounds(1),this_plot%dbounds(2) ! value bounds 
        if ( this_plot%dbounds(1)==this_plot%dbounds(2) ) then
            print *,"Automatic scaling!"
        endif


        read(42,*) this_plot%rbounds(1,1),this_plot%rbounds(1,2),this_plot%rbounds(2,1),&
                   this_plot%rbounds(2,2),this_plot%zbounds(1),this_plot%zbounds(2) ! x/y/z bounds
        
        read(instr,*) plot_type_str

        this_plot%iplot_type = -1
        this_plot%lx = -1
        
        this_plot%dodsub = .false.
        select case(plot_type_str)
        case("dot")
          print *,"Dot plot (i.e. minimum possible bin-size) of '",this_plot%plot_var,"'"
          this_plot%iplot_type = idot_plot
        case("grid")
          this_plot%iplot_type = igrid_plot
          read(instr,*) plot_type_str,this_plot%lx,this_plot%ly
          print *,this_plot%lx,"x",this_plot%ly," grid plot of'",this_plot%plot_var,"'"
        case("orad")
          this_plot%iplot_type = iradmean_plot
          print *,"Dot plot (i.e. minimum possible bin-size) of '",this_plot%plot_var,"', removing radial mean"
        case("orgd")
          this_plot%iplot_type = iradmean_plot
          read(instr,*) plot_type_str,this_plot%lx,this_plot%ly
          print *,this_plot%lx,"x",this_plot%ly," grid plot of '",this_plot%plot_var,"', removing radial mean"
        case("oz")
          this_plot%iplot_type = izmean_plot
          print *,"Dot plot (i.e. minimum possible bin-size) of '",this_plot%plot_var,"', dividing out vertical mean"
        case("ozd")
          this_plot%iplot_type = izmean_plot
          !read(instr,*) plot_type_str,this_plot%dsub
          this_plot%dodsub = .true.
          print *,"Dot plot (i.e. minimum possible bin-size) of '",this_plot%plot_var,"', subtracting vertical mean"
          !print *,"And subtracting ",this_plot%dsub," from the value everywhere"
        end select
        
        if ( this_plot%iplot_type==-1 ) then
            print *,"Plot type not found!"
            print *,"dot/grid/orad/oz/ozd"
            stop
        endif
        
        this_plot%iover_type = -1
        this_plot%doabs = .false.
        this_plot%massnorm = .false.
        this_plot%smoo3d = .false.

        select case(plot_over_str)
        case("max")
          print *,"Plotting maximum value in each pixel"
          this_plot%iover_type = imax_plot
        case("min")
          print *,"Plotting min value in each pixel"
          this_plot%iover_type = imin_plot
        case("up")
          print *,"Plotting closest point to camera in each pixel"
          this_plot%iover_type = iup_plot
          print *,"NOT YET IMPLEMENTED"
          stop
        case("sum")
          print *,"Summing all points in each pixel"
          this_plot%iover_type = isum_plot
        case("far")
          print *,"Plotting largest data value (i.e. furthest, positive or negative, from zero) in each pixel"
          this_plot%iover_type = ifar_plot
        case("surf")
          print *,"Summing all points in each pixel and dividing by physical area of pixel (i.e. column density)"
          this_plot%iover_type = isurf_plot
        case("smoo")
          print *,"SPH smoothed column density"
          this_plot%iover_type = ismooth_plot
        case("smow")
          print *,"SPH smoothed column density, divided by mass column density"
          this_plot%iover_type = ismooth_plot
          this_plot%massnorm = .true.
        case("smos")
          print *,"SPH smoothed column density, but take abs first"
          this_plot%iover_type = ismooth_plot
          this_plot%doabs = .true.
        end select
        
        if ( plot_over_str(1:4)=="sm3d" ) then
          print *,"SPH smoothed *slice*"
          this_plot%iover_type = ismooth_plot
          this_plot%smoo3d = .true.
          read(plot_over_str,*) instr(1:4),this_plot%smooslice
          print *,"Slicing at LoS z=",this_plot%smooslice
          plot_over_str = instr(1:4)
        else if ( plot_over_str(1:4)=="sw3d" ) then
          print *,"SPH smoothed *slice*, weighted by 3d density"
          this_plot%iover_type = ismooth_plot
          this_plot%smoo3d = .true.
          read(plot_over_str,*) instr(1:4),this_plot%smooslice
          print *,"Slicing at LoS z=",this_plot%smooslice
          plot_over_str = instr(1:4)
          this_plot%massnorm = .true.
        endif

        if ( this_plot%iover_type==-1 ) then
            print *,"Plot over-write type not found!"
            print *,"max/min/up/sum/far/surf/smoo/smow/smos"
            stop
        endif

        this_plot%iptype = -42

        select case(plot_ptype_str)
        case("star")
            this_plot%iptype = itstar
        case("gas")
            this_plot%iptype = itgas
        case("dark")
            this_plot%iptype = itdark
        case("any")
            this_plot%iptype = itany
!        case("feed") ! Feedback particles are a type of gas particle... maybe have a logical array so a particle can be multiple types?
!            this_plot%iptype = itfeed
        end select
        
        if ( this_plot%iptype==-42 ) then
            print *,"Particle type not found"
            print *,"star/gas/dark/any"
            stop
        endif

        this_plot%phi = this_plot%phi*0.01745329251 ! Convert to radians
        this_plot%theta = this_plot%theta*0.01745329251 ! Convert to radians
        
        write(this_plot%plot_name,"(3A)") trim(plot_type_str),trim(plot_over_str),trim(plot_ptype_str) ! Set up name of plotfile for output

    end do
    
    close(42)
end subroutine read_prams

subroutine read_stop
    use plot_prams
    implicit none
    character*512 pram_file,dumpfile
    integer :: idump,ierr
        
    write(dumpfile,"(A,'/',A,'/diskev/output/ana/ostep.dat')") trim(dirname),trim(run_name)
    open(unit=50,file=dumpfile,status='old')
    istop = 0
    istart = 1
    ierr = 0
    do while (ierr==0 )
        read(50,*,iostat=ierr) idump
        if ( ierr==0 ) then
            istop = istop + 1
        endif
    end do
    close(50)

end subroutine read_stop

subroutine init_anim
    use plot_prams, only: run_name,anim_list_file,istart,istop,nrot,anim_type,anim_file_base,pid_str
    !USE IFPORT ! needed on ifort
    implicit none
    character*128 :: anim_file
!     character*64 :: pid_str
    integer :: pid

    if (anim_type/=-1) then
        if ( istop-istart >= 1 .or. nrot>1 ) then
            ! Create blank file containing list of animation files
            pid = getpid()
            write(pid_str,"(I10.10)") pid
            !write(anim_list_file,"('animfiles',A,A,'.dat')") trim(run_name)
            !anim_list_file = "animfiles"//trim(anim_file_base)//trim(run_name)//".dat"
            anim_list_file = "animfiles"//trim(run_name)//"_"//trim(pid_str)//".dat"
            open(unit=42,file=anim_list_file,action='write')
            write(42,*)
            close(42)
        endif
    else
        anim_file = trim(anim_file_base)//trim(run_name)//".dat"
        print *,"Writing data to",anim_file
        open(unit=99,file=anim_file,action='write')
    endif
end subroutine init_anim
