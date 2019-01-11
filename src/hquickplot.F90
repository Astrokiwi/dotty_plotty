      program hquickplot
      include 'psize.inc'
      include 'pinfo.inc'
      include 'itype.inc'
      include 'rvarrays.inc'
      include 'iodefs.inc'
      include 'units.inc'
!      include 'rgbjames.inc'
      include 'rgbtipsy.inc'
!      include 'rgbredblue.inc'

	  type plot_pram_t ! Parameters for a plot
	  	integer :: iplot_type
	  	character*4 :: plot_name
	  	real :: dbounds(2)
	  	real :: rbounds(2,2)
	  	logical :: sage_flag
	  	logical :: load_star
 	    real*8, allocatable :: grid_val(:,:), grid_val2(:,:), grid_val3(:,:), grid_val4(:,:)
 	    real*8, allocatable :: rad_val(:)
	    integer :: grid_L
	    integer :: rad_L
	  end type plot_pram_t

     
      INTEGER I, IER, PGBEG, j
!      integer,parameter :: npt_max=1.e6
!      real pt(npt_max,255,2),npt(255)
!      real rbounds(2,2) ! Bounds in Mpc
!      real dbounds(2) ! Bounds in log density (units!)
      real line_x(2),line_y(2)

	  integer inrun,istart,istop,istep
	  integer iplot
	  
!	  integer iplot_type
	  integer,parameter :: istar_plot=1,itemp_plot=2,idens_plot=3,isage_plot=4,isagp_plot=5,ipage_plot=6
	  integer,parameter :: ivelx_plot=7,ively_plot=8,ivelz_plot=9,ivrad_plot=10,ivrot_plot=11
	  integer,parameter :: idspz_plot=12,idspc_plot=13,iicld_plot=14,icurl_plot=15,idsc2_plot=16
	  integer,parameter :: iover_plot=17,isfrd_plot=18,isfrg_plot=19
!	  logical :: sage_flag
	  
!	  real, allocatable :: grid_val(:,:), grid_val2(:,:)
!	  integer :: grid_L
	  character*4 junk
	  integer :: ijunk
	  integer :: grid_ij(2) ! Coordinate on grid
	  integer :: grid_ir ! Radial bin
	  real    :: r_ij(2),pix_size(2)
	  logical :: onplot
	  
	  real    :: rad, maxrad
	  real    :: rad1,rad2
	  real    :: theta1,theta2
	  
	  real :: sage_bot,sage_top,sage_min,sage_max

	  real :: this_scalar
	  logical :: maxmin_set

      character*11 hydra_in
      character*64 pram_file
      character*256 outdir
      character*300 outfile
      character*148 cloudarrayfile
      character*148 boundlistfile
      character*90 itypelongfile
      character*11 itypeshortfile
      character*512 command
      character*148 plot_type
      character*6 p_str
      character*2 p_str2
      character*4 irun_str
      character*128 movie
      character*64 time_text
      
      logical :: follow_points
      real, dimension(2) :: centrep,rsize
      integer :: read_itime
      integer :: idump, rdump
      
      real roff(2)
      
      logical make_movie,make_graphics
!
      parameter(nbuf=max(12,Nmax/8))
      real realbuf(nbuf)
      integer intbuf(nbuf),index(nmax)
      double precision dbtime
      equivalence (realbuf,intbuf),(realbuf,dbtime)
      
      integer, parameter :: n_pram_max = 6
      integer :: n_pram, i_pram
      
      integer :: my_cloud,the_cloud
      integer :: itime_cloudfile
      integer, allocatable :: itype2(:)
      integer :: biggest_n,nitype
      
      integer :: ispace
      
      integer :: ix,iy
      
      type (plot_pram_t),pointer :: this_plot ! All the data I need for this particular plot
      type (plot_pram_t),target,dimension(n_pram_max) :: all_prams
      
      integer :: COUNT_DOOM
      
      COUNT_DOOM = 0
      
      follow_points = .false.

!      dirname='/disc41/williams'
       k = 0
       rgbtable(:,1) = 1.0
       do j = 2,ncols+1
         do i = 1,3
           k = k + 1
           rgbtabletmp(i,j) = float(rgb2(k))/256.0
         end do
       end do
       rgbtable(:,ncols+2) = 0.0
!      Reverse the colour table order to be more familiar to tipsy
!       k = ncols+2
	   k = 1 ! Don't reverse
       do j = 1,ncols+2
         do i = 1,3
           rgbtable(i,j) = rgbtabletmp(i,k)
         end do
!         k = k - 1
		  k = k + 1 ! Don't reverse
       end do


	  call getarg(1,pram_file)

	  if ( pram_file=="" ) then
     	write(6,'(''Input parameter file: '',$)')
      	read(5,'(a)') pram_file
	  endif
	  
! Read in parameters
	  open(unit=42,file=pram_file)
	  read(42,"(A)") dirname
	  read(42,*) inrun,istart,istop,istep
	  read(42,"(A)") outdir ! Directory to dump plot files (e.g. scratch). Don't put anything here you want to keep!
	  read(42,"(A)") movie ! movie file name, if "nomovie" then no movie
	                       ! if "gnuplot" then no movie and output as data instead of as graphics
	  select case (movie)
	  case ( "nomovie" )
	  	make_movie = .false.
	  	make_graphics = .true.
	  	print *,"Not making a movie"
	  	print *,"******************"
	  case ( "gnuplot" )
	  	make_movie = .false.
	  	make_graphics = .false.
	  	print *,"Outputting as a gnuplot format data file"
	  	print *,"do: 'set pm3d map corners2color c1' for correct x/y plot"
	  	print *,"********************************************************"
	  case default
	  	print *,"Making a movie!"
	  	print *,"***************"
	  	make_movie = .true.
	  	make_graphics = .true.
	  end select

	  read(42,*) n_pram
	  if ( n_pram>n_pram_max ) then
	  	print *,"Too many sets of parameters :("
	  	stop
	  endif
	  if ( n_pram<=0 ) then
	  	print *,"What's the point of running if you're not going to plot anything?"
	  	print *,"Gosh!"
	  	stop
	  endif

	  do i_pram=1,n_pram
		  this_plot=>all_prams(i_pram)
		  read(42,"(A)") plot_type
		  read(42,*) this_plot%dbounds(1),this_plot%dbounds(2) ! log(scalar) bounds
		  read(42,*) this_plot%rbounds(1,1),this_plot%rbounds(1,2),this_plot%rbounds(2,1),this_plot%rbounds(2,2) ! x/y bounds
		  
		  
		  this_plot%iplot_type = -1
		  this_plot%sage_flag = .false.
		  this_plot%load_star = .false.
		  this_plot%plot_name = plot_type(1:4)
		  
		  select case(this_plot%plot_name)
		  case("dens")
		  	print *,"Density plot"
		  	this_plot%iplot_type = idens_plot
		  case ("temp")
		  	print *,"Temperature plot"
		  	this_plot%iplot_type = itemp_plot
		  case ("star")
		  	print *,"Stellar plot without age cutoff"
		  	read(plot_type,*) junk, this_plot%grid_L
		  	print *,"Grid L = ",this_plot%grid_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
		  	this_plot%iplot_type = istar_plot
		  case ("sage")
		  	print *,"Stellar plot with age cutoff"
		  	read(plot_type,*) junk, this_plot%grid_L, sage_bot, sage_top
		  	print *,"Grid L = ",this_plot%grid_L
		  	print *,"Age cuts:",sage_bot,sage_top
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
		  	this_plot%iplot_type = isage_plot
		  	this_plot%sage_flag = .true.
		  	this_plot%load_star = .true.
		  case ("sagp")
		  	print *,"Plot of stellar ages (not density)"
		  	read(plot_type,*) junk, this_plot%grid_L, sage_bot, sage_top
		  	print *,"Grid L = ",this_plot%grid_L
		  	print *,"Age range:",sage_bot,sage_top
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L))
		  	this_plot%iplot_type = isagp_plot
		  	this_plot%sage_flag = .true.
		  	this_plot%load_star = .true.
		  case ("page")
		  	print *,"Straight dot plot of stellar ages"
		  	this_plot%iplot_type = ipage_plot
		  	this_plot%sage_flag = .true.
		  	this_plot%load_star = .true.
		  case ("velx")
		  	print *,"Straight dot plot of x velocity"
		  	this_plot%iplot_type = ivelx_plot
		  case ("vely")
		  	print *,"Straight dot plot of y velocity"
		  	this_plot%iplot_type = ively_plot
		  case ("velz")
		  	print *,"Straight dot plot of z velocity"
		  	this_plot%iplot_type = ivelz_plot
		  case ("vrad")
		  	print *,"Straight dot plot of radial velocity"
		  	this_plot%iplot_type = ivrad_plot
		  case ("vrot")
		  	print *,"Straight dot plot of circular velocity"
		  	this_plot%iplot_type = ivrot_plot
		  case ("curl")
		  	print *,"Binned plot of curl of velocity field"
		  	read(plot_type,*) junk, this_plot%grid_L
		  	print *,"Grid L = ",this_plot%grid_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! curl
		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! velx
		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! mass in bin
		  	allocate (this_plot%grid_val4(this_plot%grid_L,this_plot%grid_L)) ! vely
		  	this_plot%iplot_type = icurl_plot
		  case ("dspz")
		  	print *,"Binned plot of vertical velocity dispersions"
		  	read(plot_type,*) junk, this_plot%grid_L
		  	print *,"Grid L = ",this_plot%grid_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! Cumulative velocity dispersion
		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! Mean velocity
		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! Mass in bin
		  	this_plot%iplot_type = idspz_plot
		  case ("dspc")
		  	print *,"Binned plot of circular velocity dispersions"
		  	read(plot_type,*) junk, this_plot%grid_L
		  	print *,"Grid L = ",this_plot%grid_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! Cumulative velocity dispersion
		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! Mean velocity
		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! Mass in bin
		  	this_plot%iplot_type = idspc_plot
		  case ("dsc2")
		  	print *,"Plot of circular velocity dispersions"
		  	print *,"Circular velocity is binned, but individual particle velocities are plotted with a straight dot plot"
		  	read(plot_type,*) junk, this_plot%grid_L
		  	print *,"Grid L = ",this_plot%grid_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! Cumulative velocity dispersion
		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! Mean velocity
		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! Mass in bin
		  	this_plot%iplot_type = idsc2_plot
		  case ("icld")
		  	read(plot_type,*) junk,the_cloud,itime_cloudfile
			write(cloudarrayfile,"('/disc41/williams/cloudtracks/',I4.4,'/cloudarray',I4.4,'.',I5.5)") inrun,inrun,itime_cloudfile
		  	print *,"Plotting particles contained within cloud",the_cloud,"during array dump:"
		  	print *,cloudarrayfile
		  	this_plot%iplot_type = iicld_plot
		  case ("dfol")
		  	read(plot_type,"(A4,A)") junk,boundlistfile
		  	print *,"Density plot, centred on coordinates listed in"
		  	print *,boundlistfile
		  	this_plot%iplot_type = idens_plot
		  	follow_points = .true.
		  	rsize = this_plot%rbounds(:,2)-this_plot%rbounds(:,1)
		  case ("over")
		  	read(plot_type,*) junk, this_plot%grid_L,this_plot%rad_L
		  	print *,"Non-radial stellar overdensity plot"
		  	print *,"Grid L = ",this_plot%grid_L
		  	print *,"Radial L = ",this_plot%rad_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
		  	allocate (this_plot%rad_val(this_plot%rad_L))
		  	this_plot%iplot_type = iover_plot
		  case ("sfrd")
		  	print *,"Straight dot plot of star formation rate"
		  	this_plot%iplot_type = isfrd_plot
		  	this_plot%load_star = .true.
		  case ("sfrg")
		  	print *,"Star formation surface density"
		  	read(plot_type,*) junk, this_plot%grid_L
		  	print *,"Grid L = ",this_plot%grid_L
		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
		  	this_plot%iplot_type = isfrg_plot
		  	this_plot%load_star = .true.
		  case default
		    print *,"Not a valid plot option"
		    stop
		  end select

          if ( .not. allocated(this_plot%grid_val) .and. .not.make_graphics ) then
            print *,"Gnuplot output currently only supports grid data, not point data, and ",plot_type(1:4)," is a point data plot"
            stop
          endif
		  
!		  if ( plot_type(1:4)=="dens") then
!		  	print *,"Density plot"
!		  	this_plot%iplot_type = idens_plot
!		  else if ( plot_type(1:4)=="temp") then
!		  	print *,"Temperature plot"
!		  	this_plot%iplot_type = itemp_plot
!		  else if ( plot_type(1:4)=="star") then
!		  	print *,"Stellar plot without age cutoff"
!		  	read(plot_type,*) junk, this_plot%grid_L
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
!		  	this_plot%iplot_type = istar_plot
!		  else if ( plot_type(1:4)=="sage") then
!		  	print *,"Stellar plot with age cutoff"
!		  	read(plot_type,*) junk, this_plot%grid_L, sage_bot, sage_top
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	print *,"Age cuts:",sage_bot,sage_top
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
!		  	this_plot%iplot_type = isage_plot
!		  	this_plot%sage_flag = .true.
!		  else if ( plot_type(1:4)=="sagp") then
!		  	print *,"Plot of stellar ages (not density)"
!		  	read(plot_type,*) junk, this_plot%grid_L, sage_bot, sage_top
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	print *,"Age range:",sage_bot,sage_top
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L))
!		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L))
!		  	this_plot%iplot_type = isagp_plot
!		  	this_plot%sage_flag = .true.
!		  else if ( plot_type(1:4)=="page") then
!		  	print *,"Straight dot plot of stellar ages"
!		  	this_plot%iplot_type = ipage_plot
!		  	this_plot%sage_flag = .true.
!		  else if ( plot_type(1:4)=="velx") then
!		  	print *,"Straight dot plot of x velocity"
!		  	this_plot%iplot_type = ivelx_plot
!		  else if ( plot_type(1:4)=="vely") then
!		  	print *,"Straight dot plot of y velocity"
!		  	this_plot%iplot_type = ively_plot
!		  else if ( plot_type(1:4)=="velz") then
!		  	print *,"Straight dot plot of z velocity"
!		  	this_plot%iplot_type = ivelz_plot
!		  else if ( plot_type(1:4)=="vrad") then
!		  	print *,"Straight dot plot of radial velocity"
!		  	this_plot%iplot_type = ivrad_plot
!		  else if ( plot_type(1:4)=="vrot") then
!		  	print *,"Straight dot plot of circular velocity"
!		  	this_plot%iplot_type = ivrot_plot
!		  else if ( plot_type(1:4)=="curl") then
!		  	print *,"Binned plot of curl of velocity field"
!		  	read(plot_type,*) junk, this_plot%grid_L
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! curl
!		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! velx
!		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! mass in bin
!		  	allocate (this_plot%grid_val4(this_plot%grid_L,this_plot%grid_L)) ! vely
!		  	this_plot%iplot_type = icurl_plot
!		  else if ( plot_type(1:4)=="dspz") then
!		  	print *,"Binned plot of vertical velocity dispersions"
!		  	read(plot_type,*) junk, this_plot%grid_L
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! Cumulative velocity dispersion
!		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! Mean velocity
!		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! Mass in bin
!		  	this_plot%iplot_type = idspz_plot
!		  else if ( plot_type(1:4)=="dspc") then
!		  	print *,"Binned plot of circular velocity dispersions"
!		  	read(plot_type,*) junk, this_plot%grid_L
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! Cumulative velocity dispersion
!		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! Mean velocity
!		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! Mass in bin
!		  	this_plot%iplot_type = idspc_plot
!		  else if ( plot_type(1:4)=="dsc2") then
!		  	print *,"Plot of circular velocity dispersions"
!		  	print *,"Circular velocity is binned, but individual particle velocities are plotted with a straight dot plot"
!		  	read(plot_type,*) junk, this_plot%grid_L
!		  	print *,"Grid L = ",this_plot%grid_L
!		  	allocate (this_plot%grid_val(this_plot%grid_L,this_plot%grid_L)) ! Cumulative velocity dispersion
!		  	allocate (this_plot%grid_val2(this_plot%grid_L,this_plot%grid_L)) ! Mean velocity
!		  	allocate (this_plot%grid_val3(this_plot%grid_L,this_plot%grid_L)) ! Mass in bin
!		  	this_plot%iplot_type = idsc2_plot
!		  else if ( plot_type(1:4)=="icld") then
!		  	read(plot_type,*) junk,the_cloud,itime_cloudfile
!!		  	ispace = -1
!!		  	do i=1,len(plot_type)
!!				if ( plot_type(i:i)==' ' ) then
!!					if ( ispace==-1 ) then
!!						ispace = -2
!!					else if (ispace==-2 ) then
!!						ispace = i
!!					endif
!!				endif
!!		  	end do
!!			cloudarrayfile = plot_type(ispace:)
!			write(cloudarrayfile,"('/disc41/williams/cloudtracks/',I4.4,'/cloudarray',I4.4,'.',I5.5)") inrun,inrun,itime_cloudfile
!		  	print *,"Plotting particles contained within cloud",the_cloud,"during array dump:"
!		  	print *,cloudarrayfile
!		  	this_plot%iplot_type = iicld_plot
!		  else if ( plot_type(1:4)=="dfol") then
!		  	read(plot_type,"(A4,A)") junk,boundlistfile
!		  	print *,"Density plot, centred on coordinates listed in"
!		  	print *,boundlistfile
!		  	this_plot%iplot_type = idens_plot
!		  	follow_points = .true.
!		  	rsize = this_plot%rbounds(:,2)-this_plot%rbounds(:,1)
!		  endif


		  
!		  if ( this_plot%iplot_type==-1 ) then
!		  	print *,"Not a valid plot type!"
!		  	print *,"YOU BROKE IT"
!		  	!print *,"Must be dens, temp, or star"
!		  	stop
!		  endif
		  
	  end do
	  
      write(pdirname,"(I4.4,'/')") inrun
      
      write(irun_str,"(I4.4)") inrun
      
	  if (make_movie ) then  
		  command = "rm "//outdir(1:lnblnk(outdir))//"/plot"//irun_str//"*.png"
		  print *,command
		  call system(command)
		  command = "rm "//outdir(1:lnblnk(outdir))//"/montage"//irun_str//"*.png"
		  print *,command
		  call system(command)
	  endif

      iplot = 0
      
      if ( follow_points ) then
      	open(unit=52,file=boundlistfile)
      endif
      
      print *,"Start plotting!"
!      print *,istart,istop,istep
      do intime=istart,istop,istep
      	  read_itime = -1
      	  if ( follow_points ) then
			  do while (read_itime/=intime )
				read(52,"(I6,I5,5E17.8)") read_itime,idump,rdump,rdump,rdump,centrep
				print *,read_itime,intime,centrep
			  end do
			  this_plot%rbounds(:,1) = centrep - rsize/2.
			  this_plot%rbounds(:,2) = centrep + rsize/2.
			  print *,this_plot%rbounds
		  endif
		  if ( istart/=istop ) then
			print *,"Iteration:",intime
		  endif

          write(hydra_in,"(A,I4.4,A,I5.5)") "d",inrun,".",intime

	      ier=99
	      Ndum=1
	      CALL readdata(hydra_in,ier,Ndum,rm,r,v,h,dn,e,itype,sm,dftime,dnp,sfr,this_plot%load_star)
	      if (ier/=0) then
	        write(*,*) 'Input file not found',ier
	        stop
	      endif

	      CALL readdata(hydra_in,ier,nobj,rm,r,v,h,dn,e,itype,sm,dftime,dnp,sfr,this_plot%load_star)
	      call inunit

		  ! Do conversions

	      fl=float(L)-2.*padding
	      dn = dn * nunit/fl**3
		  Kunit=eunit*2.*mum/3./kb
	      e = e * Kunit*fl**2

          write(p_str,"(I6.6)") iplot
          write(*,"('Plot',I6.6)") iplot
          iplot = iplot + 1
          
          

!**********
          do i_pram=1,n_pram ! Loop through all sets of plot parameter for this datafile
	  	  	  this_plot=>all_prams(i_pram)
          	  print *,"-> Image ",i_pram
	  	  	            	  
          	  write(p_str2,"(I2.2)") i_pram
              
              if ( make_graphics ) then
                  outfile = outdir(1:lnblnk(outdir))//"/plot"//irun_str//"p"//p_str//"p"//p_str2//".png/PNG"
                  ier = pgbeg(0,outfile,1,1)
                  IF (ier/=1) STOP "Can't open png"
            ! Set the colour spectrum
                  do j = 0,ncols+1
                   call pgscr(j,rgbtable(1,j+1),rgbtable(2,j+1),rgbtable(3,j+1))
                  end do
                  call pgscr(255,1.,1.,1.)

            ! Plot the boundary
                  call pgsci(255)			  
                  roff = 0.
                  CALL PGENV(this_plot%rbounds(1,1)+roff(1),this_plot%rbounds(1,2)+roff(1),this_plot%rbounds(2,1)+roff(2),this_plot%rbounds(2,2)+roff(2),1,0)
              endif
		      
	!	      dnmax = 0.
	!	      dnmin = 1.e10
			  maxmin_set = .false.
			  
			  print *,"zeroing data"
			  
			  if ( 	this_plot%iplot_type==istar_plot .or. &
					this_plot%iplot_type==isage_plot .or. &
					this_plot%iplot_type==isagp_plot .or. &
					this_plot%iplot_type==idspz_plot .or. &
					this_plot%iplot_type==idspc_plot .or. &
					this_plot%iplot_type==icurl_plot .or. &
					this_plot%iplot_type==idsc2_plot .or. &
					this_plot%iplot_type==iover_plot .or. &
					this_plot%iplot_type==isfrg_plot) then
				  this_plot%grid_val = 0.
			  endif
			  if ( 	this_plot%iplot_type==isagp_plot .or. &
					this_plot%iplot_type==idspz_plot .or. &
					this_plot%iplot_type==idspc_plot .or. &
					this_plot%iplot_type==icurl_plot .or. &
					this_plot%iplot_type==idsc2_plot ) then
			  	  this_plot%grid_val2 = 0.
			  endif
			  if (	this_plot%iplot_type==idspz_plot .or. &
					this_plot%iplot_type==idspc_plot .or. &
					this_plot%iplot_type==icurl_plot .or. &
					this_plot%iplot_type==idsc2_plot ) then
				  this_plot%grid_val3 = 0.
			  endif
			  if (	this_plot%iplot_type==idspz_plot .or. &
					this_plot%iplot_type==idspc_plot .or. &
					this_plot%iplot_type==icurl_plot .or. &
					this_plot%iplot_type==idsc2_plot ) then
				  this_plot%grid_val4 = 0.
			  endif

			  if (  this_plot%iplot_type==iover_plot ) then
				  this_plot%rad_val = 0.d0
				  maxrad = (0.5d0-this_plot%rbounds(1,1))*sqrt(2.)
			  endif

			  if ( this_plot%iplot_type==isage_plot .or. this_plot%iplot_type==isagp_plot .or. this_plot%iplot_type==ipage_plot) then
			  	dftime = time-dftime ! Change formation time to age
			  	dftime = dftime * 1.e4 ! Convert from 10^10 years into Myr
			  	sage_max = -1.
			  	sage_min = 1.e50
			  endif
			  
	    ! Caculate averages, if necessary
			  ! Mass-weighted velocities
			  if ( this_plot%iplot_type==idspz_plot .or. &
					this_plot%iplot_type==idspc_plot .or. &
					this_plot%iplot_type==idsc2_plot .or. &
					this_plot%iplot_type==icurl_plot ) then
					
				do iobj=1,nobj
					if ( (itype(iobj)==itgas) ) then
						onplot = .true.
			  			do ii=1,2
			  				grid_ij(ii) = ((r(ii,iobj)-this_plot%rbounds(ii,1))/(this_plot%rbounds(ii,2)-this_plot%rbounds(ii,1))*this_plot%grid_L)+1
			  				if ( grid_ij(ii)<1 .or. grid_ij(ii)>this_plot%grid_L ) onplot = .false.
			  			end do
			  			if ( onplot ) then
			  				this_plot%grid_val3(grid_ij(1),grid_ij(2)) = this_plot%grid_val3(grid_ij(1),grid_ij(2)) + rm(iobj)
			  				if ( this_plot%iplot_type==idspz_plot ) then
								this_plot%grid_val2(grid_ij(1),grid_ij(2)) = this_plot%grid_val2(grid_ij(1),grid_ij(2)) + v(3,iobj)*rm(iobj)
							else if ( this_plot%iplot_type==idspc_plot  .or. &
									  this_plot%iplot_type==idsc2_plot ) then
								this_plot%grid_val2(grid_ij(1),grid_ij(2)) = this_plot%grid_val2(grid_ij(1),grid_ij(2)) + &
								((v(2,iobj)*(r(1,iobj)-0.5d0)-v(1,iobj)*(r(2,iobj)-0.5d0))/sqrt((r(1,iobj)-.5d0)**2+(r(2,iobj)-.5d0)**2))*rm(iobj)
							else if ( this_plot%iplot_type==icurl_plot ) then
								this_plot%grid_val2(grid_ij(1),grid_ij(2)) = this_plot%grid_val2(grid_ij(1),grid_ij(2)) + v(1,iobj)*rm(iobj)
								this_plot%grid_val4(grid_ij(1),grid_ij(2)) = this_plot%grid_val4(grid_ij(1),grid_ij(2)) + v(2,iobj)*rm(iobj)
							endif
			  			endif
					endif
				end do
				! Normalise out the mass
!				rmaxgrid = -1.e40
				do ii=1,this_plot%grid_L
				 do jj=1,this_plot%grid_L
				  this_plot%grid_val2(ii,jj) = this_plot%grid_val2(ii,jj)/this_plot%grid_val3(ii,jj)
				  if ( this_plot%iplot_type==icurl_plot ) then
				  	if ( this_plot%grid_val3(ii,jj)>0 ) then
					  	this_plot%grid_val4(ii,jj) = this_plot%grid_val4(ii,jj)/this_plot%grid_val3(ii,jj)
					else
						this_plot%grid_val2(ii,jj) = sqrt(-1.)
						this_plot%grid_val4(ii,jj) = sqrt(-1.)
					endif
				  endif
!				  if ( this_plot%grid_val2(ii,jj)>rmaxgrid ) then
!					rmaxgrid = this_plot%grid_val2(ii,jj)
!					!print *,ii,jj,rmaxgrid,this_plot%grid_val3(ii,jj)
!				  endif
				 end do
				end do
			  endif ! end calculate averages
			  
			  if ( this_plot%iplot_type==icurl_plot ) then ! Do differencing
			  	do ii = 2,this_plot%grid_L-1
					do jj = 2,this_plot%grid_L-1
						! curl
						if ( isNaN(this_plot%grid_val2(ii,jj+1)) .or. isNaN(this_plot%grid_val2(ii,jj-1)) .or. isNaN(this_plot%grid_val2(ii+1,jj)) .or. isNaN(this_plot%grid_val2(ii-1,jj)) ) then
							this_plot%grid_val(ii,jj) = 0.
						else
							this_plot%grid_val(ii,jj) = (this_plot%grid_val2(ii,jj+1)-this_plot%grid_val2(ii,jj-1))! + (this_plot%grid_val4(ii+1,jj)-this_plot%grid_val4(ii-1,jj))
						endif
						
						! divergence
!						if ( isNaN(this_plot%grid_val2(ii+1,jj)) .or. isNaN(this_plot%grid_val2(ii-1,jj)) .or. isNaN(this_plot%grid_val2(ii,jj+1)) .or. isNaN(this_plot%grid_val2(ii,jj-1)) ) then
!							this_plot%grid_val(ii,jj) = 0.
!						else
!							this_plot%grid_val(ii,jj) = (this_plot%grid_val2(ii+1,jj)-this_plot%grid_val2(ii-1,jj)) + (this_plot%grid_val4(ii,jj+1)-this_plot%grid_val4(ii,jj-1))
!						endif
					end do
				end do

			  endif

			  
		      if ( this_plot%iplot_type==iicld_plot ) then
				! Open the cloud array
				open(unit=43,file=cloudarrayfile)
				read(43,*) nitype
				
				biggest_n = nitype
				if (nobj>biggest_n) then
					biggest_n = nobj !whatever is biggest
				endif
				
				! Also open the appropriate itype array
				allocate(itype2(biggest_n))
		        write(itypeshortfile,"(A,I4.4,A,I5.5)") "d",inrun,".",itime_cloudfile

				itypelongfile = dirname(1:lnblnk(dirname))//'/data/'//pdirname//itypeshortfile
				
				print *,nobj,nitype,biggest_n
				
				print *,itypelongfile
				call read_itype(itypelongfile,nitype,biggest_n,itype2)
				
			  else
			  	if (.not. allocated(itype2) ) then
				  	allocate(itype2(nobj*2)) ! waste of memory but it works with the logic below
				  	! *2 to account for star formation
				endif
			  endif

			! Plot the points - first pass (background)
		     ! (n.b. can assume make_graphics is true, it can't be false for a dot plot)
			  if ( this_plot%iplot_type==iicld_plot ) then
				  do iobj=1,nobj
!				      if ( this_plot%iplot_type==iicld_plot .and. itype2(iobj)==itgas ) then
!						  read(43,*) my_cloud
!				      endif
				      
				      if ( (itype2(iobj)==itgas .and. this_plot%iplot_type==iicld_plot) ) then
			      		! Add a point for the density/temperature of this particle
			      		if ( this_plot%iplot_type==iicld_plot ) then
							this_scalar = (log10(dn(iobj)))/5.
							if (this_scalar>0.5 ) then
								this_scalar = 0.5
							endif
!							if ( my_cloud==0 ) then
!								this_scalar = .2
!							else
!								this_scalar = -0.2
!							endif
			      		endif
			      		
			      		icol = (((this_scalar)-this_plot%dbounds(1))*(ncols-2))/(this_plot%dbounds(2)-this_plot%dbounds(1))+2
			      		if ( icol>ncols ) icol=ncols
			      		if ( .not.maxmin_set ) then
			      			maxmin_set = .true.
			      			dnmax = (this_scalar)
			      			dnmin = (this_scalar)
			      		endif
			      		if ( (this_scalar)>dnmax ) dnmax = (this_scalar)
			      		if ( (this_scalar)<dnmin ) dnmin = (this_scalar)
			      		call pgsci(icol)
			      		if ( icol>=2 ) then
							if ( 	.not.this_plot%iplot_type==ipage_plot .or. &
									(this_scalar>=this_plot%dbounds(1) .and. this_scalar<=this_plot%dbounds(2))) then
					      		call pgpt1(r(1,iobj)+roff(1),r(2,iobj)+roff(2),-1)
					      	endif
					    endif
					   endif
				  end do
			  endif
			    
!		      if ( this_plot%iplot_type==iicld_plot ) then
!		        close(43)
!				! Re-open the cloud array. This is kinda inefficient, tbh
!				open(unit=43,file=cloudarrayfile)
!				read(43,*) nitype
!			  endif
		! Plot the points - second pass (foreground) (or only pass, for most things)
		     ! (n.b. can assume make_graphics is true, it can't be false for a dot plot)
			  if ( this_plot%iplot_type/=icurl_plot ) then
				  do iobj=1,nobj
					  if ( this_plot%iplot_type==iicld_plot .and. itype2(iobj)==itgas ) then
						  read(43,*) my_cloud
					  endif

					if ( 	this_plot%iplot_type==idens_plot .or. &
							this_plot%iplot_type==itemp_plot .or. &
							this_plot%iplot_type==ipage_plot .or. &
							this_plot%iplot_type==ivelx_plot .or. &
							this_plot%iplot_type==ively_plot .or. &
							this_plot%iplot_type==ivelz_plot .or. &
							this_plot%iplot_type==ivrot_plot .or. &
							this_plot%iplot_type==ivrad_plot .or. &
							this_plot%iplot_type==iicld_plot .or. &
							this_plot%iplot_type==idsc2_plot .or. &
                            this_plot%iplot_type==isfrd_plot &
							) then
						! Straight dot plots
						if ( 		(itype(iobj)==itgas .and. .not.this_plot%iplot_type==ipage_plot) &
							.or.	(itype(iobj)==itstar .and. &
                                    (this_plot%iplot_type==ipage_plot) )&
							.or.	(itype2(iobj)==itgas .and. this_plot%iplot_type==iicld_plot) ) then
							! Add a point for the density/temperature of this particle
							if ( this_plot%iplot_type==idens_plot ) then
	!			      			if ( dn(iobj)<=0. ) then
	!			      				print *,"Density of particle",iobj," is ",dn(iobj)
	!			      				print *,"That ain't proper"
	!			      				stop
	!			      			endif
								this_scalar = log10(dn(iobj))
							else if ( this_plot%iplot_type==itemp_plot ) then
								this_scalar = log10(e(iobj))
							else if ( this_plot%iplot_type==ipage_plot ) then
								this_scalar = dftime(iobj)
							else if ( this_plot%iplot_type==ivelx_plot ) then
								this_scalar = v(1,iobj)
							else if ( this_plot%iplot_type==ively_plot ) then
								this_scalar = v(2,iobj)
							else if ( this_plot%iplot_type==ivelz_plot ) then
								this_scalar = v(3,iobj)
							else if ( this_plot%iplot_type==ivrad_plot ) then
								this_scalar = (v(1,iobj)*(r(1,iobj)-.5d0)+v(2,iobj)*(r(2,iobj)-.5d0))/sqrt((r(1,iobj)-.5d0)**2+(r(2,iobj)-.5d0)**2)
							else if ( this_plot%iplot_type==ivrot_plot ) then
								this_scalar = (v(2,iobj)*(r(1,iobj)-0.5d0)-v(1,iobj)*(r(2,iobj)-0.5d0))/sqrt((r(1,iobj)-.5d0)**2+(r(2,iobj)-.5d0)**2)
							else if ( this_plot%iplot_type==iicld_plot ) then
								if ( my_cloud == the_cloud ) then
									this_scalar = 1.
								else
									this_scalar = -1.
	!							else if (my_cloud>0 ) then
	!								this_scalar = .5
	!							else
	!								this_scalar = .0
								endif
							else if ( this_plot%iplot_type==idsc2_plot ) then
								this_scalar = (v(2,iobj)*(r(1,iobj)-0.5d0)-v(1,iobj)*(r(2,iobj)-0.5d0))/sqrt((r(1,iobj)-.5d0)**2+(r(2,iobj)-.5d0)**2)
								onplot = .true.
								do ii=1,2
									grid_ij(ii) = ((r(ii,iobj)-this_plot%rbounds(ii,1))/(this_plot%rbounds(ii,2)-this_plot%rbounds(ii,1))*this_plot%grid_L)+1
									if ( grid_ij(ii)<1 .or. grid_ij(ii)>this_plot%grid_L ) onplot = .false.
								end do
								if ( onplot ) then
	!								this_scalar = this_plot%grid_val2(grid_ij(1),grid_ij(2))
									this_scalar = this_scalar - this_plot%grid_val2(grid_ij(1),grid_ij(2))
								else
									this_scalar = 0.
								endif
							else if ( this_plot%iplot_type==isfrd_plot ) then
                                if ( sfr(iobj)>0.d0 ) then
                                    this_scalar = log10(sfr(iobj))
                                else
                                    this_scalar = -10.
                                endif
							endif
							icol = (((this_scalar)-this_plot%dbounds(1))*(ncols-2))/(this_plot%dbounds(2)-this_plot%dbounds(1))+2
							if ( icol>ncols ) icol=ncols
							!if ( icol<2 ) icol = 2
		!		      		if ( this_plot%iplot_type==ipage_plot ) then
		!		      			icol = ncols-icol+2
		!		      		endif
							if ( .not.maxmin_set ) then
								maxmin_set = .true.
								dnmax = (this_scalar)
								dnmin = (this_scalar)
							endif
							if ( (this_scalar)>dnmax ) dnmax = (this_scalar)
							if ( (this_scalar)<dnmin ) dnmin = (this_scalar)
							call pgsci(icol)
							if ( icol>=2 ) then
								if ( 	.not.this_plot%iplot_type==ipage_plot .or. &
										(this_scalar>=this_plot%dbounds(1) .and. this_scalar<=this_plot%dbounds(2))) then
									call pgpt1(r(1,iobj)+roff(1),r(2,iobj)+roff(2),-1)
								endif
							endif
						endif
					else if ( &
							this_plot%iplot_type==istar_plot .or. &
							this_plot%iplot_type==isage_plot .or. &
							this_plot%iplot_type==isagp_plot .or. &
							this_plot%iplot_type==iover_plot ) then
						! Binned star plots
						if ( itype(iobj)==itstar) then
							! Star age, density etc
							if ( (.not.(this_plot%sage_flag))&
                             .or. ( (sage_bot<=dftime(iobj) .or. sage_bot<0.) .and. (sage_top>=dftime(iobj) .or. sage_top<0.) ) ) then ! Make an age cut, if you care about stellar age
								! Add this particle to the cumulative surface density grid
								onplot = .true.
								do ii=1,2
									grid_ij(ii) = ((r(ii,iobj)-this_plot%rbounds(ii,1))/(this_plot%rbounds(ii,2)-this_plot%rbounds(ii,1))*this_plot%grid_L)+1
									if ( grid_ij(ii)<1 .or. grid_ij(ii)>this_plot%grid_L ) onplot = .false.
								end do
								if ( onplot ) then
									this_plot%grid_val(grid_ij(1),grid_ij(2)) = this_plot%grid_val(grid_ij(1),grid_ij(2)) + rm(iobj)
									if ( this_plot%iplot_type==isagp_plot ) then ! age weighted plot
										this_plot%grid_val2(grid_ij(1),grid_ij(2)) = this_plot%grid_val2(grid_ij(1),grid_ij(2)) + dftime(iobj)*rm(iobj)			  				
									endif
									if ( this_plot%iplot_type==iover_plot ) then
										rad = sqrt(sum((r(1:2,iobj)-.5d0)**2))
										grid_ir = rad/maxrad*this_plot%rad_L+1
										if ( grid_ir>0 .and. grid_ir<this_plot%rad_L+1 ) then
											this_plot%rad_val(grid_ir) = this_plot%rad_val(grid_ir)+rm(iobj)
										endif
									endif
								endif
								if ( this_plot%iplot_type==isage_plot ) then
									if ( dftime(iobj)>sage_max ) then
										sage_max = dftime(iobj)
									endif
									if ( dftime(iobj)<sage_min ) then
										sage_min = dftime(iobj)
									endif
								endif
							endif
						endif
					else if ( 	this_plot%iplot_type==idspz_plot .or. &
								this_plot%iplot_type==idspc_plot .or. &
								this_plot%iplot_type==isfrg_plot ) then
						! Binned gas plots
						if ( itype(iobj)==itgas ) then
							onplot = .true.
							do ii=1,2
								grid_ij(ii) = ((r(ii,iobj)-this_plot%rbounds(ii,1))/(this_plot%rbounds(ii,2)-this_plot%rbounds(ii,1))*this_plot%grid_L)+1
								if ( grid_ij(ii)<1 .or. grid_ij(ii)>this_plot%grid_L ) onplot = .false.
							end do
							if ( onplot ) then
								if ( this_plot%iplot_type==idspz_plot ) then
									this_plot%grid_val(grid_ij(1),grid_ij(2)) = this_plot%grid_val(grid_ij(1),grid_ij(2)) + &
															rm(iobj)*(v(3,iobj)-this_plot%grid_val2(grid_ij(1),grid_ij(2)))**2
								else if ( this_plot%iplot_type==idspc_plot ) then
									this_scalar = (v(2,iobj)*(r(1,iobj)-0.5d0)-v(1,iobj)*(r(2,iobj)-0.5d0))/sqrt((r(1,iobj)-.5d0)**2+(r(2,iobj)-.5d0)**2) ! calc circ velocity
									this_scalar = this_scalar - this_plot%grid_val2(grid_ij(1),grid_ij(2)) ! subtract average
									this_scalar = this_scalar**2 ! square
									this_scalar = this_scalar*rm(iobj)
									this_plot%grid_val(grid_ij(1),grid_ij(2)) = this_plot%grid_val(grid_ij(1),grid_ij(2)) + this_scalar
								else if ( this_plot%iplot_type==isfrg_plot ) then
									this_plot%grid_val(grid_ij(1),grid_ij(2)) = this_plot%grid_val(grid_ij(1),grid_ij(2)) + sfr(iobj)
								endif
							endif

						endif
					endif
				  end do ! end do iobj
			  endif
			  
		      if ( this_plot%iplot_type==iicld_plot ) then
				close(43)
				deallocate(itype2)
		      endif
			
			
			! Grid plot!
			  if ( 	this_plot%iplot_type==istar_plot .or. &
					this_plot%iplot_type==isage_plot .or. &
					this_plot%iplot_type==isagp_plot .or. &
					this_plot%iplot_type==idspz_plot .or. &
					this_plot%iplot_type==idspc_plot .or. &
					this_plot%iplot_type==icurl_plot .or. &
					this_plot%iplot_type==iover_plot .or. &
					this_plot%iplot_type==isfrg_plot ) then
				
				! Preprocess stuff - e.g. divide by area to get surface density
			    pix_size = (this_plot%rbounds(:,2)-this_plot%rbounds(:,1))/(this_plot%grid_L)

				if ( this_plot%iplot_type==istar_plot .or. &
				     this_plot%iplot_type==isage_plot .or. &
					 this_plot%iplot_type==iover_plot) then ! density plots, with or without age cutoff
				    this_plot%grid_val = this_plot%grid_val/(pix_size(1)*pix_size(2))
				    ! Internal units are 10^10 Msun and Mpc
				    ! So mass/dist^2 = 10^10/(10^6)^2 Msun/pc^2
				    ! = 10^(-2) Msun/pc^2
				    
				    this_plot%grid_val = this_plot%grid_val*100. ! This may be incorrect? divide instead of multiply?
				else if ( this_plot%iplot_type==isagp_plot ) then ! age weighted plot
					do ix=1,this_plot%grid_L
						do iy=1,this_plot%grid_L
							if ( this_plot%grid_val(ix,iy)>0. ) then
								this_plot%grid_val(ix,iy) = this_plot%grid_val2(ix,iy)/this_plot%grid_val(ix,iy) ! Normalise.
							else
								this_plot%grid_val(ix,iy) = 0.
							endif
						end do
					end do
				else if ( 	this_plot%iplot_type==idspz_plot .or. &
							this_plot%iplot_type==idspc_plot ) then
					! RMS type plots
					this_plot%grid_val = this_plot%grid_val/this_plot%grid_val3 ! Divide by mass (i.e. normalise)
					this_plot%grid_val = sqrt(this_plot%grid_val) ! square root the array

					!this_plot%grid_val = this_plot%grid_val2 ! What's the average?
				else if ( this_plot%iplot_type==icurl_plot ) then ! difference plots
					! NB assuming a square grid!
					this_plot%grid_val = this_plot%grid_val/pix_size(1)
					print *,pix_size(1)," should be about the same as ",pix_size(2)
					print *,"Otherwise, we're all kinds of screwed"
				else if ( this_plot%iplot_type==isfrg_plot ) then
				    this_plot%grid_val = this_plot%grid_val/(pix_size(1)*pix_size(2))
					! SFR is 10^10 Msun/ 10^10 yr = Msun/yr, conveniently 
					! Distance is Mpc, want kpc
					! mass/time/dist^2 = (10^10)/(10^10)/(10^3)^2
					! = 10^(-6) Msun/yr/kpc^2
					
				    this_plot%grid_val = this_plot%grid_val/1.e6
				endif
				
				if ( this_plot%iplot_type==iover_plot ) then
				    do grid_ir = 1,this_plot%rad_L
				        rad1 = (grid_ir-1)*maxrad/this_plot%rad_L
				        rad2 = (grid_ir)*maxrad/this_plot%rad_L
				        ! ASSUMING that the plot is centred on 0.5,0.5 !!
				        ! Account for the fact that a ring of constant radius goes "out of bounds"
				        if ( rad1>(.5-this_plot%rbounds(1,1)) ) then
				            theta1 = 4.*(pi/2. - acos(rad1/(.5-this_plot%rbounds(1,1))))
				        else
				            theta1 = 2.*pi
				        endif

				        if ( rad2>(.5-this_plot%rbounds(1,1)) ) then
				            theta2 = 4.*(pi/2. - acos(rad2/(.5-this_plot%rbounds(1,1))))
				        else
				            theta2 = 2.*pi
				        endif
				        
				        this_plot%rad_val(grid_ir) = this_plot%rad_val(grid_ir)/(theta2/2*rad2**2-theta1/2*rad1**2)
				    end do
				    this_plot%rad_val = this_plot%rad_val * 100. ! Unit conversion to Msun/pc^2
				    
    				!print *,"Prior max/min dens:",maxval(this_plot%grid_val),minval(this_plot%grid_val)

                    ! Calculate overdensity
				    do ix = 1,this_plot%grid_L
    				    do iy = 1,this_plot%grid_L
    				        grid_ij(1) = ix
    				        grid_ij(2) = iy
                            r_ij(1:2) = this_plot%rbounds(1:2,1)+((grid_ij(1:2)-.5d0)*pix_size) ! Centre of this grid cell
                            rad = sqrt(sum((r_ij(1:2)-0.5d0)**2))
                    
                            grid_ir = rad/maxrad*this_plot%rad_L+1
                            if ( grid_ir>0 .and. grid_ir<this_plot%rad_L+1 ) then
                                ! Should really interpolate?
                                this_plot%grid_val(ix,iy) = this_plot%grid_val(ix,iy)/ &
                                    this_plot%rad_val(grid_ir)
!                                this_plot%grid_val(ix,iy) = this_plot%grid_val(ix,iy) - &
!                                    this_plot%rad_val(grid_ir)
!                                this_plot%grid_val(ix,iy) = this_plot%rad_val(grid_ir)
                            else
                                print *,"Out of zone!"
                                print *,r_ij
                                print *,rad,maxrad
                                print *,grid_ir
                                stop
                            endif
    				    end do
    				end do
    				!print *,"Max/min dens:",maxval(this_plot%grid_val),minval(this_plot%grid_val)
    				!print *,"Rad max/min dens:",maxval(this_plot%rad_val),minval(this_plot%rad_val)
				endif
			    
			    
			    if ( .not.make_graphics ) then
                    ! Set up gnuplot file
!                    write(outfile,"(A4,I4.4,'_gnufile.dat')") this_plot%plot_name,inrun
                    outfile = this_plot%plot_name//irun_str//"_gnufile"//p_str//".dat"
                    open(unit=57,file=outfile)
                endif
                
                print *,"Plotting gridded data"

                dnmin = this_plot%grid_val(1,1)
                dnmax = this_plot%grid_val(1,1)
                
                print *,"dnmin/max has been set"
                do ix=1,this_plot%grid_L
                    do iy=1,this_plot%grid_L
                        grid_ij(1) = ix
                        grid_ij(2) = iy
                        this_scalar = this_plot%grid_val(grid_ij(1),grid_ij(2))
    !		      		icol = (((this_scalar)-dbounds(1))*136.)/(dbounds(2)-dbounds(1))
    !		      		if ( (this_scalar)>dnmax ) dnmax = (this_scalar)
    !		      		if ( (this_scalar)<dnmin ) dnmin = (this_scalar)
                        r_ij = this_plot%rbounds(:,1)+((grid_ij(:)-1.)*pix_size)
                        if ( this_scalar>0 .or. &
                            this_plot%iplot_type==icurl_plot ) then
                            if ( this_plot%iplot_type==istar_plot .or. &
                               this_plot%iplot_type==isage_plot .or. &
                               this_plot%iplot_type==iover_plot .or. &
                               this_plot%iplot_type==isfrg_plot ) then
                                this_scalar = log10(this_scalar)
                            endif
                            !if ( (this_scalar)>dnmax .or. dnmax==0 ) dnmax = (this_scalar)
                            !if ( (this_scalar)<dnmin .or. dnmin==0 ) dnmin = (this_scalar)
                            if ( (this_scalar)>dnmax ) dnmax = (this_scalar)
                            if ( (this_scalar)<dnmin ) dnmin = (this_scalar)
                            if ( make_graphics ) then
                                icol = (((this_scalar)-this_plot%dbounds(1))*(ncols-2))/(this_plot%dbounds(2)-this_plot%dbounds(1))+2
        !			      		icol = mod(ix+iy,137)
                                if ( icol>ncols ) icol=ncols
                                if ( icol<2 ) icol = 2
                                call pgsci(icol)
        !			      		call pgpt1(r_ij(1),r_ij(2),-1)
                                call pgrect(r_ij(1)+roff(1),r_ij(1)+pix_size(1)+roff(1),r_ij(2)+roff(2),r_ij(2)+pix_size(2)+roff(2))
                            else
                                ! Gnuplot data output
                                write(57,*) r_ij(1),r_ij(2),this_scalar
                            endif
                        else
                            !write(57,*) r_ij(1),r_ij(2),0
                            write(57,*) r_ij(1),r_ij(2),this_plot%dbounds(1)
                        endif
                    end do
                    if ( .not.make_graphics ) then
                        r_ij(1) = this_plot%rbounds(1,1)+(ix-1)*pix_size(1)
                        r_ij(2) = this_plot%rbounds(2,1)+(this_plot%grid_L)*pix_size(2)
                        !write(57,*) r_ij(1),r_ij(2),0
						write(57,*) r_ij(1),r_ij(2),this_plot%dbounds(1)
                        write(57,*)
                    endif
                end do
			  	
			  	if ( .not.make_graphics ) then
                    do iy=1,this_plot%grid_L+1
                        r_ij(1) = this_plot%rbounds(1,1)+(this_plot%grid_L)*pix_size(1)
                        r_ij(2) = this_plot%rbounds(2,1)+(iy-1)*pix_size(2)
                        write(57,*) r_ij(1),r_ij(2),0
                    end do
                    write(57,*)
                    close(57)
			  	endif
			  	
			  endif ! Finish the binned plotting routine

		      print *,"plot max,min=",dnmax,dnmin
		      if ( this_plot%iplot_type==isage_plot ) then
		      	print *,"sage max,min=",sage_max,sage_min
		      	print *,"sim time=",time*1.e4
		      endif


		! Plot the label
		      if ( make_graphics ) then
                  call pgsci(255)
                  write(time_text,'("Time = ",F10.3," Myr, it = ",I5.5)') time*1.e4,intime
                  call pgtext((this_plot%rbounds(1,1)+this_plot%rbounds(1,2))/2.+roff(1),this_plot%rbounds(2,2)+(this_plot%rbounds(2,2)-this_plot%rbounds(2,1))/20.+roff(2),time_text)      
        !		  call pgtext((this_plot%rbounds(1,1)+this_plot%rbounds(1,2))/2.,(this_plot%rbounds(2,1)+this_plot%rbounds(2,2))/2.,'TESTSADKFJDHASKLJHFDSKLFJHSDKLFJHSDLKFJHkjhsdflkjhads')      


            ! Plot axis
        !		  call pgsls(1)
                  do jj=2,ncols
                    call pgsci(jj)
        !		  	call pgslw(200)
                    line_y(1) = this_plot%rbounds(2,1)+((jj-2)*(this_plot%rbounds(2,2)-this_plot%rbounds(2,1)))/(ncols-1)
                    line_y(2) = this_plot%rbounds(2,1)+((jj-1)*(this_plot%rbounds(2,2)-this_plot%rbounds(2,1)))/(ncols-1)
                    line_x(1) = this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*0.95
                    line_x(2) = this_plot%rbounds(1,1)+(this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*1.
        !			line_x = this_plot%rbounds(1,:)
        !			line_y = this_plot%rbounds(2,:)+.02
                    call pgrect(line_x(1)+roff(1),line_x(2)+roff(1),line_y(1)+roff(2),line_y(2)+roff(2))
                  end do
                  write(time_text,'(F8.2)') this_plot%dbounds(1)
                  call pgtext(this_plot%rbounds(1,2)+roff(1),this_plot%rbounds(2,1)+roff(2),time_text)
                  write(time_text,'(F8.2)') this_plot%dbounds(2)
                  call pgtext(this_plot%rbounds(1,2)+roff(1),this_plot%rbounds(2,2)+roff(2),time_text)
            
    !			  roff(1) = this_plot%rbounds(1,2) + (this_plot%rbounds(1,2)-this_plot%rbounds(1,1))*0.3
    !			  roff(2) = 0.	  		
            
                  CALL PGEND
              endif
			  print *,"________________"
	  	end do
	  	
	  	! Do montage
	  end do
	  
	  if ( follow_points ) then
      	close(52)
      endif
	  
	  print *,"Images complete!"
	  
	  if ( make_movie ) then		
		  if ( n_pram>1 ) then
		  
		  	  print *,"Performing montages"
			  ! Do montages
			  iplot = 0
	          do intime=istart,istop,istep
	          	command = "montage"
				write(p_str,"(I6.6)") iplot
	          	do i_pram =1,n_pram
					write(p_str2,"(I2.2)") i_pram
	          		command = trim(command)//" "//outdir(1:lnblnk(outdir))//"/plot"//irun_str//"p"//p_str//"p"//p_str2//".png"
	          	end do
	          	command = trim(command)//" -geometry +0+0"
	          	command = trim(command)//" "//outdir(1:lnblnk(outdir))//"/montage"//irun_str//"p"//p_str//".png"
	          	write(*,"(A)") command
	          	call system(command)
	          	iplot = iplot + 1
	          end do
	          print *,"Making movie"
			  command = "~williams/ffmpeg/ffmpeg -r 5 -i "//outdir(1:lnblnk(outdir))//"/montage"//irun_str//"p%6d.png -y -qscale 0 "//movie
!			  command = "~williams/ffmpeg-0.5/ffmpeg -sameq -r 20 -i "//outdir(1:lnblnk(outdir))//"/montage"//irun_str//"p%6d.png -y "//movie
			  write(*,"(A)") command
			  call system(command)
      	  else  
			  print *,"Making movie"
      	  	! Just animate the one plot for all itime
			  command = "~williams/ffmpeg/ffmpeg -r 5 -i "//outdir(1:lnblnk(outdir))//"/plot"//irun_str//"p%6dp01.png -y -qscale 0 "//movie
!			  command = "~williams/ffmpeg-0.5/ffmpeg -sameq -r 20 -i "//outdir(1:lnblnk(outdir))//"/plot"//irun_str//"p%6dp01.png -y "//movie
			  write(*,"(A)") command
			  call system(command)
		  endif

	  endif
	  
	  end program hquickplot
      
