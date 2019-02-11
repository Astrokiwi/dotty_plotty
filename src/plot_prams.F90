module plot_prams
    ! Plot stuff
    type plot_pram_t ! Parameters for a plot
      integer :: iplot_type
      integer :: iover_type
      logical :: dolog
      logical :: doabs
      logical :: massnorm
      logical :: dodsub
      logical :: smoo3d
      character*4 :: plot_var
      real :: dbounds(2)
      real :: rbounds(2,2)
      real :: zbounds(2)
      real :: phi,theta
      real :: smooslice
      integer :: lx,ly ! Size of grid
      integer :: iptype
      character*256 :: plot_name
      !real :: dsub
      integer :: isca
      integer :: palette
    end type plot_pram_t
    
    real :: phi_min,phi_max
    real :: theta_min,theta_max
    real :: g_phi,g_theta ! global
    integer :: nrot
    logical :: rot_only
    logical :: xydisc ! rotate disc to xy plane
    logical :: centdisconly ! focus on disc centre but don't rotate disc
    logical :: tidedisc ! disc from frame of external galaxy

    integer, parameter :: n_runs_max = 24
    integer,save :: n_runs
    character(len=64), save, dimension(n_runs_max) :: runlist

!    integer, save :: palette_type ! Which RGB palette to select (see init_rgb_table.F90)
                                  ! Negative palette = dump ascii table
    !integer, save :: res_x,res_y ! Resolution of each plot
    real(kind=4), save :: plot_inches ! Size of the plot in inches (because that`s how PGPLOT rolls)

    
    integer, parameter :: n_pram_max = 64
    integer,save :: n_pram
    
    type (plot_pram_t), target, save, dimension(n_pram_max) :: all_prams
    
    integer, parameter ::  idot_plot=1,&
                          igrid_plot=2,&
                          iradmean_plot=3,&
                          izmean_plot=4
                          
                          
    integer, parameter ::	imax_plot=1,& ! plot max in cell
                        	imin_plot=2,& ! plot min in cell
                        	iup_plot=3,& ! closest in "z" direction overwrites
                          isum_plot=4,& ! Just sum all dots in cell
                          ifar_plot=5,& ! Plot "biggest" number (i.e. furthest from zero)
                          isurf_plot=6,& ! Sum dots & divide by cell size
                          ismooth_plot=7 ! SPH kernelize the dots & divide by cell size
                          
    
    ! I/O stuff
    integer, save :: istart,istop,istep ! Iterations to start at, stop at, and end at
    character*256,save :: dirname
    character*256,save :: run_name
    
    character*1048,save :: montage_files
!    character*1048 :: anim_files
    character*64,save :: pid_str
    character*256,save :: anim_list_file ! text file giving list of filenames of animation frames
    character*256,save :: anim_file_base ! base-name of actual animation output (.gif file) (run name is appended to this)

    integer,save :: anim_type
    character*1048,save :: anim_files

    

    logical ,save:: do_selection
    character*1048,save :: select_file

    ! rgb_table stuff
    integer, save :: ncols
    real, save, allocatable, dimension(:,:,:) :: rgbtable

    ! Options!
    
    logical, parameter :: clean_plot = .false.
    logical, parameter :: cleanish_plot = .false.
end module plot_prams
