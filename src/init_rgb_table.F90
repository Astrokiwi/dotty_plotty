subroutine init_rgb_table
    use plot_prams, only: ncols,rgbtable
    implicit none
    
    integer :: i,j,k
    real(kind=4) :: ff
    
    
    integer, dimension(125*3) :: rgb2
    data rgb2 / 	   0,  0,  0, &
                 0,  0,  13, &
                 0,  0,  26, &
                 0,  0,  38, &
                 0,  0,  51, &
                 0,  0,  64, &
                 0,  0,  77, &
                 0,  0,  89, &
                 0,  0,  102, &
                 0,  0,  115, &
                 0,  0,  128, &
                 0,  0,  140, &
                 0,  0,  153, &
                 0,  0,  166, &
                 0,  0,  179, &
                 0,  0,  191, &
                 0,  0,  204, &
                 0,  0,  217, &
                 0,  0,  230, &
                 0,  0,  242, &
                 0,  0,  255, &
                 0,  0,  255, &
                 13,  0,  255, &
                 26,  0,  255, &
                 38,  0,  255, &
                 51,  0,  255, &
                 64,  0,  255, &
                 77,  0,  255, &
                 89,  0,  255, &
                 102,  0,  255, &
                 115,  0,  255, &
                 128,  0,  255, &
                 140,  0,  255, &
                 153,  0,  255, &
                 166,  0,  255, &
                 179,  0,  255, &
                 191,  0,  255, &
                 204,  0,  255, &
                 217,  0,  255, &
                 230,  0,  255, &
                 242,  0,  255, &
                 255,  0,  255, &
                 255,  0,  243, &
                 255,  0,  231, &
                 255,  0,  219, &
                 255,  0,  206, &
                 255,  0,  194, &
                 255,  0,  182, &
                 255,  0,  170, &
                 255,  0,  158, &
                 255,  0,  146, &
                 255,  0,  134, &
                 255,  0,  121, &
                 255,  0,  109, &
                 255,  0,  97, &
                 255,  0,  85, &
                 255,  0,  73, &
                 255,  0,  61, &
                 255,  0,  49, &
                 255,  0,  36, &
                 255,  0,  24, &
                 255,  0,  12, &
                 255,  0,  0, &
                 255,  0,  0, &
                 255,  6,  0, &
                 255,  13,  0, &
                 255,  19,  0, &
                 255,  26,  0, &
                 255,  32,  0, &
                 255,  38,  0, &
                 255,  45,  0, &
                 255,  51,  0, &
                 255,  57,  0, &
                 255,  64,  0, &
                 255,  70,  0, &
                 255,  77,  0, &
                 255,  83,  0, &
                 255,  89,  0, &
                 255,  96,  0, &
                 255,  102,  0, &
                 255,  108,  0, &
                 255,  115,  0, &
                 255,  121,  0, &
                 255,  128,  0, &
                 255,  134,  0, &
                 255,  140,  0, &
                 255,  147,  0, &
                 255,  153,  0, &
                 255,  159,  0, &
                 255,  166,  0, &
                 255,  172,  0, &
                 255,  179,  0, &
                 255,  185,  0, &
                 255,  191,  0, &
                 255,  198,  0, &
                 255,  204,  0, &
                 255,  210,  0, &
                 255,  217,  0, &
                 255,  223,  0, &
                 255,  230,  0, &
                 255,  236,  0, &
                 255,  242,  0, &
                 255,  249,  0, &
                 255,  255,  0, &
                 255,  255,  12, &
                 255,  255,  24, &
                 255,  255,  36, &
                 255,  255,  49, &
                 255,  255,  61, &
                 255,  255,  73, &
                 255,  255,  85, &
                 255,  255,  97, &
                 255,  255,  109, &
                 255,  255,  121, &
                 255,  255,  134, &
                 255,  255,  146, &
                 255,  255,  158, &
                 255,  255,  170, &
                 255,  255,  182, &
                 255,  255,  194, &
                 255,  255,  206, &
                 255,  255,  219, &
                 255,  255,  231, &
                 255,  255,  243, &
                 255,  255,  255 /
    ncols = 125
    allocate(rgbtable(9,3,ncols+3))
    
    ! Tipsy style rainbow colour scheme    
    k = 0
    rgbtable(1,:,1) = 0.0
    do j = 2,ncols+1
     do i = 1,3
       k = k + 1
       rgbtable(1,i,j) = float(rgb2(k))/256.0
     end do
    end do

    rgbtable(1,:,ncols+2) = 1.0
    rgbtable(1,:,ncols+3) = 0.9

	! "Hot" colour scheme, better for printing
!    rgbtable(:,:) = 0.0	
!    rgbtable(:,1) = 1.0
!    do j=2,ncols+1
!        ff = real(j-1)/(ncols)
!        rgbtable(:,j) = [1.-ff**2,1.-sqrt(ff),0.]
!    end do
!    rgbtable(:,ncols+2) = 0.0	


    ! "Black Hot" colour scheme, for calendar
    rgbtable(2,:,:) = 1.0	
    rgbtable(2,:,1) = 0.0
    do j=2,ncols+1
        ff = real(j-1)/(ncols)
        if ( ff<.33 ) then
            rgbtable(2,:,j) = [ff/.66,0.,0.]
        else if ( ff<.66 ) then
            rgbtable(2,:,j) = [ff/.66,(ff-.33)/.68,0.]
        else
            rgbtable(2,:,j) = [1.,(ff-.33)/.68,2.94*(ff-.66)]
        endif
    end do
    rgbtable(2,:,ncols+2) = 1.0
    rgbtable(2,:,ncols+3) = 0.3

    ! Black & white colour scheme, for printing
    rgbtable(3,:,:) = 0.0 
    rgbtable(3,:,1) = 1.0
    do j=2,ncols+1
        ff = 1-real(j-1)/(ncols)
        rgbtable(3,:,j) = [ff,ff,ff]
    end do
    rgbtable(3,:,ncols+2) = 0.0
    rgbtable(3,:,ncols+3) = 0.3

    ! Double rainbow for velocities etc

    rgbtable(4,:,:) = 1.0 
    rgbtable(4,:,1) = 0.0

    do j=2,ncols+1
        ff = real(j-1)/(ncols)
        if ( ff<.5 ) then
            rgbtable(4,:,j) = [(1.-ff*1.8),0.,0.1]
        else
            rgbtable(4,:,j) = [0.,(ff-.5)*1.8+.1,0.1]
        endif
    end do
    rgbtable(4,:,ncols+2) = 1.0
    rgbtable(4,:,ncols+3) = 0.0

    ! Double rainbow for infall - negative has a shallower slope

    rgbtable(5,:,:) = .999
    rgbtable(5,:,1) = 0.0

    do j=2,ncols+1
        ff = real(j-1)/(ncols)
        if ( ff<.5 ) then
            rgbtable(5,:,j) = [(1.-ff*1.4),0.,0.5]
        else
            rgbtable(5,:,j) = [0.,(ff-.5)*1.4+.3,0.5]
        endif
    end do
    rgbtable(5,:,ncols+2) = .99
    rgbtable(5,:,ncols+3) = 0.01

    ! Don`t fade to black

    rgbtable(6,:,:) = .999
    rgbtable(6,:,1) = 0.0

    do j=2,ncols+1
        ff = real(j-1)/(ncols)
        rgbtable(6,:,j) = [ff,abs(ff-.5)*2.,1.-ff]
    end do
    rgbtable(6,:,ncols+2) = .99
    rgbtable(6,:,ncols+3) = 0.01

    ! "Black Hot" colour scheme, with better contrast
    rgbtable(7,:,:) = 1.0 
    rgbtable(7,:,1) = 0.0

    do j=2,ncols+1
        ff = real(j-1)/(ncols)
        if ( ff<.33 ) then
            rgbtable(7,:,j) = [ff/.33,0.,0.]
        else if ( ff<.66 ) then
            rgbtable(7,:,j) = [1.,(ff-.33)*3,2.*(ff-.5)]
        else
            rgbtable(7,:,j) = [1.,1.,2.*(ff-.5)]
        endif
    end do

    do j=2,ncols+1
        do i=1,3
            rgbtable(7,i,j) = max(min(1.,rgbtable(7,i,j)),0.)
        end do
    end do
    rgbtable(7,:,ncols+2) = 1.0
    rgbtable(7,:,ncols+3) = 0.3

    ! Don`t fade to black, show small deviation from zero

    rgbtable(8,:,:) = .999
    rgbtable(8,:,1) = 0.0

    do j=2,ncols+1
    if ( j<=4 ) then
        rgbtable(8,:,j) = [.5,.5,.5]
    else
        ff = real(j-1)/(ncols)
        rgbtable(8,:,j) = [ff,abs(ff-.5)*2.,1.-ff]
    endif
    end do
	rgbtable(8,:,ncols+2) = .99
    rgbtable(8,:,ncols+3) = 0.01


    ! "blue-white" colour scheme, aiming for more contrast
    rgbtable(9,:,:) = 1.0 
    rgbtable(9,:,1) = 0.0

    do j=2,ncols+1
        ff = real(j-1)/(ncols)
        if ( ff<.33 ) then
            rgbtable(9,:,j) = [0.,0.,ff/.33]
        else if ( ff<.66 ) then
            rgbtable(9,:,j) = [(ff-.33)*3,2.*(ff-.5),1.]
        else
            rgbtable(9,:,j) = [1.,2.*(ff-.5),1.]
        endif
    end do

    ! force to interval [0,1]
    do j=2,ncols+1
        do i=1,3
            rgbtable(9,i,j) = max(min(1.,rgbtable(9,i,j)),0.)
        end do
    end do
    rgbtable(9,:,ncols+2) = 1.0
    rgbtable(9,:,ncols+3) = 0.3

end subroutine init_rgb_table
