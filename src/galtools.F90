module galtools
    
    real(kind=8), dimension(3),save :: com_r
    
    real(kind=8), save :: xzgrad,yzgrad
    
    contains
    
    subroutine sort_yourself_out
        implicit none
        logical, dimension(4) :: typemask
        
        typemask = [.true.,.true.,.false.,.true.]
!        typemask = [.false.,.true.,.false.,.false.]
        print *,"Centering disc and rotating to xy plane"

        call find_com(typemask)
        call centre_com

        call find_xgrad(typemask)
        call xrot_gal
        call find_ygrad(typemask)
        call yrot_gal
        
    end subroutine sort_yourself_out
    
    subroutine centre_disc_only
        implicit none
        logical, dimension(4) :: typemask
        
        typemask = [.true.,.true.,.false.,.true.]
        print *,"Centering disc only"

        call find_com(typemask)
        call centre_com
    
    end subroutine centre_disc_only

    subroutine centre_disc_only_id(maxid)
        implicit none
        integer, intent(in) :: maxid
        logical, dimension(4) :: typemask
        
        typemask = [.true.,.true.,.false.,.true.]
        print *,"Centering disc only, p<",maxid
!        typemask = [.false.,.true.,.false.,.false.]

        call find_com_id(typemask,maxid)
        call centre_com
    end subroutine centre_disc_only_id
    
    
    subroutine find_com(type_flag)
        implicit none
        logical, dimension(4), intent(inout) :: type_flag

        integer :: max_id
        
        max_id = -1
        
        call find_com_id(type_flag,max_id)
    end subroutine find_com
    
! Find the centre of mass of all particles of the given types, but don't do anything
    subroutine find_com_id(type_flag,max_id)
        use data_vals
        implicit none
        
        ! type_flag(1) include star
        ! type_flag(2) gas
        ! type_flag(3) dark
        ! type_flag(4) feed
        ! i.e. include if type_flag(itype(ip)+1) is true
        logical, dimension(4), intent(in) :: type_flag
        
        integer :: max_id
        
        real(kind=8) :: mtot,weight
        
        integer :: ip
        
        integer :: nmatch
        
        mtot = 0.d0
        com_r = 0.d0
        
        nmatch = 0
        
        do ip=1,np
            if ( itype(ip)>2 .or. itype(ip)<0 ) then
                print *,"Bad itype",ip,itype(ip)
                stop
            endif
            if ( type_flag(itype(ip)+1) .and. (max_id==-1 .or. id_p(ip)<=max_id) ) then
                nmatch = nmatch +1
                weight = log10(sca_data(ip,2)* 6.77d-26/(sca_data(ip,1)*1.d12))+28
                !weight = log10(sca_data(ip,2)* 6.77d-26/(sca_data(ip,1)*1.d12))+24
                if ( weight>0. ) then
                    com_r = com_r + r_p(1:3,ip) * sca_data(ip,1)*weight
                    mtot = mtot + sca_data(ip,1)*weight
                endif
            endif
        end do
        
        print *,"nmatch=",nmatch,max_id
        
        com_r = com_r/mtot
        
    end subroutine find_com_id

! Translate all particles to move the previously calculated centre of mass to 0,0,0
    subroutine centre_com
        use data_vals
        implicit none
        
        integer :: ip
        
        do ip=1,np
            r_p(1:3,ip) = r_p(1:3,ip) - com_r
        end do
        
    end subroutine centre_com

! Find the x/z gradient of the given types of particles, so we can rotate the disc later
    subroutine find_xgrad(type_flag) ! find the gradient of the disc (i.e. how much the disc is rotated versus the flat plane)
        use data_vals
        implicit none
        
        ! type_flag(1) include star
        ! type_flag(2) gas
        ! type_flag(3) dark
        ! type_flag(4) feed
        ! i.e. include if type_flag(itype(ip)+1) is true
        logical, dimension(4), intent(in) :: type_flag
        
        real(kind=8) :: x2sum, xzsum
        
        
        real(kind=8) :: rad, weight
        real(kind=8), parameter :: rad_cut = .02, rad_min=.0005 ! in 100 kpc

        integer :: ip
        
        x2sum = 0.d0
        xzsum = 0.d0
        
        do ip=1,np
            if ( type_flag(itype(ip)+1) ) then
                rad = sqrt(sum(r_p(1:3,ip)**2))
                if ( rad<rad_cut .and. rad>rad_min ) then
                    !weight = log10(sca_data(ip,2)* 6.77d-26/(sca_data(ip,1)*1.d12))+28
                    weight = 1.
                    if ( weight>0 ) then
                        x2sum = x2sum + r_p(1,ip)**2 * weight
                        xzsum = xzsum + r_p(1,ip)*r_p(3,ip) * weight
                    endif
                endif
            endif
        end do
        
        xzgrad = xzsum/x2sum
    end subroutine find_xgrad

! Find the y/z gradient of the given types of particles, so we can rotate the disc later    
    subroutine find_ygrad(type_flag) ! find the gradient of the disc (i.e. how much the disc is rotated versus the flat plane)
        use data_vals
        implicit none
        
        ! type_flag(1) include star
        ! type_flag(2) gas
        ! type_flag(3) dark
        ! type_flag(4) feed!!! x not true in dotty_plotty, feedback is still gas I think?
        ! i.e. include if type_flag(itype(ip)+1) is true
        logical, dimension(4), intent(in) :: type_flag
        
        real(kind=8) :: y2sum,yzsum
        
        
        real(kind=8) :: rad, weight
        real(kind=8), parameter :: rad_cut = .02, rad_min=.0005 ! in 100 kpc
        
        integer :: ip
        
        y2sum = 0.d0
        yzsum = 0.d0
        do ip=1,np
            if ( type_flag(itype(ip)+1) ) then
                rad = sqrt(sum(r_p(1:3,ip)**2))
                if ( rad<rad_cut .and. rad>rad_min ) then
                    !weight = log10(sca_data(ip,2)* 6.77d-26/(sca_data(ip,1)*1.d12))+28
                    weight = 1.
                    if ( weight>0 ) then
                        y2sum = y2sum + r_p(2,ip)**2 * weight
                        yzsum = yzsum + r_p(2,ip)*r_p(3,ip) * weight
                    endif
                endif
            endif
        end do
        
        yzgrad = yzsum/y2sum
    end subroutine find_ygrad

! Rotate the disc according to the previously calculating x/z gradient, aiming to make this gradient flat
    subroutine xrot_gal
        use data_vals
        
        real(kind=8) :: xzgrad_c ! sqrt(1-grad**2), c for "complement" in a weird funky way
        real(kind=8) :: oldx,oldz
        integer :: ip
        
        xzgrad_c = sqrt(1.d0-xzgrad**2)
        
        do ip=1,np
            oldx = r_p(1,ip)
            oldz = r_p(3,ip)

            r_p(1,ip) = oldx*xzgrad_c + oldz*xzgrad
            r_p(3,ip) = -oldx*xzgrad + oldz*xzgrad_c
        end do
        
    end subroutine xrot_gal

! Rotate the disc according to the previously calculating y/z gradient, aiming to make this gradient flat
    subroutine yrot_gal
        use data_vals
        
        real(kind=8) :: yzgrad_c ! sqrt(1-grad**2), c for "complement" in a weird funky way
        real(kind=8) :: oldy,oldz
        integer :: ip
        
        yzgrad_c = sqrt(1.d0-yzgrad**2)
        
        do ip=1,np
            oldy = r_p(2,ip)
            oldz = r_p(3,ip)

            r_p(2,ip) = oldy*yzgrad_c + oldz*yzgrad
            r_p(3,ip) = -oldy*yzgrad + oldz*yzgrad_c
        end do
        
    end subroutine yrot_gal

end module galtools
