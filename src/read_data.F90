subroutine read_data(intime)
    use plot_prams
    implicit none
    integer, intent(in) :: intime
    integer :: idump

	integer :: ii
    character*256 :: filename
    
    ! Open the file giving the IDs for the dumps and find which one we are
    write(filename,"(A,'/',A,'/diskev/output/ana/ostep.dat')") trim(dirname),trim(run_name)

    print *,"Reading:",trim(filename)

    open(unit=50,file=filename,status='old')
    do ii=1,intime
    	read(50,*) idump
    end do
    close(50)
	
    call read_gcd(idump,dirname,run_name)
    
    
    montage_files = ""
end subroutine read_data



    ! GCD+ data is packaged according to its initial processor during the run
    ! i.e. if you ran on 128 cores, then the particles are split into 128 chunks
    ! Here we loop over the number of processors, and repackage the chunks of particles into our big scalar and vector arrays
subroutine read_gcd(idump,dirname,run_name)
    use data_vals
    use plot_prams, only: xydisc, tidedisc, do_selection, select_file, centdisconly
    use galtools, only: sort_yourself_out,centre_disc_only_id
    use tides, only: tidal_update,tidal_init,gp_x,gp_y,gp_z,gp_vx,gp_vy,gp_vz
    implicit none
    integer, intent(in) :: idump
    character*256, intent(in) :: dirname,run_name
    
    integer :: ip
    integer :: ii
    
    character*256 :: filename
    
    integer, allocatable, dimension(:) :: ng_proc,ns_proc,ndm_proc,nb_proc
    
    integer :: np_max
	integer :: ndm
	integer :: nb

    integer :: iproc,nprocs
    
    integer :: ival
    integer :: pstart, pend
	real(kind=8), allocatable, dimension(:) :: invals
	integer, allocatable, dimension(:) :: dummy	
    integer :: idummy
    real(kind=8) :: rdummy

    real(kind=8) :: rad

    integer :: nselect
    integer, dimension(:), allocatable :: select_ids
    integer, dimension(:), allocatable :: select_pn

    !integer, allocatable, dimension(:) :: id_p

    integer :: ierr
    
    logical, parameter :: read_photo = .true. ! should nicely skip over things if there are no photo files

    ! Read in selection data, if applictable

    if ( do_selection ) then
        open(unit=50,file=select_file,status='old')
        read(50,*) nselect
        allocate(select_ids(nselect))
        allocate(select_pn(nselect))
        do ii=1,nselect
            read(50,*) select_ids(ii)
        end do
        close(50)
    endif

    ! Read in data from GCD+ format
    write(filename,"(A,'/',A,'/diskev/output/data/bbvals',I6.6,'n0000')") trim(dirname),trim(run_name),idump
    
    print *,"reading file",trim(filename)
	open(unit=50,file=filename,status='old',form="unformatted")

	read(50) nb,ndm,idummy,rdummy,time
	np = nb + ndm
	print *,"nb,ndm,np = ",nb,ndm,np
	
    if ( tidedisc ) then
    	read(50) nprocs,idummy,idummy,idummy,&
    	         gp_vx,gp_vy,gp_vz,&
    	         gp_x,gp_y,gp_z
	else
	    read(50) nprocs
	endif

	allocate(ng_proc(nprocs))
	allocate(ns_proc(nprocs))
	allocate(nb_proc(nprocs))
	allocate(ndm_proc(nprocs))

	do iproc=1,nprocs
		read(50) ng_proc(iproc),ndm_proc(iproc),ns_proc(iproc)
	end do

	nb_proc = ng_proc + ns_proc ! Total number of baryon particles from each processor
	
	if ( sum(nb_proc)+sum(ndm_proc)/=np ) then
		print *,"np wrong :("
		print *,sum(nb_proc),sum(ndm_proc),np
		nb = sum(nb_proc)
		ndm = sum(ndm_proc)
        np = nb + ndm
	    print *,"forcing: nb,ndm,np = ",nb,ndm,np
	endif
	
	np_max = max(maxval(nb_proc),maxval(ndm_proc))
	
	allocate(dummy(np))
	allocate(invals(np_max))

	nvec = 1
	nsca = 18 + 14 ! 18 intrinsic values, 15 derived values
	if ( read_photo ) then
	    nsca = nsca + 2 ! lums and fluxes
	endif
	if ( .not. allocated(sca_codes) ) then
		allocate(sca_codes(nsca))
		allocate(sca_units(nsca))
        ! TO DO:
        ! Order these so that we only have the codes and data for the values we actually read
		sca_codes(1) = "m"
        sca_units(1) = "Msun"
		sca_codes(2) = "rho"
        sca_units(2) = "g/cm**3"
		sca_codes(3) = "u_p"
        sca_units(3) = "erg/g"
		sca_codes(4) = "h_p"
        sca_units(4) = "pc"
		sca_codes(5) = "divv"
        sca_units(5) = "" ! Place-holders for later
		sca_codes(6) = "alpv"
        sca_units(6) = ""
		sca_codes(7) = "alpu"
        sca_units(7) = ""
		sca_codes(8) = "myu"
        sca_units(8) = ""
		sca_codes(9) = "ZHe"
        sca_units(9) = "mass fraction"
		sca_codes(10) = "ZZ"
        sca_units(10) = "mass fraction"
		sca_codes(11) = "ZC"
        sca_units(11) = "mass fraction"
		sca_codes(12) = "ZN"
        sca_units(12) = "mass fraction"
		sca_codes(13) = "ZO"
        sca_units(13) = "mass fraction"
		sca_codes(14) = "ZNe"
        sca_units(14) = "mass fraction"
		sca_codes(15) = "ZMg"
        sca_units(15) = "mass fraction"
		sca_codes(16) = "ZSi"
        sca_units(16) = "mass fraction"
		sca_codes(17) = "ZFe"
        sca_units(17) = "mass fraction"
		sca_codes(18) = "vsq"
        sca_units(18) = "km/s"
		
		
		sca_codes(19) = "T"
        sca_units(19) = "K"
		sca_codes(20) = "cs"
        sca_units(20) = "km/s"
		sca_codes(21) = "difc"
        sca_units(21) = "g/cm/s"
		sca_codes(22) = "difv"
        sca_units(22) = "g/cm/s"
		sca_codes(23) = "difb"
        sca_units(23) = "g/cm/s"
        sca_codes(24) = "vrad" ! 3d radial velocity
        sca_units(24) = "km/s"
        sca_codes(25) = "dftv"
        sca_units(25) = "yr"
        sca_codes(26) = "dftc"
        sca_units(26) = "yr"
!         sca_codes(27) = "tcol"
!         sca_units(27) = "yr"
        sca_codes(27) = "p_p"
        sca_units(27) = "dyne/cm**2"
        sca_codes(28) = "FeO"
        sca_units(28) = "mass fraction"
        sca_codes(29) = "vr2d" ! 2d radial velocity
        sca_units(29) = "km/s"
        sca_codes(30) = "vc" ! 2d circular velocity
        sca_units(30) = "km/s"
        sca_codes(31) = "nH" ! H number density
        sca_units(31) = "/cm^3"
        sca_codes(32) = "H" ! H number, normalised for smoothing
        sca_units(32) = "N (kpc^3/cm^3)"
        
        if ( read_photo ) then
            sca_codes(33) = "lum"
            sca_units(33) = "erg/s"
            sca_codes(34) = "flx"
            sca_units(34) = "erg/s/cm^2"
        endif

		allocate(vec_codes(nvec))
		vec_codes(1) = "vel"
	endif

    if ( allocated(sca_data) ) then
    	deallocate(sca_data)
    	deallocate(vec_data)
    	deallocate(r_p)
    	deallocate(itype)
    	deallocate(ifeedtype)
    	deallocate(id_p)
    endif

    allocate(sca_data(np,nsca))
	allocate(vec_data(3,np,nvec))
    allocate(r_p(3,np))
    allocate(itype(np))
    allocate(ifeedtype(np))

    allocate(id_p(np))

	print *,"Reading integers!"
    	
	! Skip over all integer data for now - list_ap etc
	do ival=1,4
		pstart = 1
		do iproc=1,nprocs
			read(50) dummy(1:nb_proc(iproc))
            if ( ival==1) then !  .and. do_selection ! always read id_p
                pend = pstart+nb_proc(iproc)-1
                id_p(pstart:pend) = dummy(1:nb_proc(iproc))
                pstart = pend+1
            endif
			if ( ival==2 ) then
				pend = pstart+nb_proc(iproc)-1
				itype(pstart:pend) = dummy(1:nb_proc(iproc))
				pstart = pend+1
			endif
			if ( ival==4 ) then
				pend = pstart+nb_proc(iproc)-1
				ifeedtype(pstart:pend) = dummy(1:nb_proc(iproc))
				pstart = pend+1
			endif
		end do
	end do
	
    ! Convert itype from GCD+ type into my type
    
    do ip=1,np
        if ( itype(ip)<=0 ) then
            itype(ip) = itgas
        else
            itype(ip) = itstar
        endif
    end do
    
	print *,"Reading reals!"

	! Load in position, velocity, mass, density, u_p
    do ival=1,9
        pstart = 1
        do iproc=1,nprocs
            read(50) invals(1:nb_proc(iproc))
            pend = pstart+nb_proc(iproc)-1
            select case(ival)
                case (1:3) ! Position
                    r_p(ival,pstart:pend) = invals(1:nb_proc(iproc))
                case (4:6) ! Velocity
                    vec_data(ival-3,pstart:pend,1) = invals(1:nb_proc(iproc))
                case (7:9) !mass,rho,u_p, potentially all other scalars too
                    sca_data(pstart:pend,ival-6) = invals(1:nb_proc(iproc))
            end select
            pstart = pend+1
        end do
	end do
	
	close(50)
    

    ! Read in more hydro stuff
    write(filename,"(A,'/',A,'/diskev/output/data/bbhyds',I6.6,'n0000')") trim(dirname),trim(run_name),idump
    print *,"reading file",trim(filename)
    open(unit=50,file=filename,status='old',form="unformatted")
    read(50)
    read(50)
    do ival=1,5
        pstart = 1
        do iproc=1,nprocs
            read(50) invals(1:nb_proc(iproc))
            pend = pstart+nb_proc(iproc)-1
            sca_data(pstart:pend,ival+3) = invals(1:nb_proc(iproc))
            pstart = pend+1
        end do
    end do
    close(50)
    
    
    ! Read in metals
    write(filename,"(A,'/',A,'/diskev/output/data/bbmets',I6.6,'n0000')") trim(dirname),trim(run_name),idump
    print *,"reading file",trim(filename)
    open(unit=50,file=filename,status='old',form="unformatted")
    read(50)
    read(50)
    do ival=1,19
        pstart = 1
        do iproc=1,nprocs
            read(50,iostat=ierr) invals(1:nb_proc(iproc))
            if ( ierr/=0 ) exit
            pend = pstart+nb_proc(iproc)-1
            select case(ival)
                case(1:9) ! Metals
                    sca_data(pstart:pend,ival+8) = invals(1:nb_proc(iproc))
!                case(10:18) ! Initial metals, ignore                
                case(19) ! Velocity dispersion
                    sca_data(pstart:pend,18) = invals(1:nb_proc(iproc))
            end select
            pstart = pend+1
        end do
        if ( ierr/=0 ) exit
    end do
    close(50)
    
    ! read in luminosity & stuff
    if ( read_photo ) then
        write(filename,"(A,'/',A,'/diskev/output/data/bblums',I6.6,'n0000')") trim(dirname),trim(run_name),idump
    
        open(unit=50,file=filename,status='old',form="unformatted",iostat=ierr)
        if ( ierr==0 ) then ! nicely skip  this section if the file doesn't exist, which is fine
            read(50)
            read(50)

            do ival=1,2
                pstart = 1
                do iproc=1,nprocs
                    read(50,iostat=ierr) invals(1:nb_proc(iproc))
                    if ( ierr/=0 ) exit
                    pend = pstart+nb_proc(iproc)-1                    
                    sca_data(pstart:pend,ival+32) = invals(1:nb_proc(iproc))

                    pstart = pend+1
                end do
                if ( ierr/=0 ) exit
            end do
    
            close(50)
        else
            sca_data(1:np,33:34) = 0.d0 ! just zero everything if there's no photoheating file
        endif
    endif

    if ( ndm>0 ) then
    	! Read in dark matter 
    	
		write(filename,"(A,'/',A,'/diskev/output/data/bdvals',I6.6,'n0000')") trim(dirname),trim(run_name),idump
		print *,"reading file",trim(filename)
		open(unit=50,file=filename,status='old',form="unformatted")
		read(50)
		read(50)
		do ii=1,nprocs
			read(50)
		end do
		
		
		! Skip over all DM integer data for now - list_ap etc
		do ival=1,2
			do iproc=1,nprocs
				read(50) dummy(1:ndm_proc(iproc))
			end do
		end do
		
		itype(nb+1:np) = itdark
		
		! But load in the reals data
		! Load in position, velocity, mass, density, h
		do ival=1,9
			pstart = nb+1
			do iproc=1,nprocs
				read(50) invals(1:ndm_proc(iproc))
				pend = pstart+ndm_proc(iproc)-1
				select case(ival)
					case (1:3) ! Position
						r_p(ival,pstart:pend) = invals(1:ndm_proc(iproc))
					case (4:6) ! Velocity
						vec_data(ival-3,pstart:pend,1) = invals(1:ndm_proc(iproc))
					case (7:8) !mass,rho
						sca_data(pstart:pend,ival-6) = invals(1:ndm_proc(iproc))
					case (9) !h, because we don't have u_p
						sca_data(pstart:pend,4) = invals(1:ndm_proc(iproc))
				end select
				pstart = pend+1
			end do
		end do
		close(50)
    endif
	
	deallocate(invals)
	deallocate(dummy)
	deallocate(ng_proc)
	deallocate(nb_proc)
	deallocate(ns_proc)
	
    
    ! Rotate disc to xy plane, if requested
    if ( xydisc ) then
        call sort_yourself_out
    endif
    if ( centdisconly ) then
        call centre_disc_only_id(499999)
        print *,"ONLY COUNTING IDs LESS THAN 500000"
    endif

    ! Derived data
    call derived_data

    if ( tidedisc ) then
        call tidal_init
        call tidal_update
        r_p(1,:) = r_p(1,:) - gp_x*100. ! tides are in 100 kpc units, r_p is in kpc units
        r_p(2,:) = r_p(2,:) - gp_y*100.
        r_p(3,:) = r_p(3,:) - gp_z*100.
    endif

    if ( do_selection ) then

        ! Make the selection (AFTER rotating etc)

        print *,"Taking selection slice"
        
        print *,"Zeroing"

!        do ip=1,np
!            do ii=1,nselect
!                if ( id_p(ip)==select_ids(ii) ) then
!                    select_pn(ii) = ip
!                end if
!            end do
!        end do
!        
        print *,"Finding same things"
        call match_ids(select_ids,np,nselect,select_pn)

!         if ( any(select_pn==-1) ) then
!             print *,"Not all selection found!"
!             !print *,select_pn
!             !print *,select_ids
!             print *,count(select_pn==-1),nselect
!             stop
!         endif

        print *,"Copying everything across"
        !np = nselect
        ip = 0
        do while (ip<nselect)
            if ( select_pn(ip)==-1 ) then
                select_pn(ip:nselect-1) = select_pn(ip+1:nselect)
                nselect = nselect - 1
            else
                ip = ip + 1
            endif
        end do
        
        r_p(:,1:nselect) = r_p(:,select_pn)
        itype(1:nselect) = itype(select_pn)
        ifeedtype(1:nselect) = ifeedtype(select_pn)
        sca_data(1:nselect,:) = sca_data(select_pn,:)
        vec_data(:,1:nselect,:) = vec_data(:,select_pn,:)
        

        !print *,"SELECTION:",select_pn
        !print *,"SELECTION_r:",r_p(:,1:nselect)
        !print *,"SELECTION_itype",itype(1:nselect)
        !print *,"SELECTION_rho",sca_data(1:nselect,2)
        
!        print *,"TOTAL METALS:",sum(sca_data(1:np,9))
!        print *,"TOTAL MASS:",sum(sca_data(1:np,1))
!        print *,"zero metals:",count(sca_data(1:np,9)<1.e-6)
!        print *,"GAS P",count(itype(1:np)==itgas)
!        print *,"STAR P",count(itype(1:np)==itstar)
!        print *,"DARK P",count(itype(1:np)==itdark)
!        print *,"FEED P",count(itype(1:np)==itfeed)
        
        deallocate(select_pn)
        deallocate(select_ids)
    endif

end subroutine read_gcd

subroutine match_ids(select_ids,np,nselect,select_pn)
    use data_vals, only: itdark, itype,id_p
    implicit none
    integer, intent(in) :: np, nselect
    !integer, intent(in), dimension(np) :: id_p
    integer, intent(in), dimension(nselect) :: select_ids
    
    integer, intent(out), dimension(nselect) :: select_pn
    
    integer :: ip
    
    integer :: ihigh,ilow,iguess
    logical :: done_idloop

    select_pn = -1
    
    do ip=1,np
!        do ii=1,nselect
!            if ( id_p(ip)==select_ids(ii) ) then
!                select_pn(ii) = ip
!            end if
!        end do

        if ( itype(ip)/=itdark ) then
            ! Do a binary search
            ilow = 1
            ihigh = nselect
            iguess = (ihigh+ilow)/2
            done_idloop = .false.
            do while ( .not.done_idloop )
                if ( id_p(ip)==select_ids(iguess) ) then
                    done_idloop = .true.
                    select_pn(iguess) = ip
                else
                    if ( id_p(ip)>select_ids(iguess) ) then
                        ilow = iguess+1
                        iguess = (ihigh+ilow)/2
                    else
                        ihigh = iguess-1
                        iguess = (ihigh+ilow)/2
                    endif
                
                    if ( ilow>ihigh ) then
                        done_idloop = .true.
                    endif
                endif
            end do
        end if
        
    end do
    !print *,count(select_pn==-1),nselect
    
    return
end subroutine match_ids

subroutine derived_data
    use data_vals
    use gcdp_cool, only: get_cool
    implicit none

    integer :: ip
    real(kind=8) :: rad

    call get_cool

    ! This magic number combines the unit conversion and all the constants (including gamma) for calculating the temperature
    ! T = (constant) * u_p * myu_p
    sca_data(1:np,19) = 3472070.96307 * sca_data(1:np,3) * sca_data(1:np,8)
    ! cs = sqrt(10/9 * u_p)
    sca_data(1:np,20) = dsqrt(10./9. * sca_data(1:np,3) )*207.4 ! km/s

    ! diffusion coefficient (cm^2/s)
        ! sound speed
    sca_data(1:np,21) = sca_data(1:np,1) * sca_data(1:np,20) * sca_data(1:np,4) * 2089.1 ! convert 10^12 Msun/(100 kpc)^3 * km/s * (100 kpc) to g/cm /s
        ! vsq
    sca_data(1:np,22) = sca_data(1:np,1) * sca_data(1:np,18)*207.4 * sca_data(1:np,4) * 2089.1 ! convert 10^12 Msun/(100 kpc)^3 * km/s * (100 kpc) to g/cm /s
        ! vsq *and* sound speed
    sca_data(1:np,23) = sca_data(1:np,1) * sqrt(sca_data(1:np,20)**2+(sca_data(1:np,18)*207.4)**2) * sca_data(1:np,4) * 2089.1 ! convert 10^12 Msun/(100 kpc)^3 * km/s * (100 kpc) to g/cm /s
    ! diffusion time-scales
        ! vsq 
    do ip=1,np
        if ( sca_data(ip,22)>0. ) then
            sca_data(ip,25) = sca_data(ip,4)**2*sca_data(ip,1)/sca_data(ip,22)*2.04273195e14 ! convert ((100 kpc)**2*(10^12 Msun)/(100 kpc)^3)/(g/cm/s) to years
        endif
        if ( sca_data(ip,23)>0. ) then
        ! vsq *and* sound speed
            sca_data(ip,26) = sca_data(ip,4)**2*sca_data(ip,1)/sca_data(ip,23)*2.04273195e14 ! convert ((100 kpc)**2*(10^12 Msun)/(100 kpc)^3)/(g/cm/s) to years
        endif
    end do

    ! hydrogen density
    sca_data(1:np,31) = ((sca_data(1:np,1)-((sca_data(1:np,10)+sca_data(1:np,9))/1e12)) &
            /sca_data(1:np,1))*sca_data(1:np,2)*(6.77d-26/1.67265e-24) ! amu/cm^3 basically
    
    ! hydrogen fraction, weighted to give N /cm^3 in sm3d
    ! 1e12 solar masses to amu
    ! we eventually want amu/cm^3, not amu/kpc^3, so multiply by (1 cm/1 kpc)**3
    sca_data(1:np,32) = ((sca_data(1:np,1)-((sca_data(1:np,10)+sca_data(1:np,9))/1e12)))*(40771.425)


    ! Unit conversions - squish these all into a subroutine later
    
    do ip=1,np
        sca_data(ip,9:17) = (sca_data(ip,9:17))/(sca_data(ip,1)*1e12) ! Metallicities in mass fractions
    end do
    sca_data(1:np,18) = sca_data(1:np,18) * 207.4 ! velocity dispersion in km/s

    ! Ratio of Fe/O
    do ip=1,np
        if ( sca_data(ip,13)>0. ) then
            sca_data(ip,28) = sca_data(ip,17)/sca_data(ip,13)
        endif
    end do
    ! versus Solar ratio
    sca_data(1:np,28) = sca_data(1:np,28)*8.19658119658

    r_p(:,1:np) = r_p(:,1:np)*100. ! positions in *kilo*pc

    ! Radial velocities, in km/s
    do ip=1,np
        ! 3d rad
        rad = sqrt(sum(r_p(1:3,ip)**2))
        sca_data(ip,24) = sum(r_p(1:3,ip)*vec_data(1:3,ip,1))/rad * 207.4
        ! 2d rad
        rad = sqrt(sum(r_p(1:2,ip)**2))
        sca_data(ip,29) = sum(r_p(1:2,ip)*vec_data(1:2,ip,1))/rad * 207.4
        ! 2d circ
        sca_data(ip,30) = -(r_p(1,ip)*vec_data(2,ip,1)-r_p(2,ip)*vec_data(1,ip,1))/rad * 207.4
    end do

    ! P=nkT=rho/myu kT, kT=(3/2) u*myu
    ! P=(3/2)rho*u
    ! constant includes conversion factors from internal units into dyne/cm**2
    sca_data(1:np,27) = 4.37015735d-11*sca_data(1:np,2)*sca_data(1:np,3)

    
    sca_data(1:np,1) = sca_data(1:np,1) * 1.e12 ! Mass in Msun
    sca_data(1:np,2) = sca_data(1:np,2) * 6.77d-26! Density in g/cm^3
    sca_data(1:np,3) = sca_data(1:np,3) * 4.30345382e14! u_p in erg/g
    sca_data(1:np,4) = sca_data(1:np,4) * 100. ! smoothing length in kpc
    vec_data(1:3,1:np,1) = vec_data(1:3,1:np,1) * 207.4 ! velocity in km/s

    
!    sca_data(1:np,33) = sca_data(1:np,33) * 5.7526882d43 ! lum in erg/s
    sca_data(1:np,33) = sca_data(1:np,33) * 15027921107.6 ! lum in solar luminosities
    sca_data(1:np,34) = sca_data(1:np,34) * 0.00060405848d0 ! flux in erg/s/cm^2
    sca_data(1:np,34) = sca_data(1:np,34) / 0.0015859021 ! flux in Habing units
end subroutine derived_data
