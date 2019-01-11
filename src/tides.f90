module tides
    real(kind=8), save :: last_TM_tot, latest_time
    real(kind=8), save :: gp_x0, gp_y0, gp_z0
    real(kind=8), save :: gp_x, gp_y, gp_z
    real(kind=8), save :: gp_vx,gp_vy,gp_vz
    real(kind=8), save :: gp_a
    real(kind=8), save :: gp_nfw_const, gp_nfw_pot_const
    real(kind=8), save :: gp_startt
    
    real(kind=8), parameter :: TMUGYR=0.471d0, MUSM=1.0d12, LUPC=1.0d5, VUKMS=207.4d0, G=1.0d0

    real(kind=8), parameter :: min_dt = 1.e-4/TMUGYR ! 0.1 Myr
    
    real(kind=8), save :: big_mass,dwarf_mass
    
    real(kind=8), save :: tidal_tot

    contains
    subroutine tidal_init
        use gcd_data, only: dirname, run_name

        real(kind=8) :: gp_m, gp_c
        real(kind=8) :: gp_rho0 ! Never need this again!

        real(kind=8), parameter :: conc_const = 0.358 ! 200 * 4/3 * pi * rho_crit in internal units
        real(kind=8), parameter :: rho0_const = 0.0286 ! 200 /3 * rho_crit in internal units
        real(kind=8), parameter :: fourpi = 12.5663706144
        
        character(len=128) :: filename
        
        integer :: ierr
        integer :: ip

        integer :: idum
        logical :: is_pot
        
        filename = trim(dirname)//"/"//trim(run_name)//"/diskev/ini/input.dat"
        
        print *,"Getting potential from",trim(filename)
        
        open(50,file=trim(filename),status='old')
        do idum=1,10
            read(50,*,iostat=ierr)
        end do
        
        
        print *,"reading pot"
        ! Galaxy potential stuff
        read(50,*) gp_m, gp_c ! Mass and concentration parameter of halo in Msun, pc
        is_pot = .true.

        read(50,*) gp_startt ! Time (in Gyr) until external halo activates (for equilibration)
        read(50,*) gp_x0,gp_y0,gp_z0 ! Initial position of centre of halo (pc)
        read(50,*) gp_vx,gp_vy,gp_vz ! (constant) velocity of passing halo, km/s
        close(50)
        
        big_mass = gp_m
        dwarf_mass = -100

        ! Unit conversions!

        gp_m = gp_m / MUSM

        gp_x0 = gp_x0 / LUPC
        gp_y0 = gp_y0 / LUPC
        gp_z0 = gp_z0 / LUPC

        gp_x = gp_x0
        gp_y = gp_y0
        gp_z = gp_z0

        gp_vx = gp_vx / VUKMS
        gp_vy = gp_vy / VUKMS
        gp_vz = gp_vz / VUKMS

        gp_a = (gp_m / conc_const)**(1./3.) / gp_c

        gp_rho0 = rho0_const * gp_c**3/(log(1+gp_c)-gp_c/(1.+gp_c))

        gp_nfw_const = -fourpi * G * gp_rho0 * gp_a
        gp_nfw_pot_const = fourpi * G * gp_rho0 * gp_a**2
        last_TM_tot = gp_startt/TMUGYR
        latest_time = gp_startt/TMUGYR
        tidal_tot = 0.
    end subroutine tidal_init

    subroutine tidal_update(t)
        use gcd_data
        implicit none
        real(kind=8) :: fmag,ax,ay,az,dt,deltat,t
        
        deltat = time - last_TM_tot
        dt = min(deltat,min_dt)
    
        t = last_TM_tot
    
        do while ( t<time)
            ! Calculate acceleration

!            dist = sqrt((gp_x)**2+(gp_y)**2+(gp_z)**2)
!            normdist = dist/gp_a
!
!            fmag = gp_nfw_const * (normdist/(1.+normdist)-log(1.+normdist))/normdist**2
!
!            fmag = fmag/dist
!
!            ax = fmag * (gp_x)
!            ay = fmag * (gp_y)
!            az = fmag * (gp_z)

            call nfw_f(gp_x,gp_y,gp_z,ax,ay,az,fmag)
    
            ! Get force on the centre of the dwarf galaxy
            ! In frame where dwarf galaxy is stationary,
            ! the effective acceleration of the massive galaxy
            ! is the negative of this.

            ! Leapfrog!
            gp_x = gp_x + dt * gp_vx
            gp_y = gp_y + dt * gp_vy
            gp_z = gp_z + dt * gp_vz

            gp_vx = gp_vx - dt*ax
            gp_vy = gp_vy - dt*ay
            gp_vz = gp_vz - dt*az
        
            t = t + dt
        end do
            
        latest_time = t
        last_TM_tot = t
        print *,gp_x,gp_y,gp_z,"POT LOC"

    end subroutine tidal_update

    subroutine tidal_rad(dm_mass)
        use gcd_data
        implicit none

        real(kind=8), intent(in) :: dm_mass
        real(kind=8) :: dist,normdist,fmag,ax,ay,az,dt,deltat,t
        integer :: ip
        real(kind=8) :: r_j,v_esc,v_dwarf
        real(kind=8), parameter :: third = 1.d0/3.d0
        
        if ( dwarf_mass<0.d0 ) then
            dwarf_mass = 0.d0
            dwarf_mass = (sum(SCA_DATA(:,SCA_M))+dm_mass)!*MUSM ! already in solar masses
        endif

        if ( time>latest_time+min_dt ) then
            call tidal_update(t)
            call nfw_f(gp_x,gp_y,gp_z,ax,ay,az,fmag)
            dist = sqrt((gp_x)**2+(gp_y)**2+(gp_z)**2)
            v_esc = sqrt(2.*fmag*dist) * VUKMS ! this is incorrect for an NFW profile
            v_dwarf = sqrt(gp_vx**2+gp_vy**2+gp_vz**2) * VUKMS

            dist = dist*LUPC
            r_j = (dwarf_mass/(3.*big_mass))**third * dist

            write(12,"(6E12.5)",advance='no') t*4.71d8,dist,r_j,v_esc,fmag,v_dwarf
            call escape_frac(r_j,v_esc)
            call tide_work(deltat)
        endif

    end subroutine tidal_rad
    
    subroutine escape_frac(r_j,v_esc)
        use gcd_data
        implicit none
        real(kind=8) :: r_j,v_esc
        real(kind=8) :: m_tot,m_tide,m_esc
        real(kind=8) :: v_norm,rad,m
        real(kind=8), dimension(3) :: vrel_p
        
        integer :: ip
        
        m_tot = 0.d0
        m_tide = 0.d0
        m_esc = 0.d0
        
        do ip=1,np
            if ( itype(ip)==itgas .or. itype(ip)==itfeed ) then
                m = SCA_DATA(ip,SCA_M)
                m_tot = m_tot + m
                
                rad = sqrt(sum(r_p(ip,:)**2))*1.e2 !*LUPC - already in kpc
                
                if ( rad>=r_j ) then
                    m_tide = m_tide + m
                    
                    vrel_p = [gp_vx,gp_vy,gp_vz]*VUKMS ! Velocity of big galaxy relative to "stationary" dwarf galaxy
                    vrel_p = vec_data(ip,VEC_VEL,:)-vrel_p! Sum of these velocities
                    v_norm = sqrt(sum(vrel_p(:)**2))
                    
                    if ( v_norm>=v_esc ) then
                        m_esc = m_esc + m
                    endif
                endif
                
            endif
        end do
        
        write(12,"(2E12.5)",advance='no') m_tide/m_tot,m_esc/m_tot
        
    end subroutine escape_frac
    
    subroutine tide_work(deltat)
        implicit none
        real(kind=8), intent(in) :: deltat
        
        real(kind=8) :: dist,tdist
        real(kind=8) :: f1,f0
        real(kind=8) :: vdir
        
        dist = sqrt(gp_x**2 + gp_y**2 + gp_z**2)
        
        tdist = dist + 10./LUPC
        f1 = nfw_mag(tdist)
        
        tdist = dist - 10./LUPC
        f0 = nfw_mag(tdist)
        
        tidal_tot = tidal_tot + (f0 - f1) * deltat
        
        write(12,"(E12.5)") tidal_tot
        
    end subroutine tide_work

    
    subroutine nfw_particle_pot(x,y,z,pot)
        implicit none
        real(kind=8), intent(in) :: x,y,z
        real(kind=8), intent(inout) :: pot

        real(kind=8) :: dist


        dist = sqrt((x-gp_x)**2+(y-gp_y)**2+(z-gp_z)**2)
        call nfw_pot(dist,pot)

    end subroutine nfw_particle_pot

    
    subroutine nfw_pot(d,pot)
        implicit none
        real(kind=8), intent(in) :: d
        real(kind=8), intent(inout) :: pot

        real(kind=8) :: normdist

        normdist = d/gp_a
        pot = gp_nfw_pot_const * log(1+normdist)/normdist
    end subroutine nfw_pot
    
    subroutine nfw_f(x,y,z,ax,ay,az,fm)
        implicit none
        real(kind=8), intent(in) :: x,y,z
        real(kind=8), intent(out) :: ax,ay,az,fm
        
        real(kind=8) :: dist,fmag
        
        
        dist = sqrt((x)**2+(y)**2+(z)**2)
!        normdist = dist/gp_a

!        fmag = gp_nfw_const * (normdist/(1.+normdist)-log(1.+normdist))/normdist**2
        fmag = nfw_mag(dist)
        fm = fmag

        fmag = fmag/dist

        ax = fmag * (x)
        ay = fmag * (y)
        az = fmag * (z)
    end subroutine nfw_f
    
    function nfw_mag(d)
        implicit none
        real(kind=8), intent(in) :: d
        real(kind=8) :: nfw_mag
        
        real(kind=8) :: normdist

        normdist = d/gp_a

        nfw_mag = gp_nfw_const * (normdist/(1.+normdist)-log(1.+normdist))/normdist**2
        
        return
    end function nfw_mag

end module tides
