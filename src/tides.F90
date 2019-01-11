module tides
    real(kind=8), save :: last_TM_tot, latest_time
    real(kind=8), save :: gp_x0, gp_y0, gp_z0
    real(kind=8), save :: gp_x, gp_y, gp_z
    real(kind=8), save :: gp_vx,gp_vy,gp_vz
    real(kind=8), save :: gp_a
    real(kind=8), save :: gp_nfw_const, gp_nfw_pot_const
    real(kind=8), save :: gp_startt
    
    real(kind=8), parameter :: TMUGYR=0.471d0, MUSM=1.0d12, LUPC=1.0d5, VUKMS=207.4d0, G=1.0d0

!    real(kind=8), parameter :: min_dt = 1.e-4/TMUGYR ! 0.1 Myr
    real(kind=8), parameter :: min_dt = 1.e-5/TMUGYR ! 0.01 Myr
    
    real(kind=8), save :: big_mass,dwarf_mass
    
    real(kind=8), save :: tidal_tot

    contains
    subroutine tidal_init
        use plot_prams, only: dirname, run_name
        implicit none

        real(kind=8) :: gp_m, gp_c
        real(kind=8) :: gp_rho0 ! Never need this again!

        real(kind=8), parameter :: conc_const = 0.358 ! 200 * 4/3 * pi * rho_crit in internal units
        real(kind=8), parameter :: rho0_const = 0.0286 ! 200 /3 * rho_crit in internal units
        real(kind=8), parameter :: fourpi = 12.5663706144
        
        character(len=128) :: filename
        
        character(len=128), dimension(4) :: tidelines
        character(len=128) :: inline
        
        integer :: ierr
        integer :: ip

        integer :: idum
        logical :: is_pot
        
        return ! just read in stuff from bbvals file now
        
        filename = trim(dirname)//"/"//trim(run_name)//"/diskev/ini/input.dat"
        
        print *,"Getting potential from",trim(filename)
        
        open(50,file=trim(filename),status='old')
        do idum=1,4
            read(50,"(A)",iostat=ierr) tidelines(idum)
        end do
        do while (ierr==0 )
            read(50,"(A)",iostat=ierr) inline
            if ( ierr==0 ) then
                do idum=1,3
                    tidelines(idum) = tidelines(idum+1)
                end do
                tidelines(4) = inline
            endif
        end do
        close(50)
        
        print *,"reading pot"
        ! Galaxy potential stuff
        read(tidelines(1),*) gp_m, gp_c ! Mass and concentration parameter of halo in Msun, pc
        is_pot = .true.

        read(tidelines(2),*) gp_startt ! Time (in Gyr) until external halo activates (for equilibration)
        read(tidelines(3),*) gp_x0,gp_y0,gp_z0 ! Initial position of centre of halo (pc)
        read(tidelines(4),*) gp_vx,gp_vy,gp_vz ! (constant) velocity of passing halo, km/s
        
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

        print *,gp_nfw_pot_const,gp_nfw_const,gp_rho0,gp_a
    end subroutine tidal_init

    subroutine tidal_update!(t)
        use data_vals, only: time
        implicit none
        real(kind=8) :: fmag,ax,ay,az,dt,deltat,t
        
        return ! just read in stuff from bbvals file - no need to update
        
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
