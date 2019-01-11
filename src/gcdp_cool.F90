

module gcdp_cool
      integer,parameter :: MNT_CRTB=100,MNMET_CRTB=10,MNNH_CRTB=20,MNZ_CRTB=50

      integer nz_crtb,nnh_crtb,nmet_crtb,nt_crtb
      double precision z_crtb(0:MNT_CRTB),nh_crtb(0:MNNH_CRTB),met_crtb(0:MNMET_CRTB),t_crtb(0:MNT_CRTB)
      double precision, dimension(0:MNT_CRTB,0:MNMET_CRTB,0:MNNH_CRTB,0:MNZ_CRTB) :: cr_crtb,hr_crtb,myu_crtb


    contains
        subroutine read_cool
            implicit none

            integer i,j,k
            integer iz,inh,imet,it,nval,nc
            integer lmet10,dmet10,lnh10,dnh10
            integer z100,met10,nh10
            character czred*5,cmet*5,cnh*5,filen*128
            double precision crad,hrad

            double precision SI_zeor

            print *,"Reading cooling data"

            open(50,file='cool/setcool.dat',status='old')
            read(50,'(I6)') nz_crtb
            read(50,"(3I6)") nnh_crtb,lnh10,dnh10
            read(50,"(3I6)") nmet_crtb,lmet10,dmet10
            read(50,'(I6)') nt_crtb
            read(50,'(F5.2)') SI_zeor
            close(50)

            open(51,file='cool/setcoolzred.dat',status='old')
            do iz=0,nz_crtb-1
              read(51,*) z100
              z_crtb(iz)=dble(z100)/100.0d0
              write(czred,'(a1,i4.4)') 'z',z100
              do inh=0,nnh_crtb-1
                nh10=lnh10+inh*dnh10
                nh_crtb(inh)=dble(nh10)/10.0d0
                if(nh10.ge.0) then
                  write(cnh,'(a2,i3.3)') 'hp',nh10
                else
                  write(cnh,'(a2,i3.3)') 'hm',-nh10
                endif
                do imet=0,nmet_crtb-1
                  met10=lmet10+imet*dmet10
                  met_crtb(imet)=dble(met10)/10.0d0
                  if(met10.ge.0) then
                    write(cmet,'(a2,i3.3)') 'Zp',met10
                  else
                    write(cmet,'(a2,i3.3)') 'Zm',-met10
                  endif
            ! *** set filename ***
                  write(filen,'(a,a5,a5,a5,a4)') 'cool/chrate',czred,cmet,cnh,'.dat'
                  open(50,file=filen,status='old')
                  do i=1,4
                    read(50,*)
                  enddo
                  do i=0,nt_crtb-1
                    read(50,"(f5.3,f8.3,f8.3,f6.3)") t_crtb(i),cr_crtb(i,imet,inh,iz) &
                                      ,hr_crtb(i,imet,inh,iz),myu_crtb(i,imet,inh,iz)
                  enddo
                  close(50)
                enddo
              enddo
            enddo
            close(51)
            print *,"Read cooling data"

        end subroutine read_cool

    subroutine get_cool
      use data_vals
      implicit none
      double precision, parameter :: CLIMIT=1.0d-3, TUK=1.0d4, MUSM=1.0d12, DU=6.77d-26, MP=1.67265e-24, XZSOL=0.019d0, ERU=1.96d-27
      double precision, parameter :: K_MU=4.3e10, KCGS=1.381e-16, MYU=0.6d0, TPRHO=((KCGS/(MYU*MP))/K_MU)

      DOUBLE PRECISION current_z

      integer i
      double precision cradp,hradp,radp
      double precision lowmet
      double precision nhp,lognhp,logTp,metp,logmetp,logTlimK
      double precision dnh,dmet,dt
      double precision wt1,wt2,wd1,wd2,wz1,wz2,dz
      double precision wm1,wm2
      integer pn,it,id,iz,im

      integer, parameter :: SI_flagrad = 0

! lowest metallicity
      lowmet=10.0d0**met_crtb(0)
! temperature limit at K
      logTlimK=dlog10(CLIMIT*TUK)
! set delta nH, dmet, dT
      dnh=nh_crtb(2)-nh_crtb(1)
      dmet=met_crtb(2)-met_crtb(1)
      dt=t_crtb(2)-t_crtb(1)
!

! find redshift 
      if(current_z.lt.z_crtb(nz_crtb-1)) then
        do i=1,nz_crtb-1
          if(z_crtb(i).gt.current_z) then
            iz = i
            goto 70
          endif
        enddo
        iz = nz_crtb-1
   70   dz = z_crtb(iz)-z_crtb(iz-1)
! /*** weight for redshift ***/
        wz1 = (z_crtb(iz)-current_z)/dz
        wz2 = (current_z-z_crtb(iz-1))/dz
      else
        iz = nz_crtb-1
        wz1 = 0.0d0
        wz2 = 1.0d0
      endif
      
!ocl novrec
!cdir nodep
      do i =0,np-1
        pn = i + 1
        if ( itype(pn)==itgas ) then
! *** nH ***
        nhp = ((sca_data(pn,1)-((sca_data(pn,10)+sca_data(pn,9))/MUSM))/sca_data(pn,1))*sca_data(pn,2)*(DU/MP)
        Lognhp=dlog10(nhp)
! *** find metallicity ***
        metp = (sca_data(pn,9))/(sca_data(pn,1)*MUSM*XZSOL)
        if(metp.gt.lowmet) then
          logmetp=dlog10(metp)
        else
          logmetp=met_crtb(0)
        endif
        im = int((logmetp-met_crtb(0))/dmet)

        if(im.lt.0) then
! *** use the vale at met_crtb(0) ***
          im = 0
          wm1=1.0d0
          wm2=0.0d0
        else if(im.ge.nmet_crtb-1) then
! *** use the vale at met_crtb(nmet_crtb-1) ***
          im = nmet_crtb-2
          wm1=0.0d0
          wm2=1.0d0     
        else
          wm1 = (met_crtb(im+1)-logmetp)/dmet
          wm2 = (logmetp-met_crtb(im))/dmet
        endif
        id=int((lognhp-nh_crtb(0))/dnh)
        if(id.lt.0) then
! *** use the vale at nh_crtb(0) ***
          id = 0
        else if(id.ge.nnh_crtb-1) then
! *** use the vale at met_crtb(nmet_crtb-1) ***
          id = nnh_crtb-2
        endif
        wd1 = (nh_crtb(id+1)-lognhp)/dnh
        wd2 = (lognhp-nh_crtb(id))/dnh
! /*** find density, allow extrapolation ***/
        if(SI_flagrad.gt.-2) then
! *** before EoR use cooling rate lognh=0.0***
          id = int((0.0d0-nh_crtb(0))/dnh)
          wd1=(nh_crtb(id+1)-0.0)/dnh
          wd2=(0.0d0-nh_crtb(id))/dnh
        endif
! /*** temperature ***/
        logTp = dlog10(3472070.96307 * sca_data(pn,3) * sca_data(pn,8))
        if(logTp.lt.logTlimK) then
          logTp=logTlimK
        endif
! /*** find temperature ***/
        it = int((logTp-t_crtb(0))/dt)
        if(it.lt.0) then
          it = 0
        else if(it.ge.nt_crtb-1) then
          it = nt_crtb-2
        endif
        wt1 = (t_crtb(it+1)-logTp)/dt
        wt2 = (logTp-t_crtb(it))/dt

! /*** get cooling rates ***/
        cradp &
! iz-1
!   id
!     im
!       it
         =wt1*wm1*wd1*wz1*(cr_crtb(it,im,id,iz-1)) &
!       it+1
         +wt2*wm1*wd1*wz1*(cr_crtb(it+1,im,id,iz-1)) &
!     im+1
!       it
         +wt1*wm2*wd1*wz1*(cr_crtb(it,im+1,id,iz-1)) &
!       it+1
         +wt2*wm2*wd1*wz1*(cr_crtb(it+1,im+1,id,iz-1)) &
!   id+1
!     im
!       it
         +wt1*wm1*wd2*wz1*(cr_crtb(it,im,id+1,iz-1)) &
!       it+1
         +wt2*wm1*wd2*wz1*(cr_crtb(it+1,im,id+1,iz-1)) &
!     im+1
!       it
         +wt1*wm2*wd2*wz1*(cr_crtb(it,im+1,id+1,iz-1)) &
!       it+1
         +wt2*wm2*wd2*wz1*(cr_crtb(it+1,im+1,id+1,iz-1)) &
! iz
!   id
!     im
!       it
         +wt1*wm1*wd1*wz2*(cr_crtb(it,im,id,iz)) &
!       it+1
         +wt2*wm1*wd1*wz2*(cr_crtb(it+1,im,id,iz)) &
!     im+1
!       it
         +wt1*wm2*wd1*wz2*(cr_crtb(it,im+1,id,iz)) &
!       it+1
         +wt2*wm2*wd1*wz2*(cr_crtb(it+1,im+1,id,iz)) &
!   id+1
!     im
!       it
         +wt1*wm1*wd2*wz2*(cr_crtb(it,im,id+1,iz)) &
!       it+1
         +wt2*wm1*wd2*wz2*(cr_crtb(it+1,im,id+1,iz)) &
!     im+1
!       it
         +wt1*wm2*wd2*wz2*(cr_crtb(it,im+1,id+1,iz)) &
!       it+1
         +wt2*wm2*wd2*wz2*(cr_crtb(it+1,im+1,id+1,iz))

! /*** get heating rates ***/
        hradp &
! iz-1
!   id
!     im
!       it
         =wt1*wm1*wd1*wz1*(hr_crtb(it,im,id,iz-1)) &
!       it+1
         +wt2*wm1*wd1*wz1*(hr_crtb(it+1,im,id,iz-1)) &
!     im+1
!       it
         +wt1*wm2*wd1*wz1*(hr_crtb(it,im+1,id,iz-1)) &
!       it+1
         +wt2*wm2*wd1*wz1*(hr_crtb(it+1,im+1,id,iz-1)) &
!   id+1
!     im
!       it
         +wt1*wm1*wd2*wz1*(hr_crtb(it,im,id+1,iz-1)) &
!       it+1
         +wt2*wm1*wd2*wz1*(hr_crtb(it+1,im,id+1,iz-1)) &
!     im+1
!       it
         +wt1*wm2*wd2*wz1*(hr_crtb(it,im+1,id+1,iz-1)) &
!       it+1
         +wt2*wm2*wd2*wz1*(hr_crtb(it+1,im+1,id+1,iz-1)) &
! iz
!   id
!     im
!       it
         +wt1*wm1*wd1*wz2*(hr_crtb(it,im,id,iz)) &
!       it+1
         +wt2*wm1*wd1*wz2*(hr_crtb(it+1,im,id,iz)) &
!     im+1
!       it
         +wt1*wm2*wd1*wz2*(hr_crtb(it,im+1,id,iz)) &
!       it+1
         +wt2*wm2*wd1*wz2*(hr_crtb(it+1,im+1,id,iz)) &
!   id+1
!     im
!       it
         +wt1*wm1*wd2*wz2*(hr_crtb(it,im,id+1,iz)) &
!       it+1
         +wt2*wm1*wd2*wz2*(hr_crtb(it+1,im,id+1,iz)) &
!     im+1
!       it
         +wt1*wm2*wd2*wz2*(hr_crtb(it,im+1,id+1,iz)) &
!       it+1
         +wt2*wm2*wd2*wz2*(hr_crtb(it+1,im+1,id+1,iz))

        if(SI_flagrad.le.-2) then
          radp=10.0d0**(cradp-dlog10(ERU))-10.0d0**(hradp-dlog10(ERU))
          radp=radp*(10.0d0**(2.0d0*lognhp))
        else
          radp=10.0d0**(cradp+2.0d0*lognhp-dlog10(ERU))
        endif

!        ram(pn)=radp/sca_data(pn,2)
        sca_data(pn,27)=sca_data(pn,3)/(radp/sca_data(pn,2))*4.71d8
        ! tcool = u_p/ram
        endif
      end do

    end subroutine get_cool

end module gcdp_cool
