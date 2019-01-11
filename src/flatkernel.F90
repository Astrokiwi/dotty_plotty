module flatkernel
    
    integer, parameter :: nkern = 1000
    real(kind=8),dimension(nkern) :: kern_tab
    real(kind=8), save :: kern_norm

!    real(kind=8), parameter :: smoo_norm = 8./3.14159265359
    
    contains
    
    subroutine kernel_init
        implicit none
        
        integer :: ir,iz
        real(kind=8) :: r,z,x
        integer, parameter :: zsteps = 1000
        
        real(kind=8) :: norm
        
        kern_tab = 0.
!        sum = 0.
        
        do ir=1,nkern
            r = (ir-1.)/(nkern-1.)
            do iz=-zsteps,zsteps-1
                z = (iz+.5)/zsteps
                x = sqrt(r**2+z**2)
                kern_tab(ir) = kern_tab(ir) + kern(x)/zsteps
            end do
        end do

        ! normalise flattened kernel
        do ir=1,nkern
            r = (ir-1.)/(nkern-1.)
            norm = norm+kern_tab(ir)/(nkern-1.)*3.14*r
!            write(*,"(I5,3E15.5)") ir,r,kern_tab(ir),sum
        end do
        kern_tab(:) = kern_tab(:)/norm ! Normalise so it all sums to 1.
        
        ! normalise 3d kernel
        kern_norm = 0.
        do ir=1,nkern
            r = (ir-1.)/(nkern-1.)
            kern_norm = kern_norm+fkern(r)/(nkern-1.)*3.14*r**2
        end do
        
        
    end subroutine kernel_init
    
    function kern(x)
        implicit none
        real(kind=8) :: kern
        real(kind=8) :: x

        if ( x<0.5 ) then
            kern = 1-6*x**2+6*x**3
        else if ( x<=1 ) then
            kern = 2*(1-x)**3
        else
            kern = 0.
        endif
        return
    end function kern
    
    function fkern(x)
        implicit none
        
        real(kind=8) :: x

        integer :: ix
        real(kind=8) :: fw,fs
        real(kind=8) :: fkern
        
        ix = x*(nkern-1)+1
!        if ( ix<=0 ) then
!            print *,ix,x,nkern
!            stop
!        endif
        
        if ( ix>=nkern ) then
            fkern = 0.
            return
        else
        
            fw = x-(ix-1.)/(nkern-1.)
            fs = 1./(nkern-1)
        
            fkern = kern_tab(ix)+fw/fs*(kern_tab(ix+1)-kern_tab(ix))
            if ( fkern<0 ) then
                stop '<0'
            endif
            return
        endif
        
    end function fkern
    
end module flatkernel
