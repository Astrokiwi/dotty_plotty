module data_vals
!    real(kind=8),allocatable,dimension(:,:) :: r_p,v_p ! position, velocity
!    real(kind=8),allocatable,dimension(:) :: m_p,rho_p,u_p ! mass, density, internal energy
    real(kind=8),save,allocatable,dimension(:,:) :: r_p ! Everything must have position, right?

    real(kind=8),save,allocatable,dimension(:,:,:) :: vec_data
    character*3,save,allocatable,dimension(:) :: vec_codes
    integer,save :: nvec

    real(kind=8),save,allocatable,dimension(:,:) :: sca_data
    character*4,save,allocatable,dimension(:) :: sca_codes
    character*16,save,allocatable,dimension(:) :: sca_units
    integer,save :: nsca

    integer,save,allocatable,dimension(:) :: itype
    integer,save,allocatable,dimension(:) :: id_p
    integer,save,allocatable,dimension(:) :: ifeedtype
    
    integer, parameter ::   itstar = 0,&
                            itgas = 1,&
                            itdark = 2,&
                            itany = -1,&
                            itfeed = 3 ! Maybe have multiple booleans?
    
    integer,save :: np

    real(kind=8),save :: time

end module data_vals
