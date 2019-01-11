subroutine read_data(intime)
    use data_vals
    implicit none
    integer, intent(in) :: intime
    
    integer :: ip
    real(kind=8) :: rad
    
    
    ! Invent some data!
    
    np = 200
    
    nsca = 3
    allocate(sca_codes(nsca))
    allocate(sca_data(np,nsca))
    sca_codes(1) = "rho"
    sca_codes(2) = "m"
    sca_codes(3) = "u_p"
    
    nvec = 1
    allocate(vec_codes(nvec))
    allocate(vec_data(3,np,nvec))
    vec_codes(1) = "vel"
    
    allocate(r_p(3,np))
    
    call random_number(r_p)
    
    r_p = r_p*.04 + .48
    
    sca_data(:,1) = sqrt((r_p(1,:)-.5)**2+(r_p(2,:)-.5)**2)
    
    
end subroutine read_data
