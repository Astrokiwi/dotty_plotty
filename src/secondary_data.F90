module secondary_data
    
    real(kind=8),save,allocatable,dimension(:) :: sec_data ! secondary data, some function of data read in
    character*4,save,allocatable,dimension(:) :: sec_codes
    integer,save :: nsec
    
contains
    
    subroutine secondary_set_codes
        implicit none
        nsec = 2
        allocate(sec_codes(nsec))
        
        sec_codes = ["cs_p",""]
        
    end subroutine secondary_set_codes

end module secondary_data