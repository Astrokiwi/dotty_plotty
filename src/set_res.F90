    subroutine set_res
        use plot_prams, only : res_x, res_y
        implicit none

        ! Presumably your resolution is less than 100 million x 100 million
        character(len=8) :: res_x_str, res_y_str

        write(res_x_str,"(I8.8)") res_x
        write(res_y_str,"(I8.8)") res_y

        ! Just do it both ways and damn the consequences
        ! mu aha ha ha aha ha ha ha haaa
        call system("setenv PGPLOT_PNG_HEIGHT"//res_y_str)
        call system("setenv PGPLOT_PNG_WIDTH"//res_x_str)

        call system("export PGPLOT_PNG_HEIGHT="//res_y_str)
        call system("export PGPLOT_PNG_WIDTH="//res_x_str)

    end subroutine set_res