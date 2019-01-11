#FLAGS = -free -O2 -openmp -Vaxlib -mcmodel=medium -i-dynamic -lgfortran -L /home/williams/usr/local/pgplot -lpgplot -lpng -lz -L/usr/X11R6/lib -lX11 
#FLAGS = -O3 -mcmodel=medium -L /home/dwilliam/pgplot/ -lcpgplot -lpgplot -lpng -lz -L/usr/X11R6/lib -lX11
#FLAGS = -O3 -mcmodel=medium -L/home/dawil5/pgplotsrc -I/home/dawil5/pgplotsrc -lcpgplot -lpgplot -lpng -lz
#FLAGS = -g -fpe0 -mcmodel=medium -L /home/dwilliam/pgplot/ -lcpgplot -lpgplot -lpng -lz -L/usr/X11R6/lib -lX11
#FLAGS = -g -mcmodel=medium -L /home/dwilliam/pgplot/ -lcpgplot -lpgplot -lpng -lz -L/usr/X11R6/lib -lX11
# -lXmPgplot

#FLAGS = -O3 -mcmodel=medium -fpp -Dnopng


SRC=/home/djw/pgplot
LIBS=-L/usr/X11R6/lib -lX11 `$(SRC)/cpg/libgcc_path.sh` -lgcc -lm -lc
PGPLOT_LIB=-L/home/djw/pgplot/src -lpgplot -lpng -lz

FLAGS = -ffixed-line-length-none -O3 -mcmodel=medium
#FLAGS = -ffixed-line-length-none -g -fbounds-check -fbacktrace -mcmodel=medium

#FLAGS = -ffixed-line-length-none -O3 -mcmodel=medium  -L/home/djw/pgplot/src -lcpgplot -lpgplot -lpng -lz

# -lXmPgplot 
# -lpng -lz
#F77 = ifort
F77 = gfortran

.SUFFIXES: .F .F90

.F90.o:
	$(F77) $(FLAGS) -c $<  $(PGPLOT_LIB) $(LIBS)
.f.o:
	$(F77) $(FLAGS) -c $<
.F.o:
	@ if test $(SYS) = ibm ; then \
	echo "$(CPP) -P $(CPPFLAGS) $< > tmp/$*.f";\
	$(CPP) -P $(CPPFLAGS) $< > tmp/$*.f;\
	echo "$(F77) $(FLAGS) -c tmp/$*.f";\
	$(F77) $(FLAGS) -c tmp/$*.f;\
        elif test $(SYS) = linux ; then \
        echo "$(CPP) -P $(CPPFLAGS) $< > tmp/$*.f";\
        $(CPP) -P $(CPPFLAGS) $< > tmp/$*.f;\
        echo "$(F77) $(FLAGS) -c tmp/$*.f";\
        $(F77) $(FLAGS) -c tmp/$*.f;\
	else \
	echo "$(F77) $(FLAGS) $(CPPFLAGS) -c $<";\
	$(F77) $(FLAGS) $(CPPFLAGS) -c $<;\
	fi

.c.o:
	$(CC) -c  $(CPPFLAGS) $(CFLAGS) $< 

#%.o: %.mod

OBJ = flatkernel.o plot_prams.o data_vals.o galtools.o tides.o gcdp_cool.o init_rgb_table.o read_prams.o read_data.o plot_data.o dotty_plotty.o  do_montage.o do_anim.o do_strips.o

dotty_plotty: $(OBJ)
	$(F77) $(FLAGS) -o ../dotty_plotty $(OBJ)  $(PGPLOT_LIB) $(LIBS)

clean:
	/bin/rm  -rf *.o *.mod ../dotty_plotty

