# UPDATE PATHS FOR YO SELF
PGPLOT_DIR=/home/djw/pgplot
PGPLOT_LIB=-L/home/djw/pgplot -lpgplot -lpng -lz
LIBS=-L/usr/X11R6/lib -lX11 `$(PGPLOT_DIR)/cpg/libgcc_path.sh` -lgcc -lm -lc

FLAGS = -ffixed-line-length-none -O3 -mcmodel=medium
# FLAGS = -ffixed-line-length-none -g -O3 -mcmodel=medium -fbacktrace -fcheck=bounds

#F77 = ifort
F77 = gfortran

.SUFFIXES: .F90

.F90.o:
	$(F77) $(FLAGS) -c $<  $(PGPLOT_LIB)  $(LIBS)

OBJ = flatkernel.o plot_prams.o data_vals.o galtools.o tides.o gcdp_cool.o init_rgb_table.o read_prams.o read_data.o plot_data.o dotty_plotty.o do_montage.o do_anim.o do_strips.o

dotty_plotty: $(OBJ)
	$(F77) $(FLAGS) -o ../dotty_plotty $(OBJ) $(PGPLOT_LIB)  $(LIBS) 

clean:
	/bin/rm  -rf *.o *.mod ../dotty_plotty

