CC         = g++
#CFLAGS     = -pg -fopenmp
CFLAGS     = -O3 -fopenmp -g  -Wall -I/home/hchao/install/include
LIBS      = -fopenmp  -lm  -O3 -lfftw3_omp -lfftw3 -lpthread -L/home/hchao/install/lib


#############################################################################
# nothing should be changed below here

SRCS = calc_press.cpp find_blkmu.cpp calc_diff.cpp poly_density.cpp dfnloc.cpp integ_utils.cpp io_utils.cpp array_utils.cpp main.cpp read_matdisp.cpp read_input.cpp initialize.cpp dfree.cpp\
       
       
			 


OBJS = ${SRCS:.cpp=.o}

.cpp.o:
	${CC} ${CFLAGS} ${DFLAGS} -c  $<

bi_dft:  ${OBJS}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS} $(LIBS)

clean:
	rm -f *.o
	rm -f bi_dft
	rm -f *~

