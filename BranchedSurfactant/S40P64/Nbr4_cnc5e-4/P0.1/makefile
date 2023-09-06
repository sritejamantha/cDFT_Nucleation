CC         = g++
#CFLAGS     = -pg -fopenmp
CFLAGS     = -O3 -fopenmp -g -Wall 
LIBS      = -fopenmp  -lm  -O3 -lpthread 


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

