CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas
CC = g++


main: huw.cc
	${CC} ${CFLAGS} -o huw huw.cc ${LDFLAGS}

sd: sd.cc
	${CC} ${CFLAGS} -o sd sd.cc ${LDFLAGS}
	
gt: gt.cc
	${CC} ${CFLAGS} -o gt gt.cc ${LDFLAGS}

prec: precompute.cc 
	${CC} ${CFLAGS} -o precompute precompute.cc ${LDFLAGS}
