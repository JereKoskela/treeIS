CFLAGS=-Wall -Wextra -O3
LDFLAGS= -lgsl -lgslcblas -lconfig++
CC = g++


main: sd.cc 
	${CC} ${CFLAGS} -o sd sd.cc ${LDFLAGS}

gt: gt.cc 
	${CC} ${CFLAGS} -o gt gt.cc ${LDFLAGS}
