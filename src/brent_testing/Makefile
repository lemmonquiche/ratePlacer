CXX = g++
CXXFLAGS = -g

SRCS =  main.cpp minfunc.c Msg.cpp Probability.cpp RandomVariable.cpp
OBJS = ${SRCS:.cpp=.o}
HEADERS = jph.h minfunc.h Msg.hpp Probability.hpp RandomVariable.hpp

MAIN = myprog

all: ${MAIN}
	@echo   Simple compilter named myprog has been compiled

${MAIN}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${MAIN}

.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@

clean:
	${RM} ${PROGS} ${OBJS} *.o *~.
