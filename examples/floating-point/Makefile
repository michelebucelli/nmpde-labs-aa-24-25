CXXFLAGS=-std=c++17
EXE=floating_point

.PHONY: all clean distclean

all: ${EXE}

${EXE}: ${EXE}.o
	${CXX} ${CXX_FLAGS} ${LDFLAGS} $< -o $@

clean:
	rm *.o

distclean: clean
	rm ${EXE}