CPP = g++
RM = rm
CPPFLAGS = -framework Accelerate
PROGRAM = test
PROGOBJS = test.o functions.o writeout.o Sensitivity.o LeastSquaresSolver.o lapackwrapper.o

all: $(PROGRAM)

test: $(PROGOBJS)
	$(CPP) $(PROGOBJS) $(CPPFLAGS) -o test

clean:
	$(RM) $(PROGRAM) $(PROGOBJS)

test.o: test.hpp functions.hpp writeout.hpp Sensitivity.hpp

functions.o: functions.hpp

writeout.o: writeout.hpp

Sensitivity.o: Sensitivity.hpp writeout.hpp

LeastSquaresSolver.o: LeastSquaresSolver.hpp Sensitivity.hpp writeout.hpp

lapackwrapper.o: lapackwrapper.hpp