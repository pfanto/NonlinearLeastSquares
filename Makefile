CPP = g++
RM = rm
PROGRAM = test
PROGOBJS = test.o functions.o writeout.o Sensitivity.o

all: $(PROGRAM)

test: $(PROGOBJS)
	$(CPP) $(PROGOBJS) -o test

clean:
	$(RM) $(PROGRAM) $(PROGOBJS)

test.o: test.hpp functions.hpp writeout.hpp Sensitivity.hpp

functions.o: functions.hpp

writeout.o: writeout.hpp

Sensitivity.o: Sensitivity.hpp