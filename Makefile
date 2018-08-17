CPP = g++
RM = rm
PROGRAM = test
PROGOBJS = test.o functions.o writeout.o

all: $(PROGRAM)

test: $(PROGOBJS)
	$(CPP) $(PROGOBJS) -o test

clean:
	$(RM) $(PROGRAM) $(PROGOBJS)

test.o: test.hpp functions.hpp writeout.hpp

functions.o: functions.hpp

writeout.o: writeout.hpp