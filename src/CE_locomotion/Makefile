PYTHON_CONFIG ?= python3-config


#LIBS := $(shell $(PYTHON_CONFIG) --embed --libs)
#LDFLAGS := $(shell $(PYTHON_CONFIG) --ldflags)
#REMOVE=-arch arm64 -arch x86_64
#REPLACE=
#CXXFLAGS0 := $(shell $(PYTHON_CONFIG) --embed --cflags)
#CXXFLAGS := $(subst $(REMOVE),$(REPLACE),$(CXXFLAGS0))

LIBS := $(shell $(PYTHON_CONFIG) --embed --libs)
LDFLAGS := $(shell $(PYTHON_CONFIG) --embed --ldflags)
CXXFLAGS := $(shell $(PYTHON_CONFIG) --includes)
LDFLAGS += "-L/opt/homebrew/lib"
#LDFLAGS += "-L$(brew --prefix nlohmann-json)/lib"
CXXFLAGS += "-I/opt/homebrew/include"
#CXXFLAGS += "-I$(brew --prefix nlohmann-json)/include"

main: info main.o jsonUtils.o argUtils.o Worm.o WormBody.o NervousSystem.o StretchReceptor.o Muscles.o TSearch.o random.o c302NervousSystem.o c302ForW2D.o owSignalSimulatorForWorm2D.o owSignalSimulator.o
	g++ $(CXXFLAGS) $(LDFLAGS) -pthread -o main main.o jsonUtils.o argUtils.o  Worm.o WormBody.o NervousSystem.o c302NervousSystem.o c302ForW2D.o owSignalSimulatorForWorm2D.o owSignalSimulator.o StretchReceptor.o Muscles.o TSearch.o random.o $(LIBS)


info:
	@echo ' '
	@echo '  - Compiling with '
	@echo '      CXXFLAGS: ' $(CXXFLAGS) 
	@echo '      LIBS: ' $(LIBS) 
	@echo '      LDFLAGS: ' $(LDFLAGS) 
	@echo ' '



random.o: random.cpp random.h VectorMatrix.h
	g++ -c -O3 -flto random.cpp
TSearch.o: TSearch.cpp TSearch.h
	g++ -c -O3 -flto TSearch.cpp

jsonUtils.o: jsonUtils.cpp jsonUtils.h
	g++ -c -O3 -std=c++11   -flto $(CXXFLAGS) $(LDFLAGS) jsonUtils.cpp	

argUtils.o: argUtils.cpp argUtils.h
	g++ -c -O3 -std=c++11 -flto argUtils.cpp
Worm.o: Worm.cpp Worm.h
	g++ -c -O3 -flto $(CXXFLAGS) $(LDFLAGS) Worm.cpp
WormBody.o: WormBody.cpp WormBody.h
	g++ -c -O3 -flto WormBody.cpp
NervousSystem.o: NervousSystem.cpp NervousSystem.h VectorMatrix.h random.h ../neuromlModel/NervousSystemBase.h
	g++ -c -O3 -flto NervousSystem.cpp
owSignalSimulatorForWorm2D.o: ../neuromlModel/owSignalSimulatorForWorm2D.cpp ../neuromlModel/owSignalSimulatorForWorm2D.h #../neuromlModel/owSignalSimulator.h 
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) ../neuromlModel/owSignalSimulatorForWorm2D.cpp
owSignalSimulator.o: ../neuromlModel/owSignalSimulator.cpp ../neuromlModel/owSignalSimulator.h #../neuromlModel/owINeuronSimulator.h
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) ../neuromlModel/owSignalSimulator.cpp
c302NervousSystem.o: ../neuromlModel/c302NervousSystem.cpp random.h #../neuromlModel/NervousSystemBase.h
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) ../neuromlModel/c302NervousSystem.cpp
c302ForW2D.o: ../neuromlModel/c302ForW2D.cpp ../neuromlModel/owSignalSimulatorForWorm2D.h #random.h #../neuromlModel/NervousSystemBase.h
	$(CC) -c -O3 $(CXXFLAGS) $(LDFLAGS) $(EXTRA_FLAGS) ../neuromlModel/c302ForW2D.cpp
StretchReceptor.o: StretchReceptor.cpp StretchReceptor.h
	g++ -c -O3 -flto StretchReceptor.cpp
Muscles.o: Muscles.cpp Muscles.h VectorMatrix.h random.h
	g++ -c -O3 -flto Muscles.cpp
main.o: main.cpp Worm.h WormBody.h StretchReceptor.h Muscles.h TSearch.h
	g++ -c -O3 -std=c++11 -flto $(CXXFLAGS) $(LDFLAGS) main.cpp
tests.o: tests.cpp NervousSystem.h 
	g++ -c -O3 -flto tests.cpp
tests: info tests.o NervousSystem.o random.o
	g++ -pthread -o tests tests.o NervousSystem.o random.o

tests2.o: tests2.cpp NervousSystem.h random.h jsonUtils.h argUtils.h
	g++ -c -O3 -flto $(CXXFLAGS) $(LDFLAGS) tests2.cpp
tests2: tests2.o NervousSystem.o random.o jsonUtils.o argUtils.o
	g++ $(CXXFLAGS) $(LDFLAGS) -pthread -o tests2 tests2.o NervousSystem.o random.o argUtils.o jsonUtils.o $(LIBS)

clean:
	rm -f *.o ../neuromlModel/*.o main tests tests2
