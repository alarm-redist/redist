TDZDD_DIR=.
OPT = -O3 -DB_64 -I.

all: enumpart

enumpart: enumpart.cpp GraphPartitionSpec.hpp GraphPartitionSpecL31.hpp BigInteger.hpp RandomSample.hpp testGP.hpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o
	g++ enumpart.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o enumpart -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG

enumpart-debug: enumpart.cpp GraphPartitionSpec.hpp GraphPartitionSpecL31.hpp testGP.hpp
	g++ enumpart.cpp -o enumpart-debug -I$(TDZDD_DIR) -DDEBUG

sample: sample.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o
	g++ sample.cpp SAPPOROBDD/bddc.o SAPPOROBDD/BDD.o SAPPOROBDD/ZBDD.o -o sample -I$(TDZDD_DIR) -std=c++11 -O3 -DB_64 -DNDEBUG

SAPPOROBDD/bddc.o: SAPPOROBDD/bddc.c SAPPOROBDD/bddc.h
	gcc $(OPT) SAPPOROBDD/bddc.c -c -o SAPPOROBDD/bddc.o

SAPPOROBDD/BDD.o: SAPPOROBDD/BDD.cc SAPPOROBDD/BDD.h
	g++ $(OPT) SAPPOROBDD/BDD.cc -c -o SAPPOROBDD/BDD.o

SAPPOROBDD/ZBDD.o: SAPPOROBDD/ZBDD.cc SAPPOROBDD/ZBDD.h
	g++ $(OPT) SAPPOROBDD/ZBDD.cc -c -o SAPPOROBDD/ZBDD.o

clean:
	rm -f enumpart enumpart.exe enumpart-debug enumpart-debug.exe sample sample.exe *.o SAPPOROBDD/bddc.o SAPPOROBDD/ZBDD.o
