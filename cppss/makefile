CC=g++
CFLAGS= -c -Wall
LIB= `root-config --libs | sed -e 's/-lPhysics//; s/-lPostscript//'`
CFL= `root-config --cflags`

SS: sistema.o vec.o corpo.o ilMain.o
	$(CC) $(CFL) sistema.o vec.o corpo.o ilMain.o -o SS $(LIB)

graphs: sistema.o vec.o corpo.o graphs.o
	$(CC) $(CFL) sistema.o vec.o corpo.o graphs.o -o graphs $(LIB)	

cc: vec.o corpo.o corpo_test.o
	$(CC) $(CFL) vec.o corpo.o corpo_test.o -o cc $(LIB)

vettore: vec.o vec_test.o
	$(CC) vec.o vec_test.o -o vettore
	
sonda: sonda.o vec.o corpo.o sonda_test.o
	$(CC) sonda.o vec.o corpo.o sonda_test.o -o sonda
	

ilMain.o: ilMain.cpp
	$(CC) $(CFLAGS) $(CFL) ilMain.cpp $(LIB)
graphs.o: graphs.cpp
	$(CC) $(CFLAGS) $(CFL) graphs.cpp $(LIB)
sistema.o: sistema.cpp sistema.h
	$(CC) $(CFLAGS) $(CFL) sistema.cpp $(LIB)
corpo.o: corpo.cpp corpo.h
	$(CC) $(CFLAGS) $(CFL) corpo.cpp $(LIB)
vec.o: vec.cpp vec.h
	$(CC) $(CFLAGS) vec.cpp
	

corpo_test.o: corpo_test.cpp
	$(CC) $(CFLAGS) $(CFL) corpo_test.cpp $(LIB)

vec_test.o: vec_test.cpp
	$(CC) $(CFLAGS) vec_test.cpp

sonda_test.o: sonda_test.cpp
	$(CC) $(CFLAGS) sonda_test.cpp
sonda.o: sonda.cpp
	$(CC) $(CFLAGS) sonda.cpp

aa: aa.cpp
	$(CC) $(CFL) aa.cpp -o aa $(LIB)

clean:
	rm -rf *.o SS cc vettore sonda
