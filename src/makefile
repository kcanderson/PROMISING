
CC = g++
CFLAGS = -g -Wall -std=c++11 -O3

all: main.o completegraphscorer.o coreroutines.o pvaluescorer.o
	$(CC) $(CFLAGS) -o promising main.o completegraphscorer.o coreroutines.o pvaluescorer.o

main.o: main.cpp coreroutines.h CompleteGraphScorer.h PValueModuleScorer.h
	$(CC) $(CFLAGS) -c -I../include main.cpp -o main.o

completegraphscorer.o: ../include/IModuleScorer.h CompleteGraphScorer.h CompleteGraphScorer.cpp
	$(CC) $(CFLAGS) -c CompleteGraphScorer.cpp -o completegraphscorer.o

pvaluescorer.o: $(OBJDIR) PValueModuleScorer.cpp ../include/IModuleScorer.h PValueModuleScorer.h
	$(CC) $(CFLAGS) -c PValueModuleScorer.cpp -o pvaluescorer.o	

coreroutines.o: $(OBJDIR) coreroutines.cpp coreroutines.h 
	$(CC) $(CFLAGS) -c coreroutines.cpp -o coreroutines.o

clean:
	rm -f *.o
