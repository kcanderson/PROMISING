
CC = g++
CFLAGS = -g -Wall -std=c++11 -O2
OBJDIR = bin

all: $(OBJDIR)/main.o $(OBJDIR) $(OBJDIR)/completegraphscorer.o $(OBJDIR)/coreroutines.o $(OBJDIR)/pvaluescorer.o
	$(CC) $(CFLAGS) -o $(OBJDIR)/promising bin/main.o $(OBJDIR)/completegraphscorer.o $(OBJDIR)/coreroutines.o $(OBJDIR)/pvaluescorer.o

$(OBJDIR)/main.o: $(OBJDIR) src/main.cpp src/coreroutines.h src/CompleteGraphScorer.h src/PValueModuleScorer.h
	$(CC) $(CFLAGS) -c -Iinclude src/main.cpp -o $(OBJDIR)/main.o

$(OBJDIR)/completegraphscorer.o: include/IModuleScorer.h src/CompleteGraphScorer.h src/CompleteGraphScorer.cpp
	$(CC) $(CFLAGS) -c src/CompleteGraphScorer.cpp -o $(OBJDIR)/completegraphscorer.o

$(OBJDIR)/pvaluescorer.o: $(OBJDIR) src/PValueModuleScorer.cpp include/IModuleScorer.h src/PValueModuleScorer.h
	$(CC) $(CFLAGS) -c src/PValueModuleScorer.cpp -o $(OBJDIR)/pvaluescorer.o	

$(OBJDIR)/coreroutines.o: $(OBJDIR) src/coreroutines.cpp src/coreroutines.h 
	$(CC) $(CFLAGS) -c src/coreroutines.cpp -o $(OBJDIR)/coreroutines.o

clean:
	rm $(OBJDIR)/*.o
	rm $(OBJDIR)/promising

$(OBJDIR):
	mkdir -p $(OBJDIR)
