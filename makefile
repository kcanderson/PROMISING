BIN = ./bin

promising: all $(BIN)
	mv src/promising $(BIN)/promising

all:
	cd src; make

$(BIN):
	mkdir -p $(BIN)

clean:
	cd src; make clean
	rm -f $(BIN)/promising
