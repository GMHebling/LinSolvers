CC=gcc
CFLAGS=-Wall -g
LDFLAGS= -lm -lstdc++ 
OBJFILES = main.o funcoesAux.o metodosDiretos.o metodosIndiretos.o
TARGET = SolverLinear

all: $(TARGET)

$(TARGET): $(OBJFILES)
		$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~ 