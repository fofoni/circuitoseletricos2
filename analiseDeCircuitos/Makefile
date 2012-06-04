# Universidade Federal do Rio de Janeiro
#  Escola Politécnica
#  Circuitos Elétricos II / 1.2012
#  Prof. Antônio Carlos Moreirão de Queiroz
#  Alunos: Marcelle de Souza Campos
#  		   Pedro Angelo Medeiros Fonini
#  Programa de análise de circuitos no tempo para estudar os métodos de Gear

CC = gcc;
CL = gcc;

CFLAGS = -Wall
LFLAGS = -Wall

TRABALHOOBJS = myFunctions.o trabalhoCircuitosEletricos.o

TRABALHOCIRCUITOS = trabalhoCircuitosEletricos

EXECS = $(TRABALHOCIRCUITOS)

# Regra implicita
.c.o:
	$(CC) $(CFLAGS) -c $<
	
all: $(EXECS)

circuitosEletricosII: $(TRABALHOCIRCUITOS)

trabalhoCircuitosEletricos: $(TRABALHOOBJS)
	$(LD) $(LFLAGS) -o $@ $(TRABALHOOBJS)
	
clean:
	rm -f *.o *.core *.?~ $(EXECS)
	
#End Makefile 
 
