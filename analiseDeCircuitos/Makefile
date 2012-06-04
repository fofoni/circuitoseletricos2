# Universidade Federal do Rio de Janeiro
# Escola Polit�cnica
# Circuitos El�tricos II / 1.2012
# Prof. Ant�nio Carlos Moreir�o de Queiroz
# Alunos: Marcelle de Souza Campos
#         Pedro Angelo Medeiros Fonini
# Programa de an�lise de circuitos no tempo para estudar os m�todos de Gear

CC = g++
CFLAGS = -ansi -pedantic -Wall -Wextra
CPP_FILES = myFunctions.cpp trabalhoCircuitosEletricos.cpp

TRABALHOCIRCUITOS = trabalhoCircuitosEletricos
EXECS = $(TRABALHOCIRCUITOS)


### Default target

$(TRABALHOCIRCUITOS) : $(CPP_FILES)
	$(CC) $(CFLAGS) -o $(TRABALHOCIRCUITOS) $(CPP_FILES)


### Clean target

.PHONY : clean

clean :
	rm -f *.o *.core ?*~ *.exe* $(EXECS)


### Help Target

.PHONY : help

help:
	@echo "The following are the valid targets for this Makefile:"
	@echo "    $(TRABALHOCIRCUITOS) (the default if no target is provided)"
	@echo "    clean"
	@echo "    help"