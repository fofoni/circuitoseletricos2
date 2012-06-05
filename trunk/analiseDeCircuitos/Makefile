# Universidade Federal do Rio de Janeiro
# Escola Politécnica
# Circuitos Elétricos II / 1.2012
# Prof. Antônio Carlos Moreirão de Queiroz
# Alunos: Marcelle de Souza Campos
#         Pedro Angelo Medeiros Fonini
# Programa de análise de circuitos no tempo para estudar os métodos de Gear

CC = g++
CFLAGS = -ansi -pedantic -Wall -Wextra

CPP_FILES = myFunctions.cpp trabalhoCircuitosEletricos.cpp
HEADERS = circuitAnalysis.h

TRABALHOCIRCUITOS = trabalhoCircuitosEletricos
EXECS = $(TRABALHOCIRCUITOS)


### Default target

$(TRABALHOCIRCUITOS) : $(CPP_FILES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TRABALHOCIRCUITOS) $(CPP_FILES)


### Clean target

.PHONY : clean

clean :
	rm -f *.o *.core ?*~ *.exe* $(EXECS)


### Help Target

.PHONY : help

help :
	@echo "The following are the valid targets for this Makefile:"
	@echo "    $(TRABALHOCIRCUITOS) (the default if no target is provided)"
	@echo "    clean"
	@echo "    help"
