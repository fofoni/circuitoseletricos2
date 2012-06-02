/* Universidade Federal do Rio de Janeiro
   Escola Politécnica
   Circuitos Elétricos II / 1.2012
   Prof. Antônio Carlos Moreirão de Queiroz
   Alunos: Marcelle de Souza Campos
   		   Pedro Angelo Medeiros Fonini
   Programa de análise de circuitos no tempo para estudar os métodos de Gear
*/

#include "circuitsAnalyses.h"

using namespace std;

void elementsList::getElement (string line, int situation) {

	int auxiliar = 0;
	int aux = 0;
	int quantityOfSpaces = 0;
	string elementData;

	do{
		while (line [auxiliar] != ' ')
			auxiliar++;

		for (; aux ++; aux == auxiliar)
			elementData [aux] = line [aux];

		switch (quantityOfSpaces) {
			case 0:
				(*this)[elementData] = new element;
				break;
			case 1:
				( (*this)[elementData] )->originNodeOrPositiveOutputNode = atoi (elementData.c_str());
				break;
			case 2:
				( (*this)[elementData] )-> destinationNodeOrNegativeOutputNode = atoi (elementData.c_str());
				break;
			case 3:
				 switch (situation){
				 	 case (1):
				 		( (*this)[elementData] )->value = atof (elementData.c_str());
						break;
				 	 case (2):
				 		( (*this)[elementData] )->controledOriginNodeOrPositiveInputNode = atoi (elementData.c_str());
				 	 	 break;
				 	 case (3):
				 		( (*this)[elementData] )->pairsOfValues = elementData;
				 	 	 break;
				 	 case (4):
				 	 	 ( (*this)[elementData] )->nocrtlPositive = atoi (elementData.c_str());
				 	 break;
				 }
				 break;
			case 4:
				 switch (situation){
					case (1):
						( (*this)[elementData] )->inicialConditions = atof (elementData.c_str());
						break;
					case (2):
						( (*this)[elementData] )->controledDestinationNodeOrNegativeInputNode = atoi (elementData.c_str());
					 	 break;
					case (4):
					 	 ( (*this)[elementData] )->nocrtlNegative = atoi (elementData.c_str());
						 break;
				}
				break;
			case 5:
				switch (situation){
					case (2):
						( (*this)[elementData] )->value = atof (elementData.c_str());
						break;
					case (4):
						( (*this)[elementData] )->gon = atof (elementData.c_str());
						break;
				}
				break;
			case 6:
				( (*this)[elementData] )->goff = atof (elementData.c_str());
				break;
			case 7:
				( (*this)[elementData] )->vref = atof (elementData.c_str());
				break;
		}

		quantityOfSpaces ++;
		aux ++;

	} while (line [auxiliar] != '\n');
}







