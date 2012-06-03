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
				 	 case (5):
				 		 ( (*this)[elementData])->parameter = elementData;
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

int elementsList::numberOfNodes() {

	int node = 0;

	map <string, element *> :: iterator auxiliar;

	for (auxiliar = this->begin(); auxiliar != this->end(); auxiliar++) {

		if ( ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode) > node )
			node = ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode);
	}

	return (node);
}

/*
string elementsList::locateCurrent (int node1, int node2){

	string name = NULL;

	map <string, element *> :: iterator auxiliar;

	for (auxiliar = this->begin(); auxiliar != this->end(); auxiliar++) {

			if ( (((auxiliar->second)-> originNodeOrPositiveOutputNode) == node1 ) &&
			     (((auxiliar->second)-> destinationNodeOrNegativeOutputNode) == node2 ))
			     name = auxiliar->first;

	}
	return (name);
}
*/

void elementsList::buildModifiedNodalMatrix (int *matrixOrder, modifiedMatrix matrix1, modifiedMatrix matrix3){

	tensionAndCurrentName listToPrint;
	int index = numberOfNodes();

	map <string, element *> :: iterator auxiliar;

	for (auxiliar = this->begin(), index++; auxiliar != this->end(); auxiliar++, index++) {

		switch (auxiliar->first[0]) {

			case 'R': /* Resistor */
					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [(auxiliar->second)->originNodeOrPositiveOutputNode] += ( 1/((auxiliar->second)->value) );

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
							[(auxiliar->second)->destinationNodeOrNegativeOutputNode] += ( 1/((auxiliar->second)->value) );

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
						    [(auxiliar->second)->destinationNodeOrNegativeOutputNode] += ( -1/((auxiliar->second)->value) );

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [(auxiliar->second)->originNodeOrPositiveOutputNode] += ( -1/((auxiliar->second)->value) );
				break;
			case 'E': /* Fonte de tensao controlada a tensao */

					matrix1 [index]
					          [(auxiliar->second)->originNodeOrPositiveOutputNode] += -1;

					matrix1 [index]
					         [(auxiliar->second)->destinationNodeOrNegativeOutputNode] += 1;

					matrix1 [index]
					[(auxiliar->second)->controledOriginNodeOrPositiveInputNode] += ((auxiliar->second)->value);

					matrix1 [index]
					[(auxiliar->second)->controledDestinationNodeOrNegativeInputNode] += -((auxiliar->second)->value);

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					       [index] += 1;

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
							[index] += -1;

					listToPrint[index] = 'j' + (auxiliar->first);
				break;
			case 'F': /* Fonte de corrente controlada a corrente */
					matrix1 [index]
					        [(auxiliar->second)->controledOriginNodeOrPositiveInputNode] += -1;

					matrix1 [index]
					        [(auxiliar->second)->controledDestinationNodeOrNegativeInputNode] += +1;

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [index] += 1;

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [index] += -1;

					matrix1 [(auxiliar->second)->controledOriginNodeOrPositiveInputNode]
					        [index] += ((auxiliar->second)->value);

					matrix1 [(auxiliar->second)->controledDestinationNodeOrNegativeInputNode]
					        [index] += -((auxiliar->second)->value);

					listToPrint[index] = 'j' + (auxiliar->first);
									/*	(locateCurrent (((auxiliar->second)->controledOriginNodeOrPositiveInputNode),
													   ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode))); */
				break;
			case 'G': /* Fonte de corrente controlada a tensao */
					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [(auxiliar->second)->controledOriginNodeOrPositiveInputNode] += ((auxiliar->second)->value);

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [(auxiliar->second)->controledDestinationNodeOrNegativeInputNode] += -((auxiliar->second)->value);

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [(auxiliar->second)->controledOriginNodeOrPositiveInputNode] += -((auxiliar->second)->value);

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [(auxiliar->second)->controledDestinationNodeOrNegativeInputNode] += ((auxiliar->second)->value);
				break;
			case 'H': /* Fonte de tensao controlada a corrente */
					matrix1 [index]
					        [(auxiliar->second)->originNodeOrPositiveOutputNode] += -1;

					matrix1 [index]
							[(auxiliar->second)->destinationNodeOrNegativeOutputNode] += 1;

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [index] += 1;

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [index] += -1;

					listToPrint[index] = 'j' + (auxiliar->first);

					index++;

					matrix1 [index]
					        [(auxiliar->second)->controledOriginNodeOrPositiveInputNode] += -1;

					matrix1 [index]
					        [(auxiliar->second)->controledDestinationNodeOrNegativeInputNode] += +1;

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [index] += 1;

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [index] += -1;

					matrix1 [index--]
					        [index] += ((auxiliar->second)->value);

					listToPrint[index] = 'j' + (auxiliar->first);
									/*	(locateCurrent (((auxiliar->second)->controledOriginNodeOrPositiveInputNode),
					    			    ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode))); */

				break;
			case 'I': /*Fonte de corrente*/
					matrix3 [(auxiliar->second)->originNodeOrPositiveOutputNode][0] += -((auxiliar->second)->value);

					matrix3 [(auxiliar->second)->destinationNodeOrNegativeOutputNode][0] += ((auxiliar->second)->value);
				break;
			case 'V': /*Fonte de tensao */
					matrix1 [index]
					        [(auxiliar->second)->originNodeOrPositiveOutputNode] += -1;

					matrix1 [index]
					        [(auxiliar->second)->destinationNodeOrNegativeOutputNode] += +1;

					matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
					        [index] += 1;

					matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
					        [index] += -1;

					listToPrint[index] = 'j' + (auxiliar->first);
				break;
			case 'L': /*indutor*/
				break;
		}
	}

	*matrixOrder = index;

}



/**************************************************************************************************/


/* Function responsable for solving the system A x = B */
void modifiedMatrix::solveMatrixSystem (int * matrixOrder, modifiedMatrix matrix1, modifiedMatrix matrix3) {

	int order = *matrixOrder;
	int i, j;

	float A[order][order];
	float x[order];
	float B[order];


	/* Building Matrix A */
	for (i=1; i == order; i++){
	    	for (j=1; j == order; j++){
	    		A[i][j] = 0;
	    		A[i][j] += matrix1 [i][j];
	    	}
	}

	/* Building Matrix B */
	for (i=1; i== order; i++){
		B[i] = 0;
		B[i] += matrix3[i][0];
	}


	/* Solving system and finding x */


}








