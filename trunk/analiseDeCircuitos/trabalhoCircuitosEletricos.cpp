/* Universidade Federal do Rio de Janeiro
   Escola Polit�cnica
   Circuitos El�tricos II / 1.2012
   Prof. Ant�nio Carlos Moreir�o de Queiroz
   Alunos: Marcelle de Souza Campos
   		   Pedro Angelo Medeiros Fonini
   Programa de an�lise de circuitos no tempo para estudar os m�todos de Gear
*/

#include "circuitsAnalyses.h"

using namespace std;

int main (int argc, char *argv[]) {

	ifstream myFile;
 	string line;
    elementsList list;
    modifiedMatrix matrix1; /* A */
    modifiedMatrix matrix3; /* B */
    int *matrixOrder;


    myFile.open (argv[1]);

 	if (!(myFile.is_open())) {
 		cout << "Erro (" << errno << "): " << strerror (errno) << "." << endl;
 		exit (FILE_IS_NOT_OPEN);
 	}

 	getline(myFile, line);
 	if (!myFile.good())
 	        cout << "Esse arquivo nao possui netlist" << endl;
 	    else
 	    	cout << "file is good" <<endl;

 	while (myFile.good()) {

 		getline(myFile, line);
 		switch (line[0]){

                    case 'R':
                    case 'L':
                    case 'C':
 						list.getElement (line, 1);
 						break;
 					case 'E':
                    case 'F':
                    case 'G':
                    case 'H':
                    case 'O':
						list.getElement (line, 2);
 						break;
 					case 'N':
						list.getElement (line, 3);
 						break;
 					case '$':
						list.getElement (line, 4);
 						break;
 					case 'V':
 					case 'I':
 						list.getElement (line, 5);
 						break;
                    default:
                        cout << "Esse elemento " << line[0] << " nao esta implementado" << endl;
                        break;
 		}
 	}
 	myFile.close();

 	list.buildModifiedNodalMatrix(matrixOrder, matrix1, matrix3);

 	matrix1.solveMatrixSystem(matrixOrder, matrix1, matrix3);


 	cin.get();
 	return 0;
}






