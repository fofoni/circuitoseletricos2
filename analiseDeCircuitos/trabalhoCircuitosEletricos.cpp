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

int main (int argc, char *argv[]) {

	ifstream myFile;
 	string line;
    elementsList list;

    myFile.open (argv[0]);

    cout << "ola!!\n" <<endl;

 	if (!(myFile.is_open())) {
 		cout << "Erro (" << errno << "): " << strerror (errno) << "." << endl;
 		exit (FILE_IS_NOT_OPEN);
 	}

    if (!myFile.good())
        cout << "Esse arquivo nao possui netlist" << endl;
 	while (myFile.good()) {
 		getline(myFile, line);
 		switch (line[0]){

                    case 'R':
                    case 'I':
                    case 'L':
 					case 'V':
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
                    default:
                        cout << "Esse elemento " << line[0] << " nao esta implementado" << endl;
                        break;
 		}
 	}
 	myFile.close();
}







