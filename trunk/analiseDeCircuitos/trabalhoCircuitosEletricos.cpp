/* Universidade Federal do Rio de Janeiro
 * Escola Politécnica
 * Circuitos Elétricos II / 1.2012
 * Prof. Antônio Carlos Moreirão de Queiroz
 * Alunos: Marcelle de Souza Campos
 *         Pedro Angelo Medeiros Fonini
 * Programa de análise de circuitos no tempo para estudar os métodos de Gear
 */

#include "circuitsAnalyses.h"

using namespace std;

int main (int argc, char *argv[]) {

    ifstream myFile;
    string line;
    elementsList list;
    tensionAndCurrent   listToPrint;
    modifiedMatrix      matrix1; /* A */
    modifiedMatrix      matrix2; /* x */
    modifiedMatrix      matrix3; /* B */
    int *matrixOrder;   /* Ponteiro que guarda a ordem que a matriz A possuira */


    myFile.open (argv[1]);

    if (!(myFile.is_open())) {
        cout << "Erro (" << errno << "): " << strerror (errno) << "." << endl;
        exit (FILE_IS_NOT_OPEN);
    }

    getline(myFile, line);
    if (!myFile.good())
            cout << "Esse arquivo nao possui netlist" << endl;


    /* Durante a leitura do arquivo, os elementos assim como seus respectivos parametros sao guardados dentro da classe elementsList */
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

    list.buildModifiedNodalMatrix(matrixOrder, listToPrint, matrix1, matrix3);

    matrix1.solveMatrixSystem (matrixOrder, matrix1, matrix2, matrix3);

    list.printResult (argv, listToPrint, matrix2);



    cin.get();

    return 0;
}
