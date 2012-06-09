/* Universidade Federal do Rio de Janeiro
 * Escola Politécnica
 * Circuitos Elétricos II / 1.2012
 * Prof. Antônio Carlos Moreirão de Queiroz
 * Alunos: Marcelle de Souza Campos
 *         Pedro Angelo Medeiros Fonini
 * Programa de análise de circuitos no tempo para estudar os métodos de Gear
 */

#include "circuitAnalysis.h"

using namespace std;

int main (int argc, char *argv[]) {

    ifstream myFile;
    string line;
    elementsList list;
    tensionAndCurrent   listToPrint;
    modifiedMatrix      matrix1; /* A */
    modifiedMatrix      matrix2; /* x */
    modifiedMatrix      matrix3; /* B */
    int matrixOrder;   /* Ponteiro que guarda a ordem da matriz A */

    float passo;
    float tempo_final;
    int gear_order;
    int passos_internos;
    int UIC = -1;

    map<int, string> split_line;
    int qty_of_words;

    {
        cppmatrix A(4,4); cppmatrix B(4,4); cppmatrix x(4,1);
        A[1][1] =  1; A[1][2] =  1; A[1][3] =  1; A[1][4] =  1;
        A[2][1] =  1; A[2][2] = -1; A[2][3] =  1; A[2][4] = -1;
        A[3][1] =  1; A[3][2] =  1; A[3][3] = -1; A[3][4] = -1;
        A[4][1] =  1; A[4][2] = -1; A[4][3] = -1; A[4][4] =  1;
        B = A*(-1);
        x[1][1] =  1; x[2][1] =  1; x[3][1] =  2; x[4][1] =  3;
        // testa B+A=2A, B-A=0, B*A=4Id, A*x=?
    }

    if (argc < 2) {
        cerr << "Usage:" << endl;
        cerr << "  " << argv[0] << " <netlist file>" << endl;
        return FILE_IS_NOT_OPEN;
    }

    myFile.open (argv[1]);

    if (!(myFile.is_open())) {
        cerr << "Erro (" << errno << "): " << /*strerror (errno) <<*/ "." << endl;
        exit (FILE_IS_NOT_OPEN);
    }

    // pega o comentario
    getline(myFile, line);

    if (!myFile.good())
            cout << "Esse arquivo nao possui netlist" << endl;

    /* Durante a leitura do arquivo, os elementos assim como seus respectivos parametros sao guardados dentro da classe elementsList */
    while (myFile.good()) {
        getline(myFile, line);
        cout << "Lida a linha [" << line << "]" << endl;
        if (line.size() == 0) continue;
        switch (line[0]) {
          case 'R': case 'L': case 'C':
          case 'E': case 'F': case 'G': case 'H': case 'O':
          case 'N': case '$':
          case 'V': case 'I':
            list.getElement (line); break;
          case '.':
            split_line = split(line, qty_of_words);
            if (qty_of_words != 4 && qty_of_words != 5) {
                cerr << ".TRAN <passo> <tempo final> GEAR[<n>] <passos internos> [UIC]" << endl;
                exit(BAD_NETLIST);
            }
            if (split_line[0].compare(".TRAN") !=  0) {
                cerr << "Simulacao '" << split_line[0] << "' nao reconhecida." << endl;
                exit(BAD_NETLIST);
            }
            passo = atof(split_line[1].c_str());
            tempo_final = atof(split_line[2].c_str());
            if (split_line[3].size() < 4 ||
                split_line[3].substr(0, 4).compare("GEAR") != 0) {
                cerr << "Metodo '" << split_line[3] << "' nao reconhecido." << endl;
                exit(BAD_NETLIST);
            }
            if (split_line[3].size() == 4) gear_order = 2;
            else gear_order = atoi(split_line[3].substr(4).c_str());
            if (gear_order < 1 || gear_order > 8) {
                cerr << "A ordem do metodo de GEAR deve estar entre 1 e 8" << endl;
                exit(BAD_NETLIST);
            }
            passos_internos = atoi(split_line[4].c_str());
            if (qty_of_words == 4) UIC = 0;
            else {
                if (split_line[5].compare("UIC") != 0) {
                    cerr << "Comando '" << split_line[5] << "' nao reconhecido." << endl;
                    exit(BAD_NETLIST);
                }
                UIC = 1;
            }
          case '*': case '#': break; // comentarios
          default:
            cout << "Esse elemento " << line[0] << " nao esta implementado"
                 << endl;
            break;
        }
    }

    myFile.close();

    cout << endl;

    for (elementsList::iterator i = list.begin(); i != list.end(); i++) {
        cout << i->first << ":";
        (i->second)->printMyself();
        cout << endl;
    }

    if (UIC == -1) {
        cerr << "Nenhuma simulacao a ser feita." << endl;
        exit(BAD_NETLIST);
    }

    list.buildModifiedNodalMatrix(matrixOrder, listToPrint, matrix1, matrix3);

    matrix1.printMyself();
    cout << endl;
    matrix3.printMyself();

    matrix1.solveMatrixSystem (matrixOrder, matrix1, matrix2, matrix3);

    list.printResult (argv, listToPrint, matrix2);

#if !(defined(unix) || defined(__unix__) || defined(__unix))
    cin.get();
#endif

    return 0;
}