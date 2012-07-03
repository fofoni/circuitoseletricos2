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
    cppmatrix           matrix1; /* A */
    cppmatrix           matrix2; /* x */
    cppmatrix           matrix3; /* B */
    capacitor_inductor  reactiveElements;

    long double passo;
    long double tempo_final;
    int gear_order;
    int passos_internos;
    int UIC = -1;

    map<int, string> split_line;

    map<string, element*> :: iterator element = list.begin();
    capacitor_inductor :: iterator capacitorInductor = reactiveElements.begin();

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

    /* Durante a leitura do arquivo, os elementos assim como seus respectivos
      parametros sao guardados dentro da variável "list" */
    while (myFile.good()) {
        int qty_of_words;
        getline(myFile, line);
        cout << "Lida a linha [" << line << "]" << endl;
        if (line.size() == 0) continue;
        split_line = split(line, qty_of_words);
        switch (line[0]) {
          case 'R': case 'L': case 'C':
          case 'E': case 'F': case 'G': case 'H': case 'O':
          case 'N': case '$':
          case 'V': case 'I':
            list.getElement (line); break;
          case '.':
            if (qty_of_words != 4 && qty_of_words != 5) {
                cerr << ".TRAN <passo> <tempo final> GEAR[<n>] <passos internos> [UIC]" << endl;
                exit(BAD_NETLIST);
            }
            if (split_line[0].compare(".TRAN") !=  0) {
                cerr << "Simulacao '" << split_line[0] << "' nao reconhecida." << endl;
                exit(BAD_NETLIST);
            }
            passo = strtold(split_line[1].c_str(), NULL);
            tempo_final = strtold(split_line[2].c_str(), NULL);
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
          case '*': case '#': break; // comentarios da netlist
          default:
            cerr << "Esse elemento " << split_line[0] << " nao esta implementado."
                 << endl;
            exit(BAD_NETLIST);
            break;
        }
    }

    myFile.close();

    cout << endl;

//     for (elementsList::iterator i = list.begin(); i != list.end(); i++) {
//         cout << i->first << ":";
//         (i->second)->printMyself();
//         cout << endl;
//     }

    if (UIC == -1) {
        cerr << "Nenhuma simulacao a ser feita." << endl;
        exit(BAD_NETLIST);
    }

    list.buildModifiedNodalMatrix(listToPrint, matrix1, matrix3,
                                  reactiveElements, passo, gear_order, UIC);

    matrix1.printMyself();
    cout << endl;
    matrix3.printMyself();
    cout << endl;

    matrix2 = matrix1.solveMatrixSystem(matrix3);

    matrix2.printMyself();

    cout << endl;

    string header = "           t";
    for (int i = 1; i <= list.numberOfNodes(); i++) {
        char new_str[13];
        sprintf(new_str, " %11d", i);
        header = header + new_str;
    }
    for (int i = list.numberOfNodes()+1; i <= matrix1.n; i++) {
        char new_str[13];
        sprintf(new_str, " %11s", listToPrint[i].c_str());
        header = header + new_str;
    }
    cout << header << endl;
    cout << "           0";
    for (int i = 1; i <= matrix1.n; i++) {
        char new_str[13];
        sprintf(new_str, " % 11.4Lg", matrix2[i][1]);
        cout << new_str;
    }
    cout << endl;

    /*while (element != list.end() ){
        if (UIC)
            capacitorInductor->second[0] = element->second->initialConditions;
        else
            capacitorInductor->second[0] = 0;
        element++;
    }*/

    for (long int i = 1;
         i <= tempo_final*(passos_internos + 1)/passo; i++) {

        char new_str[13];
        list.buildModifiedNodalMatrix(listToPrint, matrix1, matrix3,
                                      reactiveElements, passo, gear_order, UIC);

        matrix2 = matrix1.solveMatrixSystem(matrix3);

        sprintf(new_str, " % 11.4Lg", i*passo/(passos_internos+1));
        cout << new_str;
        for (int i = 1; i <= matrix1.n; i++) {
            sprintf(new_str, " % 11.4Lg", matrix2[i][1]);
            cout << new_str;
        }
        cout << endl;
        //list.printResult (argv, listToPrint, matrix2);

        /*for (capacitorInductor = reactiveElements.begin();
            capacitorInductor != list.end();
            capacitorInductor ++){

            capacitorInductor->second[i] = capacitorInductor->second[i-1];
            capacitorInductor->second[i-1]= matrix2[capacitorInductor->second[8]] - matrix2[capacitorInductor->second[9]];
        }*/

    }

    return 0;

    //list.printResult (argv, listToPrint, matrix2);

#if !(defined(unix) || defined(__unix__) || defined(__unix))
    cin.get();
#endif

    return 0;
}
