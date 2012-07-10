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

inline
string rm_filename_extension(string filename) {
    return filename.substr(0,filename.find('.'));
}

int main (int argc, char *argv[]) {

    ifstream myFile;
    string line;
    elementsList list;
    tensionAndCurrent   listToPrint;

    // sistema linear: Ax = b
    cppmatrix           matrix1; /* A */
    cppmatrix           matrix2; /* x */
    cppmatrix           matrix3; /* b */
    capacitor_inductor  reactiveElements;

    ofstream answerFile;
    string input_filename;
    string output_filename;

    long double passo;
    long double tempo_final;
    int gear_order;
    int passos_internos;
    int UIC = -1;

    map<int, string> split_line;

    cppmatrix solucao_anterior;

    if (argc < 2) {
        cout << "Usage:" << endl;
        cout << "  " << argv[0] << " <netlist file>" << endl;
        return FILE_IS_NOT_OPEN;
    }

    input_filename = argv[1];
    myFile.open (input_filename.c_str());

    if (!(myFile.is_open())) {
        cout << "Erro (" << errno << "): " << strerror (errno) << "." << endl;
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
        if (line[line.length()-1] == '\r')
            line = line.substr(0,line.length()-1);
//         cout << "Lida a linha [" << line << "]" << endl;
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
                cout << ".TRAN <passo> <tempo final> GEAR[<n>] <passos internos> [UIC]" << endl;
                exit(BAD_NETLIST);
            }
            if (split_line[0].compare(".TRAN") !=  0) {
                cout << "Simulacao '" << split_line[0] << "' nao reconhecida." << endl;
                exit(BAD_NETLIST);
            }
            tempo_final = strtold(split_line[1].c_str(), NULL);
            passo = strtold(split_line[2].c_str(), NULL);
            if (split_line[3].size() < 4 ||
                split_line[3].substr(0, 4).compare("GEAR") != 0) {
                cout << "Metodo '" << split_line[3] << "' nao reconhecido." << endl;
                exit(BAD_NETLIST);
            }
            if (split_line[3].size() == 4) gear_order = 2;
            else gear_order = atoi(split_line[3].substr(4).c_str());
            if (gear_order < 1 || gear_order > 8) {
                cout << "A ordem do metodo de GEAR deve estar entre 1 e 8" << endl;
                exit(BAD_NETLIST);
            }
            passos_internos = atoi(split_line[4].c_str());
            if (qty_of_words == 4) UIC = 0;
            else {
                if (split_line[5].compare("UIC") != 0) {
                    cout << "Comando '" << split_line[5] << "' nao reconhecido." << endl;
                    exit(BAD_NETLIST);
                }
                UIC = 1;

            }
            break;
          case '*': case '#': break; // comentarios da netlist
          default:
            cout << "Esse elemento " << split_line[0]
                 << " nao esta implementado." << endl;
            exit(BAD_NETLIST);
            break;
        }
    }

    myFile.close();


//     cout << endl;

//     for (elementsList::iterator i = list.begin(); i != list.end(); i++) {
//         cout << i->first << ":";
//         (i->second)->printMyself();
//         cout << endl;
//     }

    if (UIC == -1) {
        cout << "Nenhuma simulacao a ser feita." << endl;
        exit(BAD_NETLIST);
    }

    for (map<string, element*>::iterator list_iter = list.begin();
         list_iter != list.end();
         list_iter++)
    {
        if ((list_iter->first[0] != 'V') or (list_iter->first[0] != 'I'))
            if (list_iter->second->parameter.compare("PULSE") == 0) {
                if (list_iter->second->t_rise <= 0)
                    list_iter->second->t_rise = passo/passos_internos;
                if (list_iter->second->t_fall <= 0)
                    list_iter->second->t_fall = passo/passos_internos;
            }
        if ((list_iter->first[0] != 'L') and (list_iter->first[0] != 'C'))
            continue;
        if (UIC==1)
            for (int k=0; k<8; k++)
                reactiveElements[list_iter->first][k] = list_iter->second->initialConditions;
        else
            for (int k=0; k<8; k++)
                reactiveElements[list_iter->first][k] = 0;
    }

#ifdef OUTPUT_MATLAB
    output_filename = rm_filename_extension(input_filename) + ".m";
#else
    output_filename = rm_filename_extension(input_filename) + ".tab";
#endif
    answerFile.open (output_filename.c_str());

    // pré-simulação
    {

        char new_str[25];
        long double erro;
        int qty_of_trials, qty_of_guesses;
        bool singular;

        matrix2.initialize(list.numberOfNodes(), 1);

        list.buildModifiedNodalMatrix(listToPrint, matrix1, matrix2, matrix3,
                                    reactiveElements,
                                    passo/passos_internos/1e9,
                                    1, 0);
        solucao_anterior = matrix1.solveMatrixSystem(matrix3, singular);

//         cout << "Esta iteracao: A = "; matrix1.printMyself();
//         cout << "               x = "; solucao_anterior.printMyself();
//         cout << "               b = "; matrix3.printMyself();
//         cout << endl;

        qty_of_trials = qty_of_guesses = 0;
        srand (time(NULL));

        NR_new_guess_t0:
        do {

            list.buildModifiedNodalMatrix(listToPrint, matrix1, solucao_anterior, matrix3,
                                          reactiveElements,
                                          passo/passos_internos/1e9, 1, 0);

            matrix2 = matrix1.solveMatrixSystem(matrix3, singular);

//             cout << "Esta iteracao: A = "; matrix1.printMyself();
//             cout << "               x = "; matrix2.printMyself();
//             cout << "               b = "; matrix3.printMyself();
//             cout << endl;

            if (singular) {
                qty_of_trials = 50;
                matrix2.initialize(matrix1.n, 1);
            }

            erro = sqrt(((matrix2-solucao_anterior).t() * (matrix2-solucao_anterior))[1][1]);

            if (++qty_of_trials >= 50) {
                long double norm;
                if (++qty_of_guesses > 50) {
                    cout << "O Método de Newton-Raphson não converge." << endl;
                    exit(TOO_MANY_NEWTON_RAPHSON);
                }
                cout << "O método de Newton-Raphson não convergiu em 10 iterações." << endl
                     << "Tentando de novo." << endl;
                qty_of_trials = 0;
                norm = sqrt((matrix2.t() * matrix2)[1][1]);
                if (norm<1) norm=1;
                for (int k=1; k<=matrix1.n; k++) {
                    long double extra = 2*((long double)(rand()))/((long double)(RAND_MAX))-1;
//                     cout << matrix2[k][1] << " + " << 100*norm << " * " << extra << " = " << solucao_anterior[k][1] + 100*norm*extra << endl;
                    if (solucao_anterior[k][1] != matrix2[k][1])
                        solucao_anterior[k][1] = matrix2[k][1] + 1000*norm*extra;
                }
                goto NR_new_guess_t0;
            }

            solucao_anterior = matrix2;

        } while (erro > 1e-9 or qty_of_trials==1);

#ifdef OUTPUT_MATLAB
        string header = "A=[%               t";
#else
        string header = "                   t";
#endif
        for (int i = 1; i <= list.numberOfNodes(); i++) {
            sprintf(new_str, " %19d", i);
            header = header + new_str;
        }
        for (int i = list.numberOfNodes()+1; i <= matrix1.n; i++) {
            sprintf(new_str, " %19s", listToPrint[i].c_str());
            header = header + new_str;
        }
        answerFile << header << endl;

        answerFile << "                   0";
        for (int i = 1; i <= matrix1.n; i++) {
            sprintf(new_str, " % 19.12Lg", matrix2[i][1]);
            answerFile << new_str;
        }
        answerFile << endl;

        for (capacitor_inductor::iterator reactive_iter = reactiveElements.begin();
             reactive_iter != reactiveElements.end();
             reactive_iter++) {
            if (reactive_iter->first[0] == 'L')
                reactive_iter->second[0] = matrix2[int(reactive_iter->second[8])][1];
            else {
                matrix2[0][1] = 0;
                reactive_iter->second[0] = matrix2[int(reactive_iter->second[8])][1]
                                         - matrix2[int(reactive_iter->second[9])][1];
            }
            reactive_iter->second[7] = reactive_iter->second[6] =
            reactive_iter->second[5] = reactive_iter->second[4] =
            reactive_iter->second[3] = reactive_iter->second[2] =
            reactive_iter->second[1] = reactive_iter->second[0];
        }

    }

    for (long int i = 1;
         i <= tempo_final*passos_internos/passo; i++)
    {

        char new_str[25];
        long double erro;
        int qty_of_trials, qty_of_guesses;
        bool singular;

//         cout << "==============>> I!:" << i << endl;

        qty_of_trials = qty_of_guesses = 0;
        NR_new_guess:
        do {

            list.buildModifiedNodalMatrix(listToPrint, matrix1, solucao_anterior, matrix3,
                                        reactiveElements, passo/passos_internos,
                                        gear_order, i*passo/passos_internos);
            matrix2 = matrix1.solveMatrixSystem(matrix3, singular);

            if (singular) {
                qty_of_trials = 50;
                matrix2.initialize(matrix1.n, 1);
            }

//             cout << "Esta iteracao: A = "; matrix1.printMyself();
//             cout << "               x = "; matrix2.printMyself();
//             cout << "               b = "; matrix3.printMyself();
//             cout << endl;

            erro = sqrt(((matrix2-solucao_anterior).t() * (matrix2-solucao_anterior))[1][1]);

            if (++qty_of_trials > 50) {
                long double norm;
                if (++qty_of_guesses > 50) {
                    cout << "O Método de Newton-Raphson não converge." << endl;
                    exit(TOO_MANY_NEWTON_RAPHSON);
                }
                cout << "O método de Newton-Raphson não convergiu em 10 iterações." << endl
                     << "Tentando de novo." << endl;
                qty_of_trials = 0;
                norm = sqrt((matrix2.t() * matrix2)[1][1]);
                if (norm<1) norm=1;
                for (int k=1; k<=matrix1.n; k++) {
                    long double extra = 2*((long double)(rand()))/((long double)(RAND_MAX))-1;
//                     cout << matrix2[k][1] << " + " << 100*norm << " * " << extra << " = " << solucao_anterior[k][1] + 100*norm*extra << endl;
                    if (solucao_anterior[k][1] != matrix2[k][1])
                        solucao_anterior[k][1] = matrix2[k][1] + 100*norm*extra;
                }
                goto NR_new_guess;
            }

            solucao_anterior = matrix2;

        } while (erro > 1e-9 or qty_of_trials==1);

        if (i % passos_internos == 0) {
            sprintf(new_str, " % 19.12Lg", i*passo/passos_internos);
            answerFile << new_str;
            for (int k = 1; k <= matrix1.n; k++) {
                sprintf(new_str, " % 19.12Lg", matrix2[k][1]);
                answerFile << new_str;
            }
            answerFile << endl;
        }

        for (capacitor_inductor::iterator reactive_iter = reactiveElements.begin();
             reactive_iter != reactiveElements.end();
             reactive_iter++) {
            reactive_iter->second[7] = reactive_iter->second[6];
            reactive_iter->second[6] = reactive_iter->second[5];
            reactive_iter->second[5] = reactive_iter->second[4];
            reactive_iter->second[4] = reactive_iter->second[3];
            reactive_iter->second[3] = reactive_iter->second[2];
            reactive_iter->second[2] = reactive_iter->second[1];
            reactive_iter->second[1] = reactive_iter->second[0];
            if (reactive_iter->first[0] == 'L')
                reactive_iter->second[0] = matrix2[int(reactive_iter->second[8])][1];
            else {
                matrix2[0][1] = 0;
                reactive_iter->second[0] = matrix2[int(reactive_iter->second[8])][1]
                                         - matrix2[int(reactive_iter->second[9])][1];
            }
        }

    }

#ifdef OUTPUT_MATLAB
    answerFile << "];";
    answerFile << "t=A(:,1); " << endl;
    for (int i = 1; i <= list.numberOfNodes(); i++) {
        char new_str[20];
        sprintf(new_str, "e%d=A(:,%d); ", i, i+1);
        answerFile << new_str;
    }
    answerFile << endl;
    for (int i = list.numberOfNodes()+1; i <= matrix1.n; i++) {
        char new_str[20];
        sprintf(new_str, "%s=A(:,%d); ", listToPrint[i].c_str(), i+1);
        answerFile << new_str;
    }
    answerFile << endl;
#endif

    answerFile.close();

    cout << "Simulação sucedida. Arquivo salvo em " << output_filename << endl;

    return 0;
}
