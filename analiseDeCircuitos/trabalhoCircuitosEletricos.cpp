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
    int qty_of_words;

    map<string, element*> :: iterator element = list.begin();
    map <string, long double[8]> :: iterator capacitorInductor = reactiveElements.begin();

    /* // coisas pra testar as funcoes:
    {
        cppmatrix A; cppmatrix B; cppmatrix b;
        A.initialize(4,4);
        b.initialize(4,1);

        A[1][1] = -12.78376391436733; A[1][2] = -8.08457973326789; A[1][3] = -1.08775891434418; A[1][4] = 9.70269659982691;
        A[2][1] =  -3.44936862205591; A[2][2] = 14.62951658126516; A[2][3] = -2.78872381570951; A[2][4] = 9.02562571827599;
        A[3][1] =  13.82450984933543; A[3][2] =  3.52830401914780; A[3][3] = 13.83591578750285; A[3][4] = 4.37974936399216;
        A[4][1] =  -5.34168532856833; A[4][2] =  5.03890307218638; A[4][3] =  4.59360560778174; A[4][4] = 4.99031594666725;
        b[1][1] =  -9.63642294769349; b[2][1] = 11.09251263328027; b[3][1] =  7.71853527401042; b[4][1] = 7.03641375842988;
        A.printMyself();
        A.solveMatrixSystem(b);

        cout << endl << endl << endl << "    =========" << endl << endl << endl;

        B[2][3] = 5;
        B[3][4] = -8;
        B[2][1] = .5;
        B.fill_out_with_zeros(5,4);
        B.printMyself();

        return 0;
    }
    */

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
      parametros sao guardados dentro da classe elementsList */
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
            cout << "Esse elemento " << split_line[0] << " nao esta implementado."
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

   /* list.buildModifiedNodalMatrix(listToPrint, matrix1, matrix3, reactiveElements, passo, gear_order, UIC);

    matrix1.printMyself();
    cout << endl;
    matrix3.printMyself();
    cout << endl;
    */


    while (element != list.end() ){
		capacitorInductor->first = element-> first;
		capacitorInductor->second[8] = element->second->originNodeOrPositiveOutputNode;
		capacitorInductor->second[9] = element->second->destinationNodeOrNegativeOutputNode;
    	if (UIC && (element->first[0] == 'L' || element->first[0] == "C" ))
    		capacitorInductor->second[0] = element->second->initialConditions;
    	else
    		capacitorInductor->second[0] = 0;
    	element++;
    }

    for (long double i = 0 ; i < (tempo_final * (passos_internos + 1) )/ passo; i++ ){
    	list.buildModifiedNodalMatrix(listToPrint, matrix1, matrix3, reactiveElements, passo, gear_order, UIC);
    	matrix2 = matrix1.solveMatrixSystem(matrix3);
    	list.printResult (argv, listToPrint, matrix2);
    	if ( i < 9 ){
    		for (capacitorInductor = reactiveElements.begin();
    			 capacitorInductor != list.end();
    			 capacitorInductor ++){

    			 capacitorInductor->second[i] = capacitorInductor->second[i-1];
    			 capacitorInductor->second[i-1]= matrix2[capacitorInductor->second[8]] - matrix2[capacitorInductor->second[9]];
    		}
    	}
    }







    //matrix1.solveMatrixSystem (matrixOrder, matrix1, matrix2, matrix3);
    matrix2 = matrix1.solveMatrixSystem(matrix3);

    matrix2.printMyself();

    return 0;

    list.printResult (argv, listToPrint, matrix2);

#if !(defined(unix) || defined(__unix__) || defined(__unix))
    cin.get();
#endif

    return 0;
}
