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

/* funcoes para manipulacao de matrizes */

void cppmatrix::initialize (int rows, int columns) {
    n=rows; m=columns;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            (*this)[i][j] = 0;
}

cppmatrix cppmatrix::operator+ (cppmatrix parcela) {
    cppmatrix resultado;
    resultado.initialize(n,m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            resultado[i][j] = (*this)[i][j] + parcela[i][j];
    return resultado;
}

cppmatrix cppmatrix::operator* (cppmatrix fator) {
    cppmatrix resultado;
    resultado.initialize(n, fator.m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= fator.m; j++)
            for (int k = 1; k <= m; k++)
                resultado[i][j] += (*this)[i][k] * fator[k][j];
    return resultado;
}

cppmatrix cppmatrix::operator* (long double factor) {
    cppmatrix resultado;
    resultado.initialize(n,m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            resultado[i][j] = factor * (*this)[i][j];
    return resultado;
}

cppmatrix cppmatrix::t () {
    cppmatrix resultado;
    resultado.initialize(m,n);
    for (int i = 1; i <= m; i++)
        for (int j = 1; j <= n; j++)
            resultado[i][j] = (*this)[j][i];
    return resultado;
}

cppmatrix cppmatrix::submatrix (int i1, int i2, int j1, int j2) {
    cppmatrix resultado;
    resultado.initialize(i2-i1+1, j2-j1+1);
    for (int i = 1; i <= i2-i1+1; i++)
        for (int j = 1; j <= j2-j1+1; j++)
            resultado[i][j] = (*this)[i+i1-1][j+j1-1];
    return resultado;
}

void cppmatrix::make_id(int size) {
    initialize(size, size);
    for (int i=1; i<=size; i++)
        (*this)[i][i] = 1;
}

template <typename T>
int sgn(T val, int zero=1) {
    return zero*(T(0) == val) + (T(0) < val) - (val < T(0));
}

void cppmatrix::subassign(int i1, int i2, int j1, int j2, cppmatrix A) {
    for (int i=1; i<=i2-i1+1; i++)
        for (int j=1; j<=j2-j1+1; j++)
            (*this)[i+i1-1][j+j1-1] = A[i][j];
}

// Solves (for the vector x) the system Ax=b, with A = *this
// Equivalent to octave/matlab's A\b
cppmatrix cppmatrix::solveMatrixSystem(cppmatrix b) {
    if (n != m) {
        cerr << "ERRO: este algoritmo nao faz decomposicao QR"
                " de matrizes nao-quadradas" << endl;
        exit(NON_SQUARE_MATRIX_QR);
    }
    cppmatrix R = *this;
    cppmatrix x;
    x.initialize(n,1);
    for (int i=1; i<=n; i++) {
        // Find H = I - 2*v*v'
        cppmatrix u = R.submatrix(i,n, i,i);
        long double normx = sqrt((u.t() * u)[1][1]);
        u[1][1] = u[1][1] + sgn(R[i][i])*normx;
        normx = sqrt((u.t() * u)[1][1]);
        if (normx == 0) {
            cerr << "WARNING: Sistema singular." << endl;
            continue;
        }
        u = u*(1/normx);
        // R := HR; b := Hb
        R.subassign(i,n, i,n, R.submatrix(i,n, i,n) +
                              ((u*u.t())*(-2))*R.submatrix(i,n, i,n));
        b.subassign(i,n, 1,1, b.submatrix(i,n, 1,1) +
                              ((u*u.t())*(-2))*b.submatrix(i,n, 1,1));
    }
    // The system is now triangular, we go on to finding the solution
    // by back-substitution
    for (int i=n; i>=1; i--) {
        if (R[i][i] == 0) {
            cerr << "ERRO: Sistema singular." << endl;
            exit(SINGULAR_LINEAR_SYSTEM);
        }
        x[i][1] = b[i][1] - (R.submatrix(i,i, i+1,n)*x.submatrix(i+1,n, 1,1))[1][1] ;
        x[i][1] = x[i][1] / R[i][i];
    }
    return x;
}

void cppmatrix::fill_out_with_zeros(int rows, int cols) {
    n = rows;
    m = cols;
    for (int i=1; i<=n; i++) {
        if (this->count(i) == 0) // if this row is empty, then...
            for (int j=1; j<=m; j++)
                (*this)[i][j] = 0;
        else
            for (int j=1; j<=m; j++)
                if ( ((*this)[i]).count(j) == 0 )
                    (*this)[i][j] = 0;
    }
}

/* funcoes para a resolucao do circuito */

map<int, string> split(string str, int &i, char delim) {
    map<int, string> result;
    const char *cstr = str.c_str();
    for (i = 0; i >= 0; i++) {
        const char *begin = cstr;
        while (*cstr != delim && *cstr)
            cstr++;
        result[i] = string(begin, cstr);
        if (0 == *cstr++)
            break;
    }
    return result;
}

void cppmatrix::printMyself() {
    for (int i=1; i<=n; i++) {
        for (int j=1; j<=m; j++)
            printf("% 22.16Lf ", (*this)[i][j]);
        cout << endl;
    }
}

// construtor
element::element () {
    originNodeOrPositiveOutputNode              = -1;
    destinationNodeOrNegativeOutputNode         = -1;
    value                                       = -1;
    controledOriginNodeOrPositiveInputNode      = -1;
    controledDestinationNodeOrNegativeInputNode = -1;
    initialConditions                           = -1;
    pairsOfValues                               = "";
    parameter                                   = "";
    nocrtlPositive                              = -1;
    nocrtlNegative                              = -1;
    gon                                         = -1;
    goff                                        = -1;
    vref                                        = -1;
}

void element::printMyself() {
    cout << endl << "PRINTMYSELF()" << endl;
    cout << "  originNodeOrPositiveOutputNode              = " <<
            originNodeOrPositiveOutputNode << endl;
    cout << "  destinationNodeOrNegativeOutputNode         = " <<
            destinationNodeOrNegativeOutputNode << endl;
    cout << "  value                                       = " <<
            value << endl;
    cout << "  controledOriginNodeOrPositiveInputNode      = " <<
            controledOriginNodeOrPositiveInputNode << endl;
    cout << "  controledDestinationNodeOrNegativeInputNode = " <<
            controledDestinationNodeOrNegativeInputNode << endl;
    cout << "  initialConditions                           = " <<
            initialConditions << endl;
    cout << "  pairsOfValues                               = [" <<
            pairsOfValues << "]" << endl;
    cout << "  parameter                                   = [" <<
            parameter << "]" << endl;
    cout << "  nocrtlPositive                              = " <<
            nocrtlPositive << endl;
    cout << "  nocrtlNegative                              = " <<
            nocrtlNegative << endl;
    cout << "  gon                                         = " <<
            gon << endl;
    cout << "  goff                                        = " <<
            goff << endl;
    cout << "  vref                                        = " <<
            vref << endl;
}

/****************************************************************************/
/* Funcao responsavel por, a partir da leitura individual de cada linha do
 * arquivo dado, distribuir os elementos do circuito assim como seus
 * respectivos parametros no objeto list, pertencente a classe elementsList.
 */

void elementsList::getElement (string line) {

    string elementName;
    int qty_of_words;
    map<int, string> split_line = split(line, qty_of_words);

//     cout << endl << "INSIDE GETELEMENT" << endl << endl;

//     for (int i=0; i<=qty_of_words; i++)
//         cout << "  word N_" << i << " : [" << split_line[i] << "]" << endl;

    elementName = split_line[0];
    (*this)[elementName] = new element;
//     (*this)[elementName]->printMyself();

//     cout << endl << "  elementData=[" << elementData << "]; " <<
//             "elementName=[" << elementName << "];" << endl << endl;

//         cout << "  elementData=[" << elementData << "]; " <<
//                 "elementName=[" << elementName << "];" << endl;

    switch (elementName[0]) {

      case 'R':
        if (qty_of_words != 3) {
            cerr << "R<nome> <no1> <no2> <resistencia>" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->value = atoi (split_line[3].c_str());
        break;

      case 'I':
        if (qty_of_words < 3) {
            cerr << "I<nome> <no+> <no-> <parametros>" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->parameter = split_line[3];
        if (split_line[3].compare("DC") == 0) {
            if (qty_of_words != 4) {
                cerr << "I<nome> <no+> <no-> DC <valor>" << endl;
                exit(BAD_NETLIST);
            }
            (*this)[elementName]->value = atof(split_line[4].c_str());
        }
        else if (split_line[3].compare("SIN") == 0) {
            cerr << "NAO TEM FONTE SENOIDAL AINDA" << endl;
            exit(BAD_NETLIST);
        }
        else if (split_line[3].compare("PULSE") == 0) {
            cerr << "NAO TEM FONTE PULSE AINDA" << endl;
            exit(BAD_NETLIST);
        }
        else {
            cerr << "Este parametro \"" << split_line[3]
                 << "\" para a fonte independente nao existe." << endl;
            exit(BAD_NETLIST);
        }
        break;

      default:
        cout << "Element " << elementName << " not implemented (yet?)"
                << endl;
    }

        /*switch (i) {
            case 1:
                ( (*this)[elementName] )->originNodeOrPositiveOutputNode =
                                                     atoi (elementData.c_str());
                break;
            case 2:
                ( (*this)[elementName] )-> destinationNodeOrNegativeOutputNode =
                                                     atoi (elementData.c_str());
                break;
            case 3:
                switch (situation){
                    case (1): // passivo
                        ( (*this)[elementName] )->value =
                                                     atof (elementData.c_str());
                        break;
                    case (2): // fonte contralada (ou amp op)
                        ( (*this)[elementName] )->controledOriginNodeOrPositiveInputNode
                                                   = atoi (elementData.c_str());
                        break;
                    case (3): // resistor nao linear
                        ( (*this)[elementName] )->pairsOfValues = elementData;
                        break;
                    case (4): // chave
                        ( (*this)[elementName] )->nocrtlPositive =
                                                     atoi (elementData.c_str());
                    break;
                    case (5): // fonte independente
                        ( (*this)[elementName])->parameter = elementData;
                    break;
                }
                break;
            case 4:
                switch (situation){
                    case (1):
                        ( (*this)[elementName] )->initialConditions =
                                                     atof (elementData.c_str());
                        break;
                    case (2):
                        ( (*this)[elementName] )->controledDestinationNodeOrNegativeInputNode
                                                   = atoi (elementData.c_str());
                        break;
                    case (4):
                        ( (*this)[elementName] )->nocrtlNegative =
                                                     atoi (elementData.c_str());
                        break;
//                     case (5):
//                         ( (*this)[elementName] )->value =
//                                                      atof (elementData.c_str());
//                         break;
                }
                break;
            case 5:
                switch (situation){
                    case (2):
                        ( (*this)[elementName] )->value = atof (elementData.c_str());
                        break;
                    case (4):
                        ( (*this)[elementName] )->gon = atof (elementData.c_str());
                        break;
                }
                break;
            case 6:
                ( (*this)[elementName] )->goff = atof (elementData.c_str());
                break;
            case 7:
                ( (*this)[elementName] )->vref = atof (elementData.c_str());
                break;
        }*/

//         (*this)[elementName]->printMyself();

//         cout << endl << "  END LOOP" << endl << endl;

//     cout << "  END GETELEMENT();" << endl << endl;

}

/****************************************************************************/
/* Função responsável em retornar a quantidade de nos que o circuito
 * analisado possui.
 */
int elementsList::numberOfNodes() {

    int node = 0;

    map <string, element *> :: iterator auxiliar;

    for (auxiliar = this->begin(); auxiliar != this->end(); auxiliar++) {

        if ( ((auxiliar->second)->originNodeOrPositiveOutputNode) > node )
            node = ((auxiliar->second)->originNodeOrPositiveOutputNode);

        if ( ((auxiliar->second)->destinationNodeOrNegativeOutputNode) > node )
            node = ((auxiliar->second)->destinationNodeOrNegativeOutputNode);

        if ( ((auxiliar->second)->controledOriginNodeOrPositiveInputNode) > node )
            node = ((auxiliar->second)->controledOriginNodeOrPositiveInputNode);

        if ( ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode) > node )
            node = ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode);

    }

    return node;
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

/***************************************************************************************************/
/* Funcao responsavel em, a partir da identificacao de cada elemento
 * pertencente ao objeto list, construir sua respectica estampa e soma-la as
 * matrizes da analise nodal modificada. As matrizes serao os objetos matrix1
 * e matrix3, que pertencem a classe cppmatrix e representam,
 * respectivamente, as matrizes A e B do sistema A x = B.
 */
void elementsList::buildModifiedNodalMatrix
    (tensionAndCurrent& listToPrint, cppmatrix& matrix1, cppmatrix& matrix3) {

    /* O index inicialmente corresponde ao numero de nos do circuito. Ao longo
     * da funcao, ele é incrementado cada vez que se faz necessario o calculo
     * de uma nova corrente no circuito. Ao final da funcao, seu valor
     * corresponde a ordem da matrix A
    */
    int index = numberOfNodes();

    stringstream node;

    map <string, element *> :: iterator auxiliar;

    /**
    *** PERGUNTAR PRA MARCELLE POR QUE QUE TEM ESSE index++ AQUI
    *** (QUE TA COMENTADO AGORA, POIS TAVA DANDO ERRADO COM ELE)
    *** acho que era pra acrescentar o nó de "terra" (mas nao precisa)
    **/
    for (auxiliar = this->begin()/*, index++*/; auxiliar != this->end(); auxiliar++) {

        /* While responsavel pela construcao da estampa dos elementos */
        switch (auxiliar->first[0]) {

          case 'R': /* Resistor */
            if ((auxiliar->second)->originNodeOrPositiveOutputNode == 0) {
                if ((auxiliar->second)->destinationNodeOrNegativeOutputNode == 0)
                    break;
                matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
                        [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
                        += 1/((auxiliar->second)->value) ;
                break;
            }
            if ((auxiliar->second)->destinationNodeOrNegativeOutputNode == 0) {
                if ((auxiliar->second)->originNodeOrPositiveOutputNode == 0)
                    break;
                matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
                        [(auxiliar->second)->originNodeOrPositiveOutputNode]
                        += 1/((auxiliar->second)->value) ;
                break;
            }
            matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
                    [(auxiliar->second)->originNodeOrPositiveOutputNode]
                    += 1/((auxiliar->second)->value) ;

            matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
                    [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
                    += 1/((auxiliar->second)->value) ;

            matrix1 [(auxiliar->second)->originNodeOrPositiveOutputNode]
                    [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
                    += -1/((auxiliar->second)->value) ;

            matrix1 [(auxiliar->second)->destinationNodeOrNegativeOutputNode]
                    [(auxiliar->second)->originNodeOrPositiveOutputNode]
                    += -1/((auxiliar->second)->value) ;
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

            index++;
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

            /* Quando a fonte depende de uma corrente, a corrente sera calculada e sera impressa no arquivo
            * sendo referida por jab, onde a e o no de origem da corrente e b e o no de destino da corrente.
            */
            node << (auxiliar->second)->controledOriginNodeOrPositiveInputNode
                    << (auxiliar->second)->controledDestinationNodeOrNegativeInputNode;
            listToPrint[index] = 'j' + node.str();
                            /*	(locateCurrent (((auxiliar->second)->controledOriginNodeOrPositiveInputNode),
                                            ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode))); */
            index++;
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

            node << (auxiliar->second)->controledOriginNodeOrPositiveInputNode
                    << (auxiliar->second)->controledDestinationNodeOrNegativeInputNode;
            listToPrint[index] = 'j' + node.str();
                /*  (locateCurrent (((auxiliar->second)->controledOriginNodeOrPositiveInputNode),
                    ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode))); */
            index++;
            break;
          case 'I': /*Fonte de corrente*/
            if ((auxiliar->second)->originNodeOrPositiveOutputNode == 0) {
                if ((auxiliar->second)->destinationNodeOrNegativeOutputNode == 0)
                    break;
                matrix3 [(auxiliar->second)->destinationNodeOrNegativeOutputNode][1]
                        += -(auxiliar->second->value);
                break;
            }
            if ((auxiliar->second)->destinationNodeOrNegativeOutputNode == 0) {
                if ((auxiliar->second)->originNodeOrPositiveOutputNode == 0)
                    break;
                matrix3 [(auxiliar->second)->originNodeOrPositiveOutputNode][1]
                        += (auxiliar->second->value);
                break;
            }
            matrix3 [(auxiliar->second)->originNodeOrPositiveOutputNode][1]
                    += ((auxiliar->second)->value);
            matrix3 [(auxiliar->second)->destinationNodeOrNegativeOutputNode][1]
                    += -((auxiliar->second)->value);
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
            index++;
            break;
          case 'L': /*indutor*/
            break;
        }
    }

    // preenche com zeros as entradas da matriz que nao foram usadas
    // (sabendo que a ordem eh `index')
    matrix1.fill_out_with_zeros(index, index);
    matrix3.fill_out_with_zeros(index, 1);

}

/****************************************************************************************************/

/* Funcao responsavel em imprimir o resultado em um arquivo */
void elementsList :: printResult (char* argv[], tensionAndCurrent listToPrint, cppmatrix matrix2) {

    map <int, string> :: iterator list;
    ofstream answerFile;
    int auxiliar = numberOfNodes();

    answerFile.open (*argv[1] + "_answer");

    /* Imprimindo a primeira linha no arquivo */
    answerFile << " t" << endl;
    for (int aux = auxiliar; aux == 0; aux++) {
        answerFile << " " << aux << endl;
        }
    for (list = listToPrint.begin(); list != listToPrint.end(); list++){
        answerFile << " " << list->second << endl;
        auxiliar++;
    }
    answerFile << "\n";

    /* Imprimindo o resto do arquivo */

    for (int aux = 1; aux == auxiliar; aux ++){
        answerFile << " " << matrix2[aux][0] ;
    }
    answerFile << "\n";

}
