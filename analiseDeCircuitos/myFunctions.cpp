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

cppmatrix::cppmatrix (int rows, int columns) {
    n=rows; m=columns;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            matrix[i][j] = 0;
}

cppmatrix cppmatrix::operator+ (cppmatrix parcela) {
    cppmatrix resultado(n,m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            resultado.matrix[i][j] = matrix[i][j] + parcela.matrix[i][j];
    return resultado;
}

cppmatrix cppmatrix::operator* (cppmatrix fator) {
    cppmatrix resultado(n,fator.m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= fator.m; j++)
            for (int k = 1; k <= m; k++)
                resultado.matrix[i][j] += matrix[i][k] + fator.matrix[k][j];
    return resultado;
}

cppmatrix cppmatrix::operator* (float factor) {
    cppmatrix resultado(n,m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            resultado.matrix[i][j] = factor * matrix[i][j];
    return resultado;
}

cppmatrix& cppmatrix::operator[] (int) {
    ;
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
            printf("% 6.6f ", matrix[i][j]);
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

void modifiedMatrix::printMyself() {
    for (modifiedMatrix::iterator i = this->begin(); i != this->end(); i++) {
        cout << "LINHA " << i->first << ": ";
        for (map<int, float>::iterator j = i->second.begin(); j != i->second.end(); j++) {
            cout << j->first << ":" << j->second << " ";
        }
        cout << endl;
    }
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
 * e matrix3, que pertencem a classe modifiedMatrix e representam,
 * respectivamente, as matrizes A e B do sistema A x = B.
 */
void elementsList::buildModifiedNodalMatrix
  (int& matrixOrder, tensionAndCurrent& listToPrint,
   modifiedMatrix& matrix1, modifiedMatrix& matrix3) {

    /* O index inicialmente corresponde ao numero de nos do circuito. Ao longo
     * da funcao, ele é incrementado cada vez que se faz necessario o calculo
     * de uma nova corrente no circuito. Ao final da funcao, seu valor
     * corresponde a ordem da matrix A
    */
    int index = numberOfNodes();

    cout << "index: " << index << endl;

    stringstream node;

    map <string, element *> :: iterator auxiliar;

    for (auxiliar = this->begin(), index++; auxiliar != this->end(); auxiliar++) {

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
                matrix3 [(auxiliar->second)->destinationNodeOrNegativeOutputNode][0]
                        += -(auxiliar->second->value);
                break;
            }
            if ((auxiliar->second)->destinationNodeOrNegativeOutputNode == 0) {
                if ((auxiliar->second)->originNodeOrPositiveOutputNode == 0)
                    break;
                matrix3 [(auxiliar->second)->originNodeOrPositiveOutputNode][0]
                        += (auxiliar->second->value);
                break;
            }
            matrix3 [(auxiliar->second)->originNodeOrPositiveOutputNode][0]
                    += ((auxiliar->second)->value);
            matrix3 [(auxiliar->second)->destinationNodeOrNegativeOutputNode][0]
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

    matrixOrder = index;

}



/**************************************************************************************************/

/* Funcao responsavel em resolver o sistema  A x = B */
void modifiedMatrix::solveMatrixSystem (int order, modifiedMatrix matrix1, modifiedMatrix& matrix2, modifiedMatrix matrix3) {

    int i, j;

    /* s/10/order/g */
    float A[10][10];
    float x[10];
    float B[10];

        /* Construindo a Matriz A */
    for (i=1; i == order; i++){
            for (j=1; j == order; j++){
                A[i][j] = 0;
                A[i][j] += matrix1 [i][j];
            }
    }

    /* Construindo a Matriz B */
    for (i=1; i== order; i++){
        B[i] = 0;
        B[i] += matrix3[i][0];
    }


    /* Resolvendo o sistema e achando x */


    /* Copiando a matriz x para matrix2 */
    for (i=1; i== order; i++){
        matrix2[i][0] = x[i];
        }
}

/****************************************************************************************************/

/* Funcao responsavel em imprimir o resultado em um arquivo */
void elementsList :: printResult (char* argv[], tensionAndCurrent listToPrint, modifiedMatrix matrix2) {

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
