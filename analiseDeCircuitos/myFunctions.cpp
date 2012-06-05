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

element::printMyself() {
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

void elementsList::getElement (string line, int situation) {

    int auxiliar = 0;
    int aux = 0;
    int quantityOfSpaces = 0;
    string elementData;

    cout << endl << "INSIDE GETELEMENT" << endl << endl;

    do {
        cout << "  auxiliar=" << auxiliar << "; aux=" << aux << "; qtySpc=" <<
                quantityOfSpaces << "; elementData=[" << elementData << "];" <<
                endl;
        cout << "  while ( line[auxiliar] != ' ') auxiliar++;" << endl;
        while (line [auxiliar] != ' ')
            auxiliar++;
        cout << "  auxiliar=" << auxiliar << "; aux=" << aux << "; qtySpc=" <<
                quantityOfSpaces << "; elementData=[" << elementData << "];" <<
                endl;

        cout << "  elementData = line.substr(" << aux << ", " << (auxiliar-aux) << ");" << endl;
        elementData = line.substr(aux, auxiliar-aux);
        cout << "  auxiliar=" << auxiliar << "; aux=" << aux << "; qtySpc=" <<
                quantityOfSpaces << "; elementData=[" << elementData << "];" <<
                endl << endl;

        switch (quantityOfSpaces) {
            case 0:
                /* Novo elemento é inserido ao objeto */
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

/****************************************************************************/
/* Função responsável em retornar a quantidade de nos que o circuito
 * analisado possui.
 */
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

/***************************************************************************************************/
/* Funcao responsavel em, a partir da identificacao de cada elemento
 * pertencente ao objeto list, construir sua respectica estampa e soma-la as
 * matrizes da analise nodal modificada. As matrizes serao os objetos matrix1
 * e matrix3, que pertencem a classe modifiedMatrix e representam,
 * respectivamente, as matrizes A e B do sistema A x = B.
 */
void elementsList::buildModifiedNodalMatrix (int *matrixOrder, tensionAndCurrent listToPrint, modifiedMatrix matrix1, modifiedMatrix matrix3){

    /* O index inicialmente corresponde ao numero de nos do circuito. Ao longo
     * da funcao, ele é incrementado cada vez que se faz necessario o calculo
     * de uma nova corrente no circuito. Ao final da funcao, seu valor
     * corresponde a ordem da matrix A
    */
    int index = numberOfNodes();
    char *node;

    map <string, element *> :: iterator auxiliar;

    for (auxiliar = this->begin(), index++; auxiliar != this->end(); auxiliar++) {

        /* While responsavel pela construcao da estampa dos elementos */
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
                    sprintf (node, "%i%i", (auxiliar->second)->controledOriginNodeOrPositiveInputNode, (auxiliar->second)->controledDestinationNodeOrNegativeInputNode);
                    listToPrint[index] = 'j' + node;
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

                    sprintf (node, "%i%i", (auxiliar->second)->controledOriginNodeOrPositiveInputNode, (auxiliar->second)->controledDestinationNodeOrNegativeInputNode);
                    listToPrint[index] = 'j' + node;
                                    /*	(locateCurrent (((auxiliar->second)->controledOriginNodeOrPositiveInputNode),
                                        ((auxiliar->second)->controledDestinationNodeOrNegativeInputNode))); */
                    index++;

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
                    index++;
                break;
            case 'L': /*indutor*/
                break;
        }
    }

    *matrixOrder = index;

}



/**************************************************************************************************/

/* Funcao responsavel em resolver o sistema  A x = B */
void modifiedMatrix::solveMatrixSystem (int * matrixOrder, modifiedMatrix matrix1, modifiedMatrix matrix2, modifiedMatrix matrix3) {

    int order = *matrixOrder;
    int i, j;

    float A[order][order];
    float x[order];
    float B[order];

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
