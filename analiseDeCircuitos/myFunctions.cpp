/* Universidade Federal do Rio de Janeiro
 * Escola Politécnica
 * Circuitos Elétricos II / 1.2012
 * Prof. Antônio Carlos Moreirão de Queiroz
 * Alunos: Marcelle de Souza Campos
 *         Pedro Angelo Medeiros Fonini
 * Programa de análise de circuitos no tempo para estudar os métodos de Gear
 */

#include "circuitAnalysis.h"

#define GET_STRING(ch) (                                                    \
    (ch=='E') ? ("E<nome> <nóV+> <nóV-> <nóv+> <nóv-> <Av>") : (            \
    (ch=='F') ? ("F<nome> <nóI+> <nóI-> <nói+> <nói-> <Ai>") : (            \
    (ch=='G') ? ("G<nome> <nóI+> <nóI-> <nóv+> <nóv-> <Gm>") : (            \
    (ch=='H') ? ("H<nome> <nóV+> <nóV-> <nói+> <nói-> <Rm>") :              \
    "Unknown error.")))                                                     \
)

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
        if (normx < 1e-9) {
            cerr << "WARNING: Sistema potencialmente singular." << endl;
            cerr << "Norma da sub-coluna: " << normx << endl;
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
    int qty_of_args;
    map<int, string> split_line = split(line, qty_of_args);

    elementName = split_line[0];
    (*this)[elementName] = new element;

    switch (elementName[0]) {

      case 'R':
        if (qty_of_args != 3) {
            cerr << "Erro na leitura do arquivo."
                 << " Maneira correta de especificar um resistor:"
                 << endl << "R<nome> <no1> <no2> <resistencia>" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->value = strtold (split_line[3].c_str(), NULL);
        break;

      case 'E':
      case 'F':
      case 'G':
      case 'H':
        if (qty_of_args != 5) {
            cerr << "Erro na leitura do arquivo."
                 << " Maneira correta de especificar a fonte controlada:"
                 << endl << GET_STRING(elementName[0]) << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->controledOriginNodeOrPositiveInputNode = atoi (split_line[3].c_str());
        (*this)[elementName]->controledDestinationNodeOrNegativeInputNode = atoi (split_line[4].c_str());
        (*this)[elementName]->value = strtold (split_line[5].c_str(), NULL);
        break;

      case 'I':
      case 'V':
        if (qty_of_args < 3) {
            cerr << "Erro na leitura do arquivo."
                 << " Maneira correta de especificar a fonte independente:"
                 << endl << elementName[0] << "<nome> <no+> <no-> <parametros>"
                 << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->parameter = split_line[3];
        if (split_line[3].compare("DC") == 0) {
            if (qty_of_args != 4) {
                cerr << elementName[0] << "<nome> <no+> <no-> DC <valor>"
                     << endl;
                exit(BAD_NETLIST);
            }
            (*this)[elementName]->value = strtold(split_line[4].c_str(), NULL);
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
        cout << "Element " << elementName << " not implemented (yet?)" << endl;
        exit(BAD_NETLIST);

    }

}

/****************************************************************************/
/* Função responsável em retornar a quantidade de nos que o circuito
 * analisado possui.
 */
int elementsList::numberOfNodes() {

    int node = 0;

    map <string, element *> :: iterator auxiliar;

    for (auxiliar = this->begin(); auxiliar != this->end(); auxiliar++) {

        if ( (auxiliar->second->originNodeOrPositiveOutputNode) > node )
            node = (auxiliar->second->originNodeOrPositiveOutputNode);

        if ( (auxiliar->second->destinationNodeOrNegativeOutputNode) > node )
            node = (auxiliar->second->destinationNodeOrNegativeOutputNode);

        if ( (auxiliar->second->controledOriginNodeOrPositiveInputNode) > node )
            node = (auxiliar->second->controledOriginNodeOrPositiveInputNode);

        if ( (auxiliar->second->controledDestinationNodeOrNegativeInputNode) > node )
            node = (auxiliar->second->controledDestinationNodeOrNegativeInputNode);

    }

    return node;
}

/*
string elementsList::locateCurrent (int node1, int node2){

    string name = NULL;

    map <string, element *> :: iterator auxiliar;

    for (auxiliar = this->begin(); auxiliar != this->end(); auxiliar++) {

            if ( ((auxiliar->second-> originNodeOrPositiveOutputNode) == node1 ) &&
                ((auxiliar->second-> destinationNodeOrNegativeOutputNode) == node2 ))
                name = auxiliar->first;

    }
    return (name);
}
*/

/******************************************************************************/
/* Funcao responsavel por, a partir da identificacao de cada elemento
 * pertencente ao objeto list, construir sua respectica estampa e soma-la as
 * matrizes da analise nodal modificada. As matrizes serao os objetos matrix1
 * e matrix3, que pertencem a classe cppmatrix e representam,
 * respectivamente, as matrizes A e B do sistema A x = B.
 */
void elementsList::buildModifiedNodalMatrix
    (tensionAndCurrent& listToPrint, cppmatrix& matrix1, cppmatrix& matrix3,
    capacitor_inductor reactiveElements, long double passo, int gear_order,
    int UIC) {

    /* O index inicialmente corresponde ao numero de nos do circuito. Ao longo
     * da funcao, ele é incrementado cada vez que se faz necessario o calculo
     * de uma nova corrente no circuito. Ao final da funcao, seu valor
     * corresponde a ordem da matrix A
    */
    int index = numberOfNodes()+1;

    stringstream node;

    map <string, element *> :: iterator auxiliar;
    capacitor_inductor :: iterator capacitorInductor = reactiveElements.begin();

    for (auxiliar = this->begin(); auxiliar != this->end(); auxiliar++) {

        /* While responsavel pela construcao da estampa dos elementos */
        switch (auxiliar->first[0]) {

          // Não precisa verificar se os nós são zero ou não. A função que
          // resolve o sistema linear ignora a linha zero e a coluna zero.

          case 'R': /* Resistor */
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [auxiliar->second->originNodeOrPositiveOutputNode]
                    += 1/(auxiliar->second->value) ;

            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    += 1/(auxiliar->second->value) ;

            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    += -1/(auxiliar->second->value) ;

            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [auxiliar->second->originNodeOrPositiveOutputNode]
                    += -1/(auxiliar->second->value) ;
            break;
          case 'E': /* Fonte de tensao controlada a tensao */
            matrix1 [index]
                    [auxiliar->second->originNodeOrPositiveOutputNode]
                    += 1;
            matrix1 [index]
                    [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    += -1;
            matrix1 [index]
                    [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    += - auxiliar->second->value;
            matrix1 [index]
                    [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    += auxiliar->second->value;
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [index]
                    += 1;
            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [index]
                    += -1;

            listToPrint[index] = "j" + auxiliar->first;

            index++;
            break;
          case 'F': /* Fonte de corrente controlada a corrente */
            matrix1 [index]
                    [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    += 1;
            matrix1 [index]
                    [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    += -1;
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [index]
                    += auxiliar->second->value;
            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [index]
                    += - auxiliar->second->value;
            matrix1 [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    [index]
                    += 1;
            matrix1 [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    [index]
                    += -1;

            /* Quando a fonte depende de uma corrente, a corrente sera
               calculada e sera impressa no arquivo sendo referida por jc_d,
               onde `c' é o no de origem da corrente e `d' é o no de destino.
            */
            {
                char helper_str[10];
                sprintf(helper_str, "j%d_%d",
                        auxiliar->second->controledOriginNodeOrPositiveInputNode,
                        auxiliar->second->controledDestinationNodeOrNegativeInputNode
                        );
                listToPrint[index] = helper_str;
            }
            index++;
            break;
          case 'G': /* Fonte de corrente controlada a tensao */
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    += auxiliar->second->value;

            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    += - auxiliar->second->value;

            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    += - auxiliar->second->value;

            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    += auxiliar->second->value;
            break;
          case 'H': /* Fonte de tensao controlada a corrente */
            // `index' refers to the curr. in the short-circuit
            // `index+1' refers to the current in the transresistance
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [index+1] += 1;
            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [index+1] += -1;
            matrix1 [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    [index] += 1;
            matrix1 [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    [index] += -1;

            matrix1 [index]
                    [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    += 1;
            matrix1 [index]
                    [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    += -1;

            matrix1 [index+1]
                    [auxiliar->second->originNodeOrPositiveOutputNode]
                    += 1;
            matrix1 [index+1]
                    [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    += -1;
            matrix1 [index+1]
                    [index]
                    += - auxiliar->second->value;

            {
                char helper_str[10];
                sprintf(helper_str, "j%d_%d",
                        auxiliar->second->controledOriginNodeOrPositiveInputNode,
                        auxiliar->second->controledDestinationNodeOrNegativeInputNode
                        );
                listToPrint[index] = helper_str;
            }
            index++;

            listToPrint[index] = "j" + auxiliar->first;
            index++;

            break;
          case 'I': /* Fonte de corrente */
            matrix3 [auxiliar->second->originNodeOrPositiveOutputNode][1]
                    += - (auxiliar->second->value);
            matrix3 [auxiliar->second->destinationNodeOrNegativeOutputNode][1]
                    += (auxiliar->second->value);
            break;
          case 'V': /* Fonte de tensao */
            matrix1 [index]
                    [auxiliar->second->originNodeOrPositiveOutputNode]
                    += 1;
            matrix1 [index]
                    [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    += -1;
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [index]
                    += 1;
            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [index]
                    += -1;
            matrix3 [index][1]
                    += auxiliar->second->value;
            listToPrint[index] = "j" + (auxiliar->first);
            index++;
            break;
          case 'L': case 'C': /* indutor e capacitor*/
        	  /*capacitorInductor->first = auxiliar ->first; // essa linha ta dando erro
        	  this->gearMethod (auxiliar, reactiveElements, passo, gear_order, UIC);
        	  if (auxiliar->second->originNodeOrPositiveOutputNode == 0) {
        	       if (auxiliar->second->destinationNodeOrNegativeOutputNode == 0)
        	       break;
        	       matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
        	               [auxiliar->second->destinationNodeOrNegativeOutputNode]
        	               += 1/(auxiliar->second->impedance) ;
        	       break;
        	  }
        	  if (auxiliar->second->destinationNodeOrNegativeOutputNode == 0) {
        	      if (auxiliar->second->originNodeOrPositiveOutputNode == 0)
        	      break;
        	      matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
        	      [auxiliar->second->originNodeOrPositiveOutputNode]
        	      += 1/(auxiliar->second->impedance) ;
        	      break;
        	  }
        	  matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
        	          [auxiliar->second->originNodeOrPositiveOutputNode]
        	          += 1/(auxiliar->second->impedance) ;

        	  matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
        	          [auxiliar->second->destinationNodeOrNegativeOutputNode]
        	          += 1/(auxiliar->second->impedance) ;

        	  matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
        	          [auxiliar->second->destinationNodeOrNegativeOutputNode]
        	          += -1/(auxiliar->second->impedance) ;

        	  matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
        	          [auxiliar->second->originNodeOrPositiveOutputNode]
        	          += -1/(auxiliar->second->impedance) ;
        	  if (auxiliar->first[0] == 'C') {
        		  capacitorInductor->second[8] = auxiliar->second->originNodeOrPositiveOutputNode;
        		  capacitorInductor->second[9] = auxiliar->second->destinationNodeOrNegativeOutputNode;
        		  if (auxiliar->second->originNodeOrPositiveOutputNode == 0) {
        			  if (auxiliar->second->destinationNodeOrNegativeOutputNode == 0)
        			  break;
        			  matrix3 [auxiliar->second->destinationNodeOrNegativeOutputNode][1]
        			          += -(auxiliar->second->reactiveValue);
        	          break;
        	      }
        	      if (auxiliar->second->destinationNodeOrNegativeOutputNode == 0) {
        	    	  if (auxiliar->second->originNodeOrPositiveOutputNode == 0)
        	          break;
        	          matrix3 [auxiliar->second->originNodeOrPositiveOutputNode][1]
        	                  += (auxiliar->second->reactiveValue);
        	         break;
        	      }
        	      matrix3 [auxiliar->second->originNodeOrPositiveOutputNode][1]
        	               += (auxiliar->second->reactiveValue);
        	      matrix3 [auxiliar->second->destinationNodeOrNegativeOutputNode][1]
        	               += -(auxiliar->second->reactiveValue);
        	  }

        	  else{
        		  capacitorInductor->second[8] = index;

        		  matrix1 [index]
        		          [auxiliar->second->originNodeOrPositiveOutputNode] += -1;

        		  matrix1 [index]
        		          [auxiliar->second->destinationNodeOrNegativeOutputNode] += +1;

        		  matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
        		           [index] += 1;

        		  matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
        		          [index] += -1;

        		  matrix3 [index][1] += -(auxiliar->second->reactiveValue);

        		  listToPrint[index] = 'j' + (auxiliar->first);
        		  index++;
        	  }*/

        	  break;
        }
    }

    // preenche com zeros as entradas da matriz que nao foram usadas
    // (sabendo que a ordem eh `index-1')
    matrix1.fill_out_with_zeros(index-1, index-1);
    matrix3.fill_out_with_zeros(index-1, 1);

}

/*****************************************************************************************************/
/* Funcao responsavel em substituir indutor/capacitor pelo respectivo modelo */

void elementsList :: gearMethod (iterator elementPosition, capacitor_inductor reactiveElements,
								 long double passo, int gear_order, int UIC) {

	map <string, element *> :: iterator auxiliar = elementPosition;
	map <string, long double[10]> :: iterator element = reactiveElements.begin();

	cppmatrix matrixA, matrixB, matrixC;
    long double fontValue = 0;

	for (int i=0; i <8; i++){
		matrixC[i][0] = 0;
		for (int j=0; j<8; j++){
			matrixA[i][j] = 0;
		}
	}

	for (int i=0; i == gear_order; i++){
		matrixC [i][0] = 1;
		for (int j=0; j < gear_order; j++){
			matrixA [i][j]= (-j)^(gear_order);
		}
		matrixA [i][gear_order] = gear_order;
	}

	while (auxiliar->first != element->first)
		element++;

	matrixB = matrixA.solveMatrixSystem(matrixC);

	for (int i=0; i < gear_order; i++)
		fontValue+= matrixB[i][0] * element->second[i];

	auxiliar->second->impedance =  matrixB[gear_order][0] * passo * 1/(auxiliar->second->value);

	if (auxiliar->first[0] == 'L' )
		auxiliar->second->reactiveValue = fontValue * auxiliar->second->impedance;
	else
		auxiliar->second->reactiveValue = fontValue / auxiliar->second->impedance;

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

    answerFile.close();
}
