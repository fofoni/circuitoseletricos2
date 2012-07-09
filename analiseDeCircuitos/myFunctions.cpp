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
cppmatrix cppmatrix::operator- (cppmatrix parcela) {
    cppmatrix resultado;
    resultado.initialize(n,m);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            resultado[i][j] = (*this)[i][j] - parcela[i][j];
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
            cerr << "         Norma da sub-coluna: " << normx << endl;
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
        if (this->count(i) == 0) // if this row is empty...
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
        else if (split_line[3].compare("PULSE") == 0) {
            if (qty_of_args != 11) {
                cerr << "Erro na leitura do arquivo."
                << " Maneira correta de especificar a fonte pulsada:"
                << endl << elementName << " <no+> <no-> PULSE <amplit. 1> <amplt. 2> <atraso> <tempo subida> <tempo ligada> <período> <num. de ciclos>"
                << endl;
                exit(BAD_NETLIST);
            }
            (*this)[elementName]->ampl = strtold(split_line[4].c_str(), NULL);
            (*this)[elementName]->ampl2 = strtold(split_line[5].c_str(), NULL);
            (*this)[elementName]->atraso = strtold(split_line[6].c_str(), NULL);
            (*this)[elementName]->t_rise = strtold(split_line[7].c_str(), NULL);
            (*this)[elementName]->t_fall = strtold(split_line[8].c_str(), NULL);
            (*this)[elementName]->t_ligada = strtold(split_line[9].c_str(), NULL);
            (*this)[elementName]->periodo = strtold(split_line[10].c_str(), NULL);
            (*this)[elementName]->num_de_ciclos = strtold(split_line[11].c_str(), NULL);
        }
        else if (split_line[3].compare("SIN") == 0) {
            if (qty_of_args != 10) {
                cerr << "Erro na leitura do arquivo."
                << " Maneira correta de especificar a fonte senoidal:"
                << endl << elementName << " <no+> <no-> SIN <nível DC> <amplitude> <freq(Hz)> <atraso> <atenuação> <ângulo> <num. de ciclos>"
                << endl;
                exit(BAD_NETLIST);
            }
            (*this)[elementName]->dc_level = strtold(split_line[4].c_str(), NULL);
            (*this)[elementName]->ampl = strtold(split_line[5].c_str(), NULL);
            (*this)[elementName]->freq = strtold(split_line[6].c_str(), NULL);
            (*this)[elementName]->atraso = strtold(split_line[7].c_str(), NULL);
            (*this)[elementName]->atenuacao = strtold(split_line[8].c_str(), NULL);
            (*this)[elementName]->angulo = strtold(split_line[9].c_str(), NULL);
            (*this)[elementName]->num_de_ciclos = strtold(split_line[10].c_str(), NULL);
        }
        else {
            cerr << "Este parametro \"" << split_line[3]
                 << "\" para a fonte independente nao existe." << endl;
            exit(BAD_NETLIST);
        }
        break;

      case 'O':
        if (qty_of_args != 4) {
            cerr << "Erro na leitura do arquivo."
                 << " Maneira correta de especificar um amp.op.:"
                 << endl << elementName << " <out+> <out-> <in+> <in->" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->controledOriginNodeOrPositiveInputNode = atoi (split_line[3].c_str());
        (*this)[elementName]->controledDestinationNodeOrNegativeInputNode = atoi (split_line[4].c_str());
        break;

      case 'C': case 'L':
        if ((qty_of_args < 3) or (qty_of_args > 4)) {
            cerr << "Erro na leitura do arquivo."
                 << " Maneira correta de especificar um elemento reativo:"
                 << endl << elementName
                 << " <nó1> <nó2> <valor> [IC=condição inicial]" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi (split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi (split_line[2].c_str());
        (*this)[elementName]->value = strtold(split_line[3].c_str(), NULL);
        if (qty_of_args == 4) {
            if (split_line[4].size() <= 3) {
                cerr << "Erro na leitura do arquivo."
                     << " Maneira correta de especificar um elemento reativo:"
                     << endl << elementName
                     << " <nó1> <nó2> <valor> [IC=condição inicial]" << endl;
                exit(BAD_NETLIST);
            }
            (*this)[elementName]->initialConditions = strtold(split_line[4].substr(3).c_str(), NULL);
//            cerr << "READ IC=" << (*this)[elementName]->initialConditions << endl;
        }
        break;

      case 'N':
        if (qty_of_args != 10) {
            cerr << "Erro na leitura do arquivo."
                << " Maneira correta de especificar um resistor nao-linear:"
                << endl << elementName
                << " <nó1> <nó2> <quatro pares de valores vi ji>" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi(split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi(split_line[2].c_str());
        (*this)[elementName]->tensoes[0] = strtold(split_line[3].c_str(), NULL);
        (*this)[elementName]->correntes[0] = strtold(split_line[4].c_str(), NULL);
        (*this)[elementName]->tensoes[1] = strtold(split_line[5].c_str(), NULL);
        (*this)[elementName]->correntes[1] = strtold(split_line[6].c_str(), NULL);
        (*this)[elementName]->tensoes[2] = strtold(split_line[7].c_str(), NULL);
        (*this)[elementName]->correntes[2] = strtold(split_line[8].c_str(), NULL);
        (*this)[elementName]->tensoes[3] = strtold(split_line[9].c_str(), NULL);
        (*this)[elementName]->correntes[3] = strtold(split_line[10].c_str(), NULL);
        for (int k=0;k<3;k++) {
            (*this)[elementName]->condutancia[k] = ((*this)[elementName]->correntes[k+1] - (*this)[elementName]->correntes[k])/
                                                   ((*this)[elementName]->tensoes[k+1] - (*this)[elementName]->tensoes[k]);
            (*this)[elementName]->fonte_corrente[k] = (*this)[elementName]->correntes[k] -
                                                      (*this)[elementName]->condutancia[k] * (*this)[elementName]->tensoes[k];
        }

//         cout << "=============== RES NÃO LINEAR =================" << endl;
//         cout << (*this)[elementName]->originNodeOrPositiveOutputNode; cout << endl;
//         cout << (*this)[elementName]->destinationNodeOrNegativeOutputNode; cout << endl;
//         cout << (*this)[elementName]->tensoes[0]; cout << endl;
//         cout << (*this)[elementName]->correntes[0]; cout << endl;
//         cout << (*this)[elementName]->tensoes[1]; cout << endl;
//         cout << (*this)[elementName]->correntes[1]; cout << endl;
//         cout << (*this)[elementName]->tensoes[2]; cout << endl;
//         cout << (*this)[elementName]->correntes[2]; cout << endl;
//         cout << (*this)[elementName]->tensoes[3]; cout << endl;
//         cout << (*this)[elementName]->correntes[3]; cout << endl;
//         for (int k=0;k<3;k++) {
//             cout << (*this)[elementName]->condutancia[k]; cout << endl;
//             cout << (*this)[elementName]->fonte_corrente[k]; cout << endl;
//         }
//         cout << "=============== RES NÃO LINEAR =================" << endl;

        break;

      case '$':
        if (qty_of_args < 4 or qty_of_args > 7) {
            cerr << "Erro na leitura do arquivo."
                << " Maneira correta de especificar uma chave:"
                << endl << elementName
                << " <nó1> <nó2> <nó ctrl+> <nó ctrl-> [<gon> <gof> <vref>]" << endl;
            exit(BAD_NETLIST);
        }
        (*this)[elementName]->originNodeOrPositiveOutputNode = atoi(split_line[1].c_str());
        (*this)[elementName]->destinationNodeOrNegativeOutputNode = atoi(split_line[2].c_str());
        (*this)[elementName]->nocrtlPositive = atoi(split_line[3].c_str());
        (*this)[elementName]->nocrtlNegative = atoi(split_line[4].c_str());
        switch (qty_of_args) {
          case 4:
            (*this)[elementName]->gon = 1e20;
            (*this)[elementName]->goff = 0;
            (*this)[elementName]->vref = 0;
            break;
          case 5:
            (*this)[elementName]->gon = strtold(split_line[5].c_str(), NULL);
            (*this)[elementName]->goff = 0;
            (*this)[elementName]->vref = 0;
            break;
          case 6:
            (*this)[elementName]->gon = strtold(split_line[5].c_str(), NULL);
            (*this)[elementName]->goff = strtold(split_line[6].c_str(), NULL);
            (*this)[elementName]->vref = 0;
            break;
          case 7:
            (*this)[elementName]->gon = strtold(split_line[5].c_str(), NULL);
            (*this)[elementName]->goff = strtold(split_line[6].c_str(), NULL);
            (*this)[elementName]->vref = strtold(split_line[7].c_str(), NULL);
            break;
        }
        break;

      default:
        cout << "Element " << elementName << " not implemented." << endl;
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
    (tensionAndCurrent& listToPrint, cppmatrix& matrix1, cppmatrix matrix2, cppmatrix& matrix3,
    capacitor_inductor& reactiveElements, long double passo, int gear_order,
    int UIC, long double time_inst) {

    /* O index inicialmente corresponde ao numero de nos do circuito. Ao longo
     * da funcao, ele é incrementado cada vez que se faz necessario o calculo
     * de uma nova corrente no circuito. Ao final da funcao, seu valor
     * corresponde a ordem da matrix A
    */
    int index = numberOfNodes()+1;

    stringstream node;

    map <string, element *> :: iterator auxiliar;
//     capacitor_inductor :: iterator capacitorInductor = reactiveElements.begin();

    matrix1.clear();
    matrix3.clear();

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
          case '$': /* chave */
            {
                long double g_chave;
                long double queda_de_tensao = matrix2[auxiliar->second->nocrtlPositive][1]-matrix2[auxiliar->second->nocrtlNegative][1];
                if (queda_de_tensao >= auxiliar->second->vref) g_chave = auxiliar->second->gon;
                else g_chave = auxiliar->second->goff;
                matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                        [auxiliar->second->originNodeOrPositiveOutputNode]
                        += g_chave;
                matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        += g_chave;
                matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                        [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        += -g_chave;
                matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        [auxiliar->second->originNodeOrPositiveOutputNode]
                        += -g_chave;
            }
            break;
          case 'N':
            {
              long double queda_de_tensao = matrix2[auxiliar->second->originNodeOrPositiveOutputNode][1]-matrix2[auxiliar->second->destinationNodeOrNegativeOutputNode][1];
              int regiao;
              if (queda_de_tensao < auxiliar->second->tensoes[1]) regiao=0;
              else if (queda_de_tensao > auxiliar->second->tensoes[2]) regiao=2;
              else regiao=1;
//               cout << "no+=" << matrix2[auxiliar->second->originNodeOrPositiveOutputNode][1] << " no-=" << matrix2[auxiliar->second->destinationNodeOrNegativeOutputNode][1] << " regiao:"<<regiao<<endl;
              matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                      [auxiliar->second->originNodeOrPositiveOutputNode]
                      += auxiliar->second->condutancia[regiao];
              matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                      [auxiliar->second->destinationNodeOrNegativeOutputNode]
                      += auxiliar->second->condutancia[regiao];
              matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                      [auxiliar->second->destinationNodeOrNegativeOutputNode]
                      += - auxiliar->second->condutancia[regiao];
              matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                      [auxiliar->second->originNodeOrPositiveOutputNode]
                      += - auxiliar->second->condutancia[regiao];
              matrix3 [auxiliar->second->originNodeOrPositiveOutputNode][1]
                      += - auxiliar->second->fonte_corrente[regiao];
              matrix3 [auxiliar->second->destinationNodeOrNegativeOutputNode][1]
                      += auxiliar->second->fonte_corrente[regiao];
            }
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
            {
                long double valor_da_fonte;
                if (auxiliar->second->parameter.compare("DC") == 0)
                    valor_da_fonte = auxiliar->second->value;
                else if (auxiliar->second->parameter.compare("SIN") == 0) {
                    if ((time_inst < auxiliar->second->atraso) or (time_inst > auxiliar->second->atraso + auxiliar->second->num_de_ciclos/auxiliar->second->freq))
                        valor_da_fonte = auxiliar->second->dc_level + auxiliar->second->ampl * sin(auxiliar->second->angulo * M_PI / 180.0 );
                    else
                        valor_da_fonte = auxiliar->second->dc_level + auxiliar->second->ampl * exp(-(time_inst-auxiliar->second->atraso)*auxiliar->second->atenuacao) * sin(2 * M_PI * auxiliar->second->freq * (time_inst-auxiliar->second->atraso) + auxiliar->second->angulo * M_PI / 180.0);
                }
                else if (auxiliar->second->parameter.compare("PULSE") == 0) {
                    if ((time_inst < auxiliar->second->atraso) or (time_inst > auxiliar->second->atraso + auxiliar->second->num_de_ciclos*auxiliar->second->periodo))
                        valor_da_fonte = auxiliar->second->ampl;
                    else {
                        long double t_mod_per = (time_inst-auxiliar->second->atraso) - auxiliar->second->periodo * floor((time_inst-auxiliar->second->atraso)/auxiliar->second->periodo);
                        if (t_mod_per < auxiliar->second->t_rise)
                            valor_da_fonte = auxiliar->second->ampl + (auxiliar->second->ampl2 - auxiliar->second->ampl)*t_mod_per/auxiliar->second->t_rise;
                        else if (t_mod_per <= auxiliar->second->t_rise + auxiliar->second->t_ligada)
                            valor_da_fonte = auxiliar->second->ampl2;
                        else if (t_mod_per < auxiliar->second->t_rise + auxiliar->second->t_ligada + auxiliar->second->t_fall)
                            valor_da_fonte = auxiliar->second->ampl2 - (auxiliar->second->ampl2 - auxiliar->second->ampl)*(t_mod_per-(auxiliar->second->t_rise + auxiliar->second->t_ligada))/auxiliar->second->t_fall;
                        else
                            valor_da_fonte = auxiliar->second->ampl;
                    }
                }
                else exit(UNKNOWN_ERROR);
                matrix3 [auxiliar->second->originNodeOrPositiveOutputNode][1]
                        += - valor_da_fonte;
                matrix3 [auxiliar->second->destinationNodeOrNegativeOutputNode][1]
                        += valor_da_fonte;
            }
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
            if (auxiliar->second->parameter.compare("DC") == 0)
                matrix3 [index][1] += auxiliar->second->value;
            else if (auxiliar->second->parameter.compare("SIN") == 0) {
                if ((time_inst < auxiliar->second->atraso) or (time_inst > auxiliar->second->atraso + auxiliar->second->num_de_ciclos/auxiliar->second->freq))
                    matrix3 [index][1] += auxiliar->second->dc_level + auxiliar->second->ampl * sin(auxiliar->second->angulo * M_PI / 180.0 );
                else
                    matrix3 [index][1] += auxiliar->second->dc_level + auxiliar->second->ampl * exp(-(time_inst-auxiliar->second->atraso)*auxiliar->second->atenuacao) * sin(2 * M_PI * auxiliar->second->freq * (time_inst-auxiliar->second->atraso) + auxiliar->second->angulo * M_PI / 180.0);
            }
            else if (auxiliar->second->parameter.compare("PULSE") == 0) {
                //cout << time_inst << "; atraso=" << auxiliar->second->atraso << " numciclos=" << auxiliar->second->num_de_ciclos << " per=" << auxiliar->second->periodo << endl;
                //cout << ;
                if ((time_inst < auxiliar->second->atraso) or (time_inst > auxiliar->second->atraso + auxiliar->second->num_de_ciclos*auxiliar->second->periodo))
                    matrix3 [index][1] += auxiliar->second->ampl;
                else {
                    long double t_mod_per = (time_inst-auxiliar->second->atraso) - auxiliar->second->periodo * floor((time_inst-auxiliar->second->atraso)/auxiliar->second->periodo);
                    if (t_mod_per < auxiliar->second->t_rise)
                        matrix3 [index][1] += auxiliar->second->ampl + (auxiliar->second->ampl2 - auxiliar->second->ampl)*t_mod_per/auxiliar->second->t_rise;
                    else if (t_mod_per <= auxiliar->second->t_rise + auxiliar->second->t_ligada)
                        matrix3 [index][1] += auxiliar->second->ampl2;
                    else if (t_mod_per < auxiliar->second->t_rise + auxiliar->second->t_ligada + auxiliar->second->t_fall)
                        matrix3 [index][1] += auxiliar->second->ampl2 - (auxiliar->second->ampl2 - auxiliar->second->ampl)*(t_mod_per-(auxiliar->second->t_rise + auxiliar->second->t_ligada))/auxiliar->second->t_fall;
                    else
                        matrix3 [index][1] += auxiliar->second->ampl;
                }
            }
            else exit(UNKNOWN_ERROR);
            listToPrint[index] = "j" + (auxiliar->first);
            index++;
            break;
          case 'O': /* Operacional */
            matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                    [index]
                    += 1;
            matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                    [index]
                    += -1;
            matrix1 [index]
                    [auxiliar->second->controledOriginNodeOrPositiveInputNode]
                    += -1;
            matrix1 [index]
                    [auxiliar->second->controledDestinationNodeOrNegativeInputNode]
                    += 1;
            listToPrint[index] = "j" + auxiliar->first;
            index++;
            break;
          case 'L': case 'C': /* indutor e capacitor*/
            gearMethod (auxiliar->first, reactiveElements,
                        passo, gear_order, UIC);
            if (auxiliar->first[0] == 'C') {

                matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                        [auxiliar->second->originNodeOrPositiveOutputNode]
                        += 1/auxiliar->second->impedance;
                matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        += 1/auxiliar->second->impedance;
                matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                        [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        += -1/auxiliar->second->impedance;
                matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        [auxiliar->second->originNodeOrPositiveOutputNode]
                        += -1/auxiliar->second->impedance;

                reactiveElements[auxiliar->first][8] = auxiliar->second->originNodeOrPositiveOutputNode;
                reactiveElements[auxiliar->first][9] = auxiliar->second->destinationNodeOrNegativeOutputNode;
                matrix3 [auxiliar->second->originNodeOrPositiveOutputNode][1]
                        += (auxiliar->second->reactiveValue);
                matrix3 [auxiliar->second->destinationNodeOrNegativeOutputNode][1]
                        += -(auxiliar->second->reactiveValue);

            }
            else {

                reactiveElements[auxiliar->first][8] = index;

                matrix1 [index]
                        [auxiliar->second->originNodeOrPositiveOutputNode]
                        += -1;
                matrix1 [index]
                        [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        += 1;
                matrix1 [index]
                        [index]
                        += auxiliar->second->impedance;
                matrix1 [auxiliar->second->originNodeOrPositiveOutputNode]
                        [index] += 1;
                matrix1 [auxiliar->second->destinationNodeOrNegativeOutputNode]
                        [index] += -1;
                matrix3 [index][1] += auxiliar->second->reactiveValue;

                listToPrint[index] = 'j' + auxiliar->first;
                index++;
            }

            break;
        }
    }

    // preenche com zeros as entradas da matriz que nao foram usadas
    // (sabendo que a ordem eh `index-1')
    matrix1.fill_out_with_zeros(index-1, index-1);
    matrix3.fill_out_with_zeros(index-1, 1);

}

void print_state(capacitor_inductor reactiveElements) {
    char new_str[13];
    string header = "";
    cout << "================================================" << endl;
    cout << "==  ";
    for (capacitor_inductor::iterator it = reactiveElements.begin();
         it != reactiveElements.end();
         it++) {
        sprintf(new_str, " %11s", it->first.c_str());
        header = header + new_str;
    }
    cout << header << endl;
    for (int k=0; k<8; k++) {
        cout << "==  ";
        for (capacitor_inductor::iterator it = reactiveElements.begin();
             it != reactiveElements.end(); it++) {
            sprintf(new_str, " % 11.4Lg", it->second[k]);
            cout << new_str;
        }
        cout << endl;
    }
    cout << "== == == == == == == == == ==" << endl << "==  ";
    for (capacitor_inductor::iterator it = reactiveElements.begin();
         it != reactiveElements.end(); it++) {
        sprintf(new_str, " % 11.4Lg", it->second[8]);
        cout << new_str;
    }
    cout << endl << "==  ";
    for (capacitor_inductor::iterator it = reactiveElements.begin();
         it != reactiveElements.end(); it++) {
        if (it->first[0] == 'C')
            sprintf(new_str, " % 11.4Lg", it->second[9]);
        else
            sprintf(new_str, "            ");
        cout << new_str;
    }
    cout << endl << "================================================" << endl;
}

/*****************************************************************************/
/* Funcao responsavel em substituir indutor/capacitor pelo respectivo modelo */

void elementsList::gearMethod(string element_name,
                              capacitor_inductor reactiveElements,
                              long double passo, int gear_order, int UIC) {

    long double sourceValue = 0;

    // gear_coef[8] is the term that multiplies the time step.
    long double gear_coef[9];

    switch (gear_order) {
      case 1:
        gear_coef[8] = 1.0;
        gear_coef[0] = 1.0;
        gear_coef[1] = 0.0;
        gear_coef[2] = 0.0;
        gear_coef[3] = 0.0;
        gear_coef[4] = 0.0;
        gear_coef[5] = 0.0;
        gear_coef[6] = 0.0;
        gear_coef[7] = 0.0;
        break;
      case 2:
        gear_coef[8] =  ((long double)(2)) / ((long double)(3));
        gear_coef[0] =  ((long double)(4)) / ((long double)(3));
        gear_coef[1] = -((long double)(1)) / ((long double)(3));
        gear_coef[2] = 0.0;
        gear_coef[3] = 0.0;
        gear_coef[4] = 0.0;
        gear_coef[5] = 0.0;
        gear_coef[6] = 0.0;
        gear_coef[7] = 0.0;
        break;
      case 3:
        gear_coef[8] =  ((long double)( 6)) / ((long double)(11));
        gear_coef[0] =  ((long double)(18)) / ((long double)(11));
        gear_coef[1] = -((long double)( 9)) / ((long double)(11));
        gear_coef[2] =  ((long double)( 2)) / ((long double)(11));
        gear_coef[3] = 0.0;
        gear_coef[4] = 0.0;
        gear_coef[5] = 0.0;
        gear_coef[6] = 0.0;
        gear_coef[7] = 0.0;
        break;
      case 4:
        gear_coef[8] =  ((long double)(12)) / ((long double)(25));
        gear_coef[0] =  ((long double)(48)) / ((long double)(25));
        gear_coef[1] = -((long double)(36)) / ((long double)(25));
        gear_coef[2] =  ((long double)(16)) / ((long double)(25));
        gear_coef[3] = -((long double)( 3)) / ((long double)(25));
        gear_coef[4] = 0.0;
        gear_coef[5] = 0.0;
        gear_coef[6] = 0.0;
        gear_coef[7] = 0.0;
        break;
      case 5:
        gear_coef[8] =  ((long double)( 60)) / ((long double)(137));
        gear_coef[0] =  ((long double)(300)) / ((long double)(137));
        gear_coef[1] = -((long double)(300)) / ((long double)(137));
        gear_coef[2] =  ((long double)(200)) / ((long double)(137));
        gear_coef[3] = -((long double)( 75)) / ((long double)(137));
        gear_coef[4] =  ((long double)( 12)) / ((long double)(137));
        gear_coef[5] = 0.0;
        gear_coef[6] = 0.0;
        gear_coef[7] = 0.0;
        break;
      case 6:
        gear_coef[8] =  ((long double)( 60)) / ((long double)(147));
        gear_coef[0] =  ((long double)(360)) / ((long double)(147));
        gear_coef[1] = -((long double)(450)) / ((long double)(147));
        gear_coef[2] =  ((long double)(400)) / ((long double)(147));
        gear_coef[3] = -((long double)(225)) / ((long double)(147));
        gear_coef[4] =  ((long double)( 72)) / ((long double)(147));
        gear_coef[5] = -((long double)( 10)) / ((long double)(147));
        gear_coef[6] = 0.0;
        gear_coef[7] = 0.0;
        break;
      case 7:
        gear_coef[8] =  0.3856749311292203;
        gear_coef[0] =  2.6997245179004445;
        gear_coef[1] = -4.0495867768486979;
        gear_coef[2] =  4.4995408631693161;
        gear_coef[3] = -3.3746556473786473;
        gear_coef[4] =  1.6198347107421940;
        gear_coef[5] = -0.4499540863173528;
        gear_coef[6] =  0.0550964187327434;
        gear_coef[7] =  0.0;
        break;
      case 8:
        gear_coef[8] =  0.3679369250676245;
        gear_coef[0] =  2.9434954005987457;
        gear_coef[1] = -5.1511169509341315;
        gear_coef[2] =  6.8681559346235712;
        gear_coef[3] = -6.4388961887049154;
        gear_coef[4] =  4.1208935607669224;
        gear_coef[5] = -1.7170389836515461;
        gear_coef[6] =  0.4204993429348310;
        gear_coef[7] = -0.0459921156334772;
        break;
    }

//     cout << "==================GEAR===================" << endl <<
//             gear_coef[0] << endl << gear_coef[1] << endl <<
//             gear_coef[2] << endl << gear_coef[3] << endl <<
//             gear_coef[4] << endl << gear_coef[5] << endl <<
//             gear_coef[6] << endl << gear_coef[7] << endl <<
//             element_name << endl <<
//             "=========================================" << endl;

    for (int i=0; i < gear_order; i++)
        sourceValue += gear_coef[i] * reactiveElements[element_name][i];

    if (element_name[0] == 'L') {
        (*this)[element_name]->impedance = ((*this)[element_name]->value)/(gear_coef[8] * passo);
        (*this)[element_name]->reactiveValue = sourceValue * (*this)[element_name]->impedance;
//         cout << "IMP:  " << (*this)[element_name]->impedance << endl;
//         cout << "REAC: " << (*this)[element_name]->reactiveValue << endl;
    }
    else {
        (*this)[element_name]->impedance = gear_coef[8] * passo / ((*this)[element_name]->value);
        (*this)[element_name]->reactiveValue = sourceValue / (*this)[element_name]->impedance;
//         cout << "IMP:  " << (*this)[element_name]->impedance << endl;
//         cout << "REAC: " << (*this)[element_name]->reactiveValue << endl;
    }

}


/****************************************************************************************************/

/* Função responsavel em imprimir o resultado em um arquivo */
// void elementsList :: printResult (string netlist_filename, tensionAndCurrent listToPrint, cppmatrix matrix2) {
//
//     map <int, string> :: iterator list;
//     ofstream answerFile;
//     int auxiliar = numberOfNodes();
//
//     answerFile.open (netlist_filename + "_answer.txt");
//
//     /* Imprimindo a primeira linha no arquivo */
//     answerFile << " t" << endl;
//     for (int aux = auxiliar; aux == 0; aux++) {
//         answerFile << " " << aux << endl;
//         }
//     for (list = listToPrint.begin(); list != listToPrint.end(); list++){
//         answerFile << " " << list->second << endl;
//         auxiliar++;
//     }
//     answerFile << "\n";
//
//     /* Imprimindo o resto do arquivo */
//
//     for (int aux = 1; aux == auxiliar; aux ++){
//         answerFile << " " << matrix2[aux][0] ;
//     }
//     answerFile << "\n";
//
//     answerFile.close();
// }
