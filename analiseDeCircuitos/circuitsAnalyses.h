/* Universidade Federal do Rio de Janeiro
 * Escola Politécnica
 * Circuitos Elétricos II / 1.2012
 * Prof. Antônio Carlos Moreirão de Queiroz
 * Alunos: Marcelle de Souza Campos
 *         Pedro Angelo Medeiros Fonini
 * Programa de análise de circuitos no tempo para estudar os métodos de Gear
 */

#include <iostream>
#include <istream>
#include <fstream>
#include <cerrno>
#include <string>
#include <map>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#ifndef     CIRCUITSANALYSES_H_
#define     CIRCUITSANALYSES_H_

#define     FILE_IS_NOT_OPEN    1
#define     ARMA_DONT_USE_BLAS



/* Classe que possui os parametros de todos os
   possiveis elementos dentro de um netlist */
class element {
    public:
        int     originNodeOrPositiveOutputNode;
        int     destinationNodeOrNegativeOutputNode;
        float   value;
        int     controledOriginNodeOrPositiveInputNode;
        int     controledDestinationNodeOrNegativeInputNode;
        float   inicialConditions;
        string  pairsOfValues;
        string  parameter;
        int     nocrtlPositive;
        int     nocrtlNegative;
        float   gon;
        float   goff;
        float   vref;
    };

/* Classe que vai conter as matrizes da analise modificada,
   formada a partir das estampas dos elementos */
                               /* row, column, content */
class modifiedMatrix : public map<int, map<int, float> > {
  public:
    void solveMatrixSystem (int*, modifiedMatrix, modifiedMatrix, modifiedMatrix);
};

/* Classe que numera todas as variaveis do circuitos, sendo elas tensoes no nos e devidas correntes */
class tensionAndCurrent : public map <int, string>{
};

/* Container cuja chave e o nome do elemento do circuito formada pela classe element */
class elementsList : public map<string, element*> {
    public:
    void    getElement (string, int);
    int     numberOfNodes();
    void    buildModifiedNodalMatrix (int*, tensionAndCurrent, modifiedMatrix, modifiedMatrix);
    void    printResult (char**, tensionAndCurrent, modifiedMatrix);
    /* string locateCurrent (int, int);*/

};

#endif
