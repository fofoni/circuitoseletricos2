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
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

#ifndef     CIRCUITSANALYSES_H_
#define     CIRCUITSANALYSES_H_

#define     FILE_IS_NOT_OPEN        1
#define     BAD_NETLIST             2
#define     NON_SQUARE_MATRIX_QR    3
#define     SINGULAR_LINEAR_SYSTEM  4

#define     tensionAndCurrent   map <int, string>

map<int, string> split(string str, int &i, char delim = ' ');

/* Classe que possui os parametros de todos os
   possiveis elementos dentro de um netlist */
class element {
    public:
        int     originNodeOrPositiveOutputNode;
        int     destinationNodeOrNegativeOutputNode;
        long double   value;
        int     controledOriginNodeOrPositiveInputNode;
        int     controledDestinationNodeOrNegativeInputNode;
        long double   initialConditions;
        string  pairsOfValues;
        string  parameter;
        int     nocrtlPositive;
        int     nocrtlNegative;
        long double   gon;
        long double   goff;
        long double   vref;
        void    printMyself();
        element();
    };

/* Classe que vai conter as matrizes da analise modificada,
   formada a partir das estampas dos elementos */
                               /* row, column, content */
class modifiedMatrix : public map<int, map<int, long double> > {
  public:
    void solveMatrixSystem (int, modifiedMatrix, modifiedMatrix&, modifiedMatrix);
    void printMyself();
};


/* Container cuja chave é o nome do elemento do circuito
   formada pela classe element */
class elementsList : public map<string, element*> {
  public:
    void    getElement (string);
    int     numberOfNodes();
    void    buildModifiedNodalMatrix (int&, tensionAndCurrent&,
                                      modifiedMatrix&, modifiedMatrix&);
    void    printResult (char**, tensionAndCurrent, modifiedMatrix);
    /* string locateCurrent (int, int);*/
};

class cppmatrix : public map<int, map<int, long double> > {
  public:
    int n,m;
    cppmatrix operator + (cppmatrix);
    cppmatrix operator * (cppmatrix);
    cppmatrix operator * (long double);
    cppmatrix t(); // transpose
    cppmatrix submatrix(int, int, int, int);
    void subassign(int,int, int,int, cppmatrix);
    cppmatrix solve(cppmatrix b);
    void printMyself();
    void initialize(int, int);
    void make_id(int);
};

#endif
