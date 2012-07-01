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
        long double   reactiveValue;
        long double   impedance;
        void    printMyself();
        element();
    };

/* Classe que vai conter as matrizes da analise modificada,
   formada a partir das estampas dos elementos */
                               /* row, column, content */
class cppmatrix : public map<int, map<int, long double> > {
  public:
    int n,m;
    map <int, bool> colunas_removidas;
    map <int, bool> linhas_removidas;
    cppmatrix operator + (cppmatrix);
    cppmatrix operator * (cppmatrix);
    cppmatrix operator * (long double);
    cppmatrix t(); // transpose
    cppmatrix submatrix(int, int, int, int);
    void subassign(int,int, int,int, cppmatrix);
    cppmatrix solveMatrixSystem(cppmatrix);
    void printMyself();
    void initialize(int, int);
    void make_id(int);
    void fill_out_with_zeros(int, int);
    cppmatrix remove_opamp_entries();
};

/* Classe que contem os valores das correntes anteriores dos indutores e das tensoes
 * anteriores dos capacitores nas 8 primeiras posicoes. Nas 2 ultimas posicoes contem
 * os nos de origem e destino dos capacitores e na penultima posicao o valor da linha
 * na qual a corrente dos indutors encontra-se.
 */
class capacitor_inductor : public map<string, long double [10]> {
};


/* Container cuja chave eh o nome do elemento do circuito
   formada pela classe element */
class elementsList : public map<string, element*> {
  public:
    void    getElement (string);
    int     numberOfNodes();
    void	gearMethod (iterator, capacitor_inductor, long double, int, int);
    void    buildModifiedNodalMatrix (tensionAndCurrent&,
                                      cppmatrix&, cppmatrix&,
                                      capacitor_inductor, long double, int, int);
    void    printResult (char**, tensionAndCurrent, cppmatrix);
    /* string locateCurrent (int, int);*/
};

#endif
