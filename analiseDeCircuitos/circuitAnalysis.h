/* Universidade Federal do Rio de Janeiro
 * Escola Politécnica
 * Circuitos Elétricos II / 1.2012
 * Prof. Antônio Carlos Moreirão de Queiroz
 * Alunos: Marcelle de Souza Campos
 *         Pedro Angelo Medeiros Fonini
 * Programa de análise de circuitos no tempo para estudar os métodos de Gear
 */

// Se, durante o newton raphson, o sistema ficar singular, faz o
// solveMatrixSystem retornar um bool por parametro dizendo que deu singular,
// e ai tem que testar outro chute inicial

#ifndef     CIRCUITSANALYSES_H_
#define     CIRCUITSANALYSES_H_

#define _XOPEN_SOURCE 601

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

#define     FILE_IS_NOT_OPEN        1
#define     BAD_NETLIST             2
#define     NON_SQUARE_MATRIX_QR    3
#define     SINGULAR_LINEAR_SYSTEM  4
#define     UNKNOWN_ERROR           5

#define CHANGE_STRTOLD
#define OUTPUT_MATLAB

#ifdef CHANGE_STRTOLD
# define strtold(x,y) ((long double)(strtod((x),(y))))
#endif

#ifndef M_PI
# define M_PI ((long double)(3.14159265358979))
#endif

#define     tensionAndCurrent       map <int, string>

/* Classe que contem os valores das correntes anteriores dos indutores e das tensoes
 * anteriores dos capacitores nas 8 primeiras posicoes. Nas 2 ultimas posicoes contem
 * os nos de origem e destino dos capacitores e na penultima posicao o valor da linha
 * na qual a corrente dos indutors encontra-se.
 */
#define     capacitor_inductor      map<string, map<int, long double> >

void print_state(capacitor_inductor reactiveElements);

map<int, string> split(string str, int &i, char delim = ' ');

string rm_filename_extension(string filename);

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
        string  parameter;
        int     nocrtlPositive;
        int     nocrtlNegative;
        long double   gon;
        long double   goff;
        long double   vref;

        // res nao linear
        long double   tensoes[4];
        long double   correntes[4];
        long double   condutancia[3];
        long double   fonte_corrente[3];

        // elementos reativos
        long double   reactiveValue;
        long double   impedance;

        void    printMyself();
        element();

        // fonte independente
        long double dc_level;
        long double ampl;
        long double freq;
        long double atraso;
        long double atenuacao;
        long double angulo;
        long double num_de_ciclos;

        // fonte pulse
        long double ampl2;
        long double t_rise;
        long double t_fall;
        long double t_ligada;
        long double periodo;
    };

/* Classe que vai conter as matrizes da analise modificada,
   formada a partir das estampas dos elementos */
                          /* row,  column, content */
class cppmatrix : public map<int, map<int, long double> > {
  public:
    int n,m;
    cppmatrix operator + (cppmatrix);
    cppmatrix operator - (cppmatrix);
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
};

/* Container cuja chave eh o nome do elemento do circuito
   formada pela classe element */
class elementsList : public map<string, element*> {
  public:
    void    getElement (string);
    int     numberOfNodes();
    void	gearMethod (string element_name,
                        capacitor_inductor, long double, int, int);
    void    buildModifiedNodalMatrix (tensionAndCurrent&,
                                      cppmatrix&, cppmatrix, cppmatrix&,
                                      capacitor_inductor&,
                                      long double, int, int, long double);
    //void    printResult (char**, tensionAndCurrent, cppmatrix);
    /* string locateCurrent (int, int);*/
};

#endif
