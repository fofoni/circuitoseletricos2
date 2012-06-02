/* Universidade Federal do Rio de Janeiro
   Escola Politécnica
   Circuitos Elétricos II / 1.2012
   Prof. Antônio Carlos Moreirão de Queiroz
   Alunos: Marcelle de Souza Campos
   		   Pedro Angelo Medeiros Fonini
   Programa de análise de circuitos no tempo para estudar os métodos de Gear
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

#ifndef MYFUNCTIONS_H_
#define MYFUNCTIONS_H_

#define     FILE_IS_NOT_OPEN    1

class tensionAndCurrentName : public map <int, string>{
	public:
};

class element {
 	public:
 		int 	originNodeOrPositiveOutputNode;
 		int 	destinationNodeOrNegativeOutputNode;
 		float	value;
 		int 	controledOriginNodeOrPositiveInputNode;
 		int		controledDestinationNodeOrNegativeInputNode;
 		float 	inicialConditions;
 		string	pairsOfValues;
 		string  parameter;
		int		nocrtlPositive;
 		int  	nocrtlNegative;
 		float 	gon;
 		float	goff;
 		float	vref;
 	};

class elementsList : public map<string, element*> {
    public:
	   void	getElement (string, int);
	   int numberOfNodes();
	   void buildModifiedNodalMatrix ();
	   string locateCurrent (int, int);

};

class modifiedMatrix : public map<int, map<int, float> > {   /* line, columns, content */
	public:
		void solveMatrixSystem ();
};



#endif
