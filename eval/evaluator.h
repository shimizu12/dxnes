#include <iostream>
#include <cmath>
#include <vector>

#include "../calVec/calcurateVector.h"

#ifndef EVALUATOR_H
#define EVALUATOR_H

using namespace std;

class eval{
private:
   // static int dimension;
   // static double evaluationPoint;

public:
   static double evalFunc( vector<double> input );    //
   static vector<double> evalFunc( vector<vector<double> > input );

   static double Sphare( vector<double> input );      // Sphare
   static double k_tablet( vector<double> input );    // k-tablet(k=dim/4)
   static double Ellipsoid( vector<double> input );   // Ellipsoid
   static double Rosenbrock( vector<double> input );  // Rosenbrock
   static double Ackley( vector<double> input );      // Ackley
   static double Bohachevsky( vector<double> input ); // Bohachevsky
   static double Schaffer( vector<double> input );    // Schaffer
   static double Rastrigin( vector<double> input );   // Rastrigin

};

#endif
