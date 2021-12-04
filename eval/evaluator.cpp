#include "evaluator.h"

using namespace std;

/* === evalFunc ==================================
評価関数
最小化問題にするように設計する
evaluationPoint = f(x)
=============================================== */

double eval::evalFunc( vector<double> input ){
   double evaluationPoint = 1000000;

   evaluationPoint = eval::Sphare( input );

   return evaluationPoint;
}

vector<double> eval::evalFunc( vector<vector<double> > input ){
   vector<vector<double> > vtmp = calVec::transpose( input );
   int n = vtmp.size();
   vector<double> eval( n, 1000 );

   for( int i = 0; i < n; ++i ){
      eval[i] = eval::evalFunc( vtmp[i] );
   }

   return eval;
}



/* === benchmark ==================================
================================================ */

double eval::Sphare( vector<double> input ){
   double evaluationPoint;
   int dimension = input.size();

   for( int i = 0; i < dimension; ++i ){
      evaluationPoint += pow( input[i], 2 );
   }

   return evaluationPoint;
}

double eval::k_tablet( vector<double> input ){
   double evaluationPoint;
   int dimension = input.size();
   int k = dimension / 4;

   for( int i = 0; i < k; ++i ){
      evaluationPoint += pow( input[i], 2 );
   }
   for( int i = k; i < dimension; ++i ){
      evaluationPoint += pow( 100 * input[i], 2 );
   }

   return evaluationPoint;
}

double eval::Ellipsoid( vector<double> input ){
   double evaluationPoint, tmpValue;
   int dimension = input.size();

   dimension = dimension - 1;
   for( int i = 0; i <= dimension; ++i ){
      tmpValue = pow( 1000, ( i / dimension ) ) * input[i];
      evaluationPoint += ( tmpValue * tmpValue );
   }

   return evaluationPoint;
}

double eval::Rosenbrock( vector<double> input ){
   double evaluationPoint, tmpValue;
   int dimension = input.size();

   dimension = dimension - 1;
   for( int i = 0; i < dimension; ++i ){
      tmpValue = pow( input[i+1] - input[i] * input[i], 2 ) + pow( input[i] - 1, 2 );
      evaluationPoint += ( 100 * tmpValue );
   }

   return evaluationPoint;
}

double eval::Ackley( vector<double> input ){
   double evaluationPoint, tmpValue;
   int dimension = input.size();

   evaluationPoint = 20 + exp( 1 );
   for( int i = 0; i < dimension; ++i ){
      tmpValue += input[i] * input[i];
   }
   tmpValue = -0.2 * sqrt( tmpValue / dimension );
   evaluationPoint -= 20 * exp( tmpValue );
   tmpValue = 0;
   for( int i = 0; i < dimension; ++i ){
      tmpValue += cos( 2 * M_PI * input[i] );
   }
   evaluationPoint -= exp( tmpValue / dimension );

   return evaluationPoint;
}

double eval::Bohachevsky( vector<double> input ){
   double evaluationPoint, tmpValue;
   int dimension = input.size();

   dimension = dimension - 1;
   for( int i = 0; i < dimension; ++i ){
      tmpValue = 0.7 + pow( input[i], 2 ) + ( 2 * pow( input[i+1], 2 ) );
      evaluationPoint += tmpValue - ( 0.3 * cos( 3 * M_PI * input[i] ) ) - ( 0.4 * cos( 4 * M_PI * input[i] ) );
   }

   return evaluationPoint;
}

double eval::Schaffer( vector<double> input ){
   double evaluationPoint, tmpValue;
   int dimension = input.size();

   dimension = dimension - 1;
   for( int i = 0; i < dimension; ++i ){
      tmpValue = pow( input[i], 2 ) + pow( input[i+1], 2 );
      evaluationPoint += pow( tmpValue, 0.25 ) * ( pow( sin( 50 * pow( tmpValue, 0.1 ) ), 2 ) + 1.0 );
   }

   return evaluationPoint;
}

double eval::Rastrigin( vector<double> input ){
   double evaluationPoint;
   int dimension = input.size();

   evaluationPoint = 10 * dimension;
   for( int i = 0; i < dimension; ++i ){
      evaluationPoint += pow( input[i], 2 ) - 10 * cos( 2 * M_PI * input[i] );
   }

   return evaluationPoint;
}
