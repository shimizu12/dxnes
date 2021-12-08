#ifndef CALCURATEVECTOR_H
#define CALCURATEVECTOR_H

#include <iostream>

#include <vector>
#include <random>
#include <Eigen/Core>


using namespace std;

class calVec{
private:

public:
   static vector<double> randn( int dim );
   static vector<vector<double> > randn( int dim, int lambda );

   static vector<vector<double> > repmat( vector<double> x, int horizon, int vertical );
   static vector<vector<double> > repmat( vector<vector<double> > x, int horizon, int vertical );

   static vector<vector<double> > eye( int x );

   static double norm( vector<double> vec );

   static vector<double> add( vector<double> v, double x );
   static vector<double> add( vector<double> v1, vector<double> v2 );
   static vector<vector<double> > add( vector<vector<double> > v, double x );
   static vector<vector<double> > add( vector<vector<double> > v1, vector<vector<double> > v2 );

   static double times( vector<double> v1, vector<double> v2 );
   static vector<double> times( vector<double> v1, double x );
   static vector<double> times( vector<vector<double> > v1, vector<double> v2 );
   static vector<vector<double> > times( vector<vector<double> > v1, double x );
   static vector<vector<double> > times( vector<vector<double> > v1, vector<vector<double> > v2 );

   static vector<vector<double> > transpose( vector<double> v1 );
   static vector<vector<double> > transpose( vector<vector<double> > v1 );

   static double trace( vector<vector<double> > v );

   static vector<vector<double> >  cat( vector<vector<double> > x, vector<vector<double> > y, int mode );

   static double sum( vector<double> v );
   static vector<double> sum( vector<vector<double> > v, int dim );

   static vector<double> sqrtVec( vector<double> v );

   static vector<vector<double> > expm( vector<vector<double> > vec );

   static vector<double> expVec( vector<double> v );
   static vector<vector<double> > expVec( vector<vector<double> > v );

   static vector<vector<double> > inverseVec( vector<vector<double> > vec );

};

#endif // SAMPLE_H
