#include "calcurateVector.h"

using namespace std;
using namespace Eigen;

/* --- randn -------------------------------------------------------------------
正規化個体を生成
正規分布N(0, 1)に従うランダム値を dim*lambda 個だけ生成する

input:
int dim        : 列数
int lambda     : 行数

output:
vector<vector<double> > normalDistribution   :
( vector<double> normalDistribution           : ) */
vector<vector<double> > calVec::randn( int dim, int lambda ){
   vector<vector<double> > normalDistribution( dim, vector<double>( lambda, 0 ) );
   for( int i = 0; i < dim; ++i ){
      normalDistribution[i] = randn( lambda );
   }
   return normalDistribution;
}
vector<double> calVec::randn( int dim ){
   vector<double> normalDistribution( dim, 0 );
   double ave = 0.0;
   double dispersion = 1.0;

   random_device seed_gen;
   default_random_engine engine( seed_gen() );
   normal_distribution<> dist( ave, dispersion );

   for( int i = 0; i < dim; ++i ){
      normalDistribution[i] = dist( engine );
   }

   return normalDistribution;
}



/* --- repmat ------------------------------------------------------------------
matlab.repmatの再現
ただ matlab とは第二引数と第三引数が逆になっている

input:
vector<vector<double> > x       :

output:
vector<vector<double> > values  :  */
vector<vector<double> > calVec::repmat( vector<vector<double> > x, int horizon, int vertical ){
   vector<vector <double> > values = x;

   --vertical;
   for( int i = 0; i < vertical; ++i ){
      for( int j = 0; j < x.size(); ++j ){
         for( int k = 0; k < x[j].size(); ++k ){
            values[j].push_back( x[j][k] );
         }
      }
   }

   int num = values.size();

   --horizon;
   for( int i = 0; i < horizon; ++i ){
      for( int j = 0; j < num; ++j ){
         values.push_back( values[j] );
      }
   }

   return values;
}
vector<vector<double> > calVec::repmat( vector<double> x, int horizon, int vertical ){
   vector<vector <double> > values( x.size(), vector<double>(1, 0) );
   int n = x.size();
   for( int i = 0; i < n; ++i ){
      values[i][0] = x[i];
   }
   values = repmat( values, horizon, vertical );
   return values;
}



/* --- eye ---------------------------------------------------------------------
dim*dim単位行列を生成する

input:
int dim                         : 単位行列の次元数

output:
vector<vector<double> > value   : 生成した単位行列 */
vector<vector<double> > calVec::eye( int x ){
   vector<vector<double> > values;
   //行列のリサイズ
   values.resize(x);
   for( int i = 0; i < x; ++i ){
      values[i].resize(x);
   }
   //単位行列となるように代入
   for( int i = 0; i < x; ++i ){
      for( int j = 0; j < x; ++j ){
         if( i == j ){
            values[i][j] = 1;
         }else{
            values[i][j] = 0;
         }
      }
   }
   return values;
}



/* --- norm --------------------------------------------------------------------
ベクトルのノルムを算出

input:
vector<double> x      : ノルム計算したいベクトル

output:
double norm_value     :  ベクトルのpノルム */
double calVec::norm( vector<double> vec ){
   int n = vec.size();
   double norm_value = 0;
   VectorXd M( n );
   for( int i = 0; i < n; ++i ){
      M(i) = vec[i];
   }

   norm_value = M.norm();
   // std::cout << "norm_value = " << norm_value << '\n';
   return norm_value;
}


/* ---  add(1) -----------------------------------------------------------------
行列+行列を計算する

input:
vector<vector<double> > v1      : 計算前のベクトル1
vector<vector<double> > v2      : 計算前のベクトル2

output:
vector<double> values           : 計算後のベクトル  */
vector<vector<double> > calVec::add( vector<vector<double> > v1, vector<vector<double> > v2 ){
   if( ( v1.size() != v2.size() )
   &&  ( v1[0].size() != v2[0].size() ) ){
      std::cout << "行列同士の和算ができません!" << '\n';
      std::cout << "mat[" << v1.size() << "," << v1[0].size() << "] + mat[" << v2.size() << "," << v2[0].size() << "]" << endl;
      exit(1);
   }
   vector<vector<double> > values = v1;
   for( int i = 0; i < values.size(); ++i ){
      for( int j = 0; j < values[0].size(); ++j ){
         values[i][j] = v1[i][j] + v2[i][j];
      }
   }
   return values;
}
//ベクトル+ベクトル
vector<double> calVec::add( vector<double> v1, vector<double> v2 ){
   if( v1.size() != v2.size() ){
      std::cout << "ベクトル同士の和算ができません!" << '\n';
      std::cout << "vec[" << v1.size() << "] + vec[" << v2.size() << "]" << endl;
      exit(1);
   }
   vector<double> values = v1;
   for( int i = 0; i < values.size(); ++i ){
      values[i] = v1[i] + v2[i];
   }
   return values;
}



/* ---  add(2) -----------------------------------------------------------------
ベクトル+スカラーを計算する

input:
vector<double> v                : 計算前のベクトル
( vector<vector<double> > v       : 計算前の行列 )
<double> x                      : 加算するスカラー

output:
vector<double> values           : 計算後のベクトル
( vector<vector<double> > values  : 計算後の行列 ) */
vector<double> calVec::add( vector<double> v, double x ){
   vector<double> values = v;
   for( int i = 0; i < v.size(); ++i ){
      values[i] += x;
   }
   return values;
}
//行列+スカラー
vector<vector<double> > calVec::add( vector<vector<double> > v, double x ){
   vector<vector<double> > values = v;
   for( int i = 0; i < v.size(); ++i ){
      values[i] = add( v[i], x );
   }
   return values;
}



/* --- times(1) ----------------------------------------------------------------
行列*行列を計算する

input:
vector<vector<double> > v1      : 左から掛ける行列
vector<vector<double> > v2      : 右から掛ける行列

output:
<vector<double> > v             : 計算後の行列 */
vector<vector<double> > calVec::times( vector<vector<double> > v1, vector<vector<double> > v2 ){
   vector<vector<double> > v( v1.size(), vector<double>( v2[0].size(), 0 ) );
   if( v1[0].size() != v2.size() ){
      std::cout << "行列の積ができません!" << endl;
      std::cout << "[" << v1.size() << ", " << v1[0].size() << "] * [" << v2.size() << ", " << v2[0].size() << "]" << endl;
      exit(1);
   }
   for( int i = 0; i < v.size(); ++i ){
      for( int j = 0; j < v[0].size(); ++j ){
         for( int k = 0; k < v2.size(); ++k ){
            v[i][j] += ( v1[i][k] * v2[k][j] );
         }
      }
   }
   return v;
}
//行列*ベクトル
vector<double> calVec::times( vector<vector<double> > v1, vector<double> v2 ){
   vector<double> v( v1.size(), 0 );
   if( v1[0].size() != v2.size() ){
      std::cout << "行列の積ができません!" << endl;
      std::cout << "(" << v1.size() << "," << v1[0].size() << ") * (" << v2.size() << ")" << endl;
      exit(1);
   }
   for( int i = 0; i < v.size(); ++i ){
      for( int k = 0; k < v2.size(); ++k ){
         v[i] += ( v1[i][k] * v2[k] );
      }
   }
   return v;
}
//ベクトル*ベクトル
double calVec::times( vector<double> v1, vector<double> v2 ){
   // if(v1.size()!=v2.size()){
   //   std::cout << "ベクトル同士の和算ができません!" << '\n';
   //   exit(1);
   // }
   double values = 0;
   for( int i = 0; i < v1.size(); ++i ){
      values += ( v1[i] + v2[i] );
   }
   return values;
}


/* --- times(2) ----------------------------------------------------------------
行列*スカラーを計算する

input:
vector<vector<double> > v1      : 掛けたい行列
( vector<double> v1               : 掛けたいベクトル )
double x                        : 掛けたいスカラー

output:
<vector<double> > v             : 計算後の行列 */
vector<vector<double> > calVec::times( vector<vector<double> > v1, double x ){
   vector<vector<double> > v( v1.size(), vector<double>( v1[0].size(), 0 ) );
   for( int i = 0; i < v1.size(); ++i ){
      for( int j = 0; j < v1[0].size(); ++j ){
         v[i][j] = v1[i][j] * x;
      }
   }
   return v;
}
//ベクトル*スカラー
vector<double> calVec::times( vector<double> v1, double x ){
   vector<double> v = v1;
   int n = v1.size();
   for( int i = 0; i < n; ++i ){
      v[i] = v1[i] * x;
   }
   return v;
}



/* --- transpose ---------------------------------------------------------------
行列を転置する

input:
vector<vector<double> > v1       : 転置前の行列

output:
vector<vector<double> > transVec : 転置後の行列

*/
vector<vector<double> > calVec::transpose( vector<vector<double> > v1 ){
   vector<vector<double> > transVec( v1[0].size(), vector<double>( v1.size(), 0 ) );
   int n = v1.size(), m = v1[0].size();
   for( int i = 0; i < n; ++i ){
      for( int j = 0; j < m; ++j ){
         transVec[j][i] = v1[i][j];
      }
   }
   return transVec;
}
//ベクトルの転置
vector<vector<double> > calVec::transpose( vector<double> v1 ){
   vector<vector<double> > transVec( 1, vector<double>( v1.size(), 0 ) );
   int n = v1.size();
   for( int i = 0; i < n; ++i ){
      transVec[0][i] = v1[i];
   }
   return transVec;
}


/* --- trace -------------------------------------------------------------------
行列の対角成分を合計

input:
vector<vector<double> > v     :

output:
double values                 : 合計値

*/
double calVec::trace( vector<vector<double> > v ){
   double values = 0;
   int n = v.size();
   if( n != v[0].size() ){
      std::cout << "trace()使えないよ" << endl;
      exit(1);
   }
   for( int i = 0; i < n; ++i ){
      values += v[i][i];
   }
   return values;
}



/* --- cat ---------------------------------------------------------------------
行列同士を連結する
mode = 1 : xの右にyを連結
mode = 2 : xの下にyを連結

*/
vector<vector<double> >  calVec::cat( vector<vector<double> > x, vector<vector<double> > y, int mode ){
   vector<vector <double> > values = x;
   int n = y.size(), m = y[0].size();
   if( mode == 2 ){ //xの右にyを連結
      for( int j = 0; j < n; ++j ){
         for( int k = 0; k < m; ++k ){
            values[j].push_back( y[j][k] );
         }
      }
   }else if( mode == 1 ){ //xの下にyを連結
      for( int i = 0; i < n; ++i ){
         values.push_back( y[i] );
      }
   }

   return values;
}



/* --- sum ---------------------------------------------------------------------
ベクトルの各要素の合計を返す
行列ならば指定した次元の合計を返す

*/
double calVec::sum( vector<double> v ){
   double s = 0, n = v.size();
   for( int i = 0; i < n; ++i ){
      s += v[i];
   }
   return s;
}

vector<double> calVec::sum( vector<vector<double> > v, int dim ){
   int n = v.size();
   int m = v[0].size();
   vector<double> s;
   vector<vector<double> > vtmp = v;
   if( dim == 1 ){
      s = vector<double>( m, 0 );
      vtmp = transpose( vtmp );
      n = m;

   }else if( dim == 2 ){
      s = vector<double>( n, 0 );

   }
   for( int i = 0; i < n; ++i ){
      s[i] = sum( vtmp[i] );
   }

   return s;
}


/* --- sqrtVec -----------------------------------------------------------------
ベクトルの各要素のルートをとったものを返す

*/
vector<double> calVec::sqrtVec( vector<double> v ){
   vector<double> values( v.size(), 0 );
   int n = v.size();
   for( int i = 0; i < n; ++i ){
      values[i] = sqrt( v[i] );
   }

   return values;
}



/* --- expm --------------------------------------------------------------------
行列指数を計算する
固有値と固有ベクトルはEigenを使って計算している

*/
vector<vector<double> > calVec::expm( vector<vector<double> > vec ){
   int n = vec.size();
   int m = vec[0].size();
   vector<double> ssVec( n, 0 );
   vector<vector<double> > ssMat( n, vector<double>( n, 0 ) );
   //行列指数の格納
   vector<vector<double> > values( n, vector<double>( n, 0 ) );
   //固有値の格納
   vector<double> eigenVal( n, 0 );
   //が固有ベクトル(各列が固有ベクトル)
   vector<vector<double> > eigenVec( n, vector<double>( n, 0 ) );

   //行列をEigenで使用できる型に変換
   MatrixXd M( n, m );
   for( int i = 0; i < n; ++i ){
      for( int j = 0; j < m; ++j ){
         M( i, j ) = vec[i][j];
      }
   }
   //固有値，固有ベクトルを計算
   SelfAdjointEigenSolver<MatrixXd> ES;
   ES.compute( M );
   for( int i = 0; i < n; ++i ){
      for( int j = 0; j < n; ++j ){
         if ( i == j ){
            eigenVal[i] = ES.eigenvalues()( i );
         }
         eigenVec[i][j] = ES.eigenvectors()( i, j );
      }
   }

   ssVec = expVec( eigenVal );
   for( int i = 0; i < n; ++i ){
      ssMat[i][i] = ssVec[i];
   }
   values = times( times( eigenVec, ssMat ), inverseVec( eigenVec ) );


   return values;
}




/* --- expVec ------------------------------------------------------------------
行列やベクトルの各要素の指数を返す

*/
vector<double> calVec::expVec( vector<double> v ){
   vector<double> values = v;
   int n = v.size();
   for( int i = 0; i < n; ++i ){
      values[i] = exp( v[i] );
   }
   return values;
}
vector<vector<double> > calVec::expVec( vector<vector<double> > v ){
   vector<vector<double> > values = v;
   int n = v.size();
   for( int i = 0; i < n; ++i ){
      values[i] = expVec( v[i] );
   }
   return values;
}



/* --- inverseVec --------------------------------------------------------------
逆行列だす

*/
vector<vector<double> > calVec::inverseVec( vector<vector<double> > vec ){
   int n = vec.size();
   vector<vector<double> > invVec = eye( n );
   double buf;

   //掃き出し法
   for( int i = 0; i < n; ++i ){
      buf = 1 / vec[i][i];

      for( int j = 0; j < n; ++j ){
         vec[i][j] *= buf;
         invVec[i][j] *= buf;
      }
      for( int j = 0; j < n; ++j ){
         if( i != j ){
            buf = vec[j][i];
            for( int k = 0; k < n; ++k ){
               vec[j][k] -= ( vec[i][k] * buf );
               invVec[j][k] -= ( invVec[i][k] * buf );
            }
         }
      }

   }

   return invVec;
}
