#include "dxnes.h"

using namespace std;

/*
コンストラクタ
パラメータ初期化
*/
dxnes::dxnes(){
   // 世代数
   g = 0;

   // 学習継続の判定
   condition = true;

   // 学習中のログを表示する true( or false )
   isDispLog = false;
   isDispResult = true;
}



/*
デストラクタ
メモリの解放
*/
dxnes::~dxnes(){
   vector<double>().swap( mu );
   vector<double>().swap( expZ );
   vector<double>().swap( xopt );
   vector<double>().swap( sumUZ );
   vector<double>().swap( ssVec );
   vector<double>().swap( wRank );
   vector<double>().swap( wDist );
   vector<double>().swap( Gdelta );
   vector<double>().swap( psigma );
   vector<double>().swap( wiRank );
   vector<double>().swap( wiDist );
   vector<double>().swap( etaBest );
   vector<double>().swap( etaSigmaSet );
   vector<vector<double> >().swap( A );
   vector<vector<double> >().swap( B );
   vector<vector<double> >().swap( I );
   vector<vector<double> >().swap( x );
   vector<vector<double> >().swap( z );
   vector<vector<double> >().swap( GB );
   vector<vector<double> >().swap( GM );
   vector<vector<double> >().swap( ssMat );

}


vector<double> dxnes::dxnes_start( vector<double> initMu ){
   return dxnes::dxnes_start( initMu.size(), 100, initMu, 1.0, 3000 );
}

vector<double> dxnes::dxnes_start( int initDim, int pop_size ){
   return dxnes::dxnes_start( initDim, pop_size, vector<double>( 1, 10 ), 1.0, 3000 );
}

vector<double> dxnes::dxnes_start( int initDim, int pop_size, vector<double> initMu, double initSigma ){
   return dxnes::dxnes_start( initDim, pop_size, initMu, initSigma, 3000 );
}

vector<double> dxnes::dxnes_start( int initDim, int pop_size, vector<double> initMu, double initSigma, int maxGens ){

   // --- Algorithm parameters -------------------------------------------------
   // lambda : 個体数
   dim = initDim;
   lam = (pop_size < 2)? 2: 2*floor( pop_size / 2 );

   // learning rates : 学習率
   etaMu = 1.0;
   ssSca = 1.0 + lam / ( lam + ( 2 * dim ) );
   etaSigmaSet = vector<double>( 3, 0 );
   etaSigmaSet[0] = 3.0;
   etaSigmaSet[1] = 0.5 * ssSca;
   etaSigmaSet[2] = ssSca;

   etaBest = vector<double>( 2, 0 );
   ssSca = lam + ( 2 * dim * dim ) + 100;
   etaBest[1] = lam / ssSca;
   ssSca = ( lam + 2 * dim ) / ssSca;
   etaBest[0] = ssSca;
   ssSca = min( 1.0, sqrt( lam / dim ) );
   etaBest[0] = etaBest[0] * ssSca;

   // utility function : 重み
   wiRank = vector<double>( lam, 0 );
   wRank = vector<double>( lam, 0 );
   ssSca = log( ( lam / 2 ) + 1 );
   for( int i = 0; i < lam; ++i ){
      wiRank[i] = max( 0.0, ssSca - log( i + 1 ) );
   }

   ssSca = 1 / calVec::sum( wiRank );
   wRank = calVec::add( calVec::times( wiRank, ssSca ), - ( 1 / lam ) );
   alpha = 0.9 + 0.15 * log( dim );

   // evolution path : 進化パス
   ssSca = 1 - 1 / ( 4 * dim );
   ssSca = ssSca + 1 / ( 21 * dim * dim );
   base = sqrt( dim ) * ssSca;

   ssVec = calVec::add( wRank, 1 / lam );
   for( int i = 0, size = ssVec.size(); i < size; ++i ){
      ssVec[i] = pow( ssVec[i], 2 );
   }
   mueff = 1 / calVec::sum( ssVec );

   ssSca = ( mueff + 2.0 ) / ( dim + mueff + 5.0 );
   csigma = ssSca / sqrt( dim );

   // --- Algorithm start ------------------------------------------------------
   // initialize
   if( initMu.size() != dim ){
      vector<double> tmp( dim, 1 );
      for( int i = 0; i < dim; ++i ){
         tmp[i] = initMu[0];
      }
      initMu = tmp;
   }
   mu = initMu;
   sigma = initSigma;
   I = calVec::eye( dim );
   B = I;
   wiDist = vector<double>( lam, 0 );
   wDist = vector<double>( lam, 0 );
   psigma = vector<double>( dim, 0 );
   xopt = mu;
   fopt = eval::evalFunc( mu );

   cout << "初期個体 ";
   show( initMu );

   // ----- mainloop -----------------------------------------------------------
   while( condition ){
      g ++;
      A = calVec::times( B, sigma );

      // 正規個体z， 評価個体x の生成
      z = calVec::randn( dim, ( lam / 2 ) );
      z = calVec::cat( z, calVec::times( z, - 1 ), 2 );
      x = calVec::add( calVec::times( A, z ), calVec::repmat( mu, 1, lam ) );

      // 評価個体の評価
      fx = eval::evalFunc( x );
      // 評価値により個体をソート
      x = calVec::transpose( x );
      z = calVec::transpose( z );
      vector<double> tmpX, tmpZ;
      for( int i = 0, size = fx.size(); i < size; ++i ){
         for( int j = i; j < size; ++j ){
            if( fx[i] > fx[j] ){
               tmpX = x[i]; x[i] = x[j]; x[j] = tmpX;
               tmpZ = z[i]; z[i] = z[j]; z[j] = tmpZ;
               ssSca = fx[i]; fx[i] = fx[j]; fx[j] = ssSca;
            }
         }
      }


      // もっとも評価の高い個体x[0]とその評価値fx[0]の更新
      if( fx[0] <= fopt ){
         eval_gen = g;
         xopt = x[0];
         fopt = fx[0];
      }
      x = calVec::transpose( x );
      z = calVec::transpose( z );

      // 進化パスpsigmaの更新
      ssMat = calVec::repmat( calVec::transpose( wRank ), dim, 1 );
      for( int i = 0, size_i = ssMat.size(); i < size_i; ++i ){
         for( int j = 0, size_j = ssMat[0].size(); j < size_j; ++j ){
            ssMat[i][j] *= z[i][j];
         }
      }
      sumUZ = calVec::sum( ssMat, 2 );
      ssSca = csigma * ( 2 - csigma ) * mueff;
      ssVec = calVec::times( sumUZ, sqrt( ssSca ) );
      psigma = calVec::add( calVec::times( psigma, ( 1- csigma ) ), ssVec );
      rate = calVec::norm( psigma ) / base;

      // 重みwと学習率etaの更新
      // 移動期
      GM = vector<vector<double> >( dim, vector<double>( dim, 0 ) );
      if( 1.0 <= rate ){
         ssMat = z;
         for( int i = 0, size_i = ssMat.size(); i < size_i; ++i ){
            for( int j = 0, size_j = ssMat[0].size(); j < size_j; ++j ){
               ssMat[i][j] *= z[i][j];
            }
         }
         ssVec = calVec::sqrtVec( calVec::sum( ssMat, 1 ) );
         expZ = calVec::expVec( calVec::times( ssVec, alpha ) );
         for( int i = 0; i < lam; ++i ){
            wiDist[i] = expZ[i] * wiRank[i];
         }
         ssSca = 1 / calVec::sum( wiDist );
         wDist = calVec::add( calVec::times( wiDist, ssSca ), - ( 1 / lam ) );
         ssMat = calVec::repmat( calVec::transpose( wDist ), dim, 1 );
         for( int i = 0, size_i = ssMat.size(); i < size_i; ++i ){
            for( int j = 0, size_j = ssMat[0].size(); j < size_j; ++j ){
               ssMat[i][j] *= z[i][j];
            }
         }
         Gdelta = calVec::sum( ssMat, 2 );
         ssSca = calVec::sum( wDist );
         ssMat = calVec::times( ssMat, calVec::transpose( z ) );
         GM = calVec::add( ssMat, calVec::times( calVec::times( I, ssSca ), - 1 ) );
         etaSigma = etaSigmaSet[0];
         etaB = etaBest[0];

      // 停滞期
      }else if( rate < 1.0 ){
         Gdelta = sumUZ;
         ssMat = calVec::repmat( calVec::transpose( wRank ), dim, 1 );
         for( int i = 0, size_i = ssMat.size(); i < size_i; ++i ){
            for( int j = 0, size_j = ssMat[0].size(); j < size_j; ++j ){
               ssMat[i][j] *= z[i][j];
            }
         }
         ssSca = calVec::sum( wRank );
         ssMat = calVec::times( ssMat, calVec::transpose( z ) );
         GM = calVec::add( ssMat, calVec::times( calVec::times( I, ssSca ), - 1 ) );
         etaSigma = etaSigmaSet[1];
         etaB = etaBest[1];

         // 収束期
         if( rate < 0.1 ){
            etaSigma = etaSigmaSet[2];
         }

      }

      // --- update the parameters ---------------------------------------------
      // 平均ベクトルmu，共分散行列Bの更新
      ssMat = calVec::times( A, etaMu );
      mu = calVec::add( mu, calVec::times( ssMat, Gdelta ) );
      Gsigma = calVec::trace( GM ) / dim;
      GB = calVec::add( GM, calVec::times( calVec::times( I, Gsigma ), - 1 ) );
      ssSca = ( etaSigma / 2 ) * Gsigma;
      sigma = sigma * exp( ssSca );
      ssMat = calVec::times( GB, ( etaB / 2 ) );
      B = calVec::times( B, calVec::expm( ssMat ) );


      // 学習終了条件
      if(   ( sigma < 10e-15 )
         || ( maxGens <= g )
         || ( fx[0] < 10e-10 )
      ){
         condition = false;
      }

      // 学習ログの表示
      if( isDispLog ){ dispLog(); }

   } // --- mainloop終了 --------------------------------------------------------

   // リザルト表示
   if( isDispResult ){ dispResult(); }

   return xopt;
}


/* dispLog ---------------------------------------------------------------------
学習中のログを標準出力

*/
void dxnes::dispLog(){
   cout << endl << "------------------------------------------------------------------------------" << endl;
   if( rate < 0.1 ){
      cout << "---- 収束期";
   }else if( rate < 1.0 ){
      cout << "---- 停滞期";
   }else{
      cout << "---- 移動機";
   }
   cout << " " << right << setw( 5 ) << g << " 回目の世界線 " << "----------------------------------------------- " << endl;

   // cout << left << setw( 2 ) << "wRank ";
   // show( wRank );
   // cout << left << setw( 2 ) << "wDist ";
   // show( wDist );
   // cout << left << setw( 2 ) << "A ";
   // show( A );
   // cout << left << setw( 2 ) << "z" << " = ";
   // show( z );
   // cout << left << setw( 2 ) << "x ";
   // show( x );
   // cout << left << setw( 2 ) << "x[0] ";
   // show( transpose( x )[0] );
   // cout << "-------------------------------- " << endl;
   cout << left << setw( 8 ) << "fx[0]" << " = " << fx[0] << endl;
   // cout << "-------------------------------- " << endl;
   // cout << left << setw( 2 ) << "fx" << " = ";
   // show( fx );
   // cout << left << setw( 2 ) << "xopt ";
   // show( xopt );
   // cout << left << setw( 8 ) << "fopt" << " = " << fopt << endl;
   // cout << left << setw( 4 ) << "sumUZ ";
   // show( sumUZ );
   // cout << left << setw( 4 ) << "psigma ";
   // show( psigma );
   // cout << left << setw( 8 ) << "rate" << " = " << rate << endl;
   // cout << left << setw( 2 ) << "GM ";
   // show( GM );
   // cout << left << setw( 4 ) << "Gdelta ";
   // show( Gdelta );
   // cout << left << setw( 4 ) << "mu ";
   // show( mu );
   // cout << left << setw( 8 ) << "fmu" << " = " << evalfanc( mu ) << endl;
   // cout << left << setw( 2 ) << "GB ";
   // show( GB );
   // cout << left << setw( 8 ) << "Gsigma" << " = " << Gsigma << endl;
   // cout << left << setw( 8 ) << "sigma" << " = " << sigma << endl;
   // cout << left << setw( 2 ) << "B ";
   // show( B );

}

/* dispResult ------------------------------------------------------------------
学習結果を標準出力

*/
void dxnes::dispResult(){
   cout << endl << "==============================================================================" << endl;
   cout << "==== result ==================================================================" << endl;

   // 何世代まで探索したか
   cout << "最終世代数 : " << g << endl;
   // 探索した最後の世代 * 一世代当たりの評価回数
   cout << "総評価回数 : " << g * lam << " ( = " << g << " × " << lam << " )" << endl;
   // 最小評価値のある世代
   cout << endl << "最高評価世代 : " << eval_gen << endl;
   // 最小評価値
   cout << "最小評価値 : " << fopt << endl;
   // 最優個体
   cout << endl << "最優個体 ";
   show( xopt );

}


/* show ------------------------------------------------------------------------
行列やベクトルの各要素を表示する

*/
void dxnes::show( vector<vector<double> > v ){
   int var = v.size(), hol = v[0].size();
   cout << "mat[" << var << ", " << hol << "]" << endl;
   for( int i = 0; i < var; ++i ){
      cout << right << setw( 4 ) << i << " |";
      for( int j = 0; j < hol; ++j ){
         cout << setw( 14 ) << v[i][j] << " ";
      }
      cout << endl;
   }
}
void dxnes::show( vector<double> v ){
   int var = v.size();
   cout << "vec[" << var << "]" << endl;
   cout << right << setw( 4 ) << " = [ ";
   for( int i = 0; i < var; ++i ){
      cout << setw( 14 ) << v[i] << " ";
   }
   cout << "]" << endl;
}
