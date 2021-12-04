#ifndef DXNES_H
#define DXNES_H

#include <iostream>

#include <vector>
#include <iomanip>

#include "../eval/evaluator.h"
#include "../calVec/calcurateVector.h"

using namespace std;

class dxnes{
private:
   int g;
   int dim;
   int lam;
   double sigma;
   vector<double> mu;
   // 学習率
   double etaB;
   double etaMu;
   double etaSigma;
   vector<double> etaSigmaSet;
   vector<double> etaBest;
   // 重み
   double alpha;
   vector<double> wiRank;
   vector<double> wiDist;
   vector<double> wRank;
   vector<double> wDist;
   // 個体生成
   vector<vector<double> > z;
   vector<vector<double> > x;
   vector<vector<double> > A;
   vector<vector<double> > B;
   // 評価値
   vector<double> fx;
   // 進化パス
   double base;
   double mueff;
   double csigma;
   vector<double> psigma;
   // 出力
   int eval_gen;
   double fopt;
   vector<double> xopt;
   // 学習継続の判定
   bool condition;

   // 個体生成のためのパラメータ
   double rate;
   double Gsigma;
   vector<double> expZ;
   vector<double> Gdelta;
   vector<double> sumUZ;
   vector<vector<double> > I;
   vector<vector<double> > GM;
   vector<vector<double> > GB;

   // どこでも使っていい，一時退避用の中間変数
   double ssSca;                    // スカラー
   vector<double> ssVec;            // ベクトル
   vector<vector<double> > ssMat;   // マトリクス

   // 学習中にログを表示する true ( or false )
   bool isDispLog;
   bool isDispResult;


public:
   dxnes();
   ~dxnes();

   vector<double> dxnes_start( vector<double> initMu );
   vector<double> dxnes_start( int initDim, int pop_size );
   vector<double> dxnes_start( int initDim, int pop_size, vector<double> initMu, double initSigma );
   vector<double> dxnes_start( int initDim, int pop_size, vector<double> initMu, double initSigma, int maxGens );

   void dispLog();
   void dispResult();

   static void show( vector<double> v );
   static void show( vector<vector<double> > v );

};

#endif // SAMPLE_H
