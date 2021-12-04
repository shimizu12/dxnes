#include "dxnes/dxnes.h"

int main() {

   vector<double> v(5, 10);
   // ランダムな初期解を生成
   v = calVec::times( calVec::randn( v.size() ), 10000 );


   dxnes dx;
   v = dx.dxnes_start( v );

   return 0;
}
