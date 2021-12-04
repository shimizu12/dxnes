#include "dxnes/dxnes.h"

int main() {

   vector<double> v(5, 10);
   v = calVec::times( calVec::randn( v.size() ), 10000 );


   dxnes dx;
   dx.dxnes_start( v );
   // v = dx.dxnes_start( v );
   // dxnes::show( v );

   return 0;
}
