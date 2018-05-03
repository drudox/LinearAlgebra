# include "Matrix.H"
# include "SystemGauss.H"


using namespace std;

using namespace mg::numeric::algebra ;

int main (){
  

  Gauss<Matrix,double> g1("input17.dat", "inputB1.dat");

  g1.gauss() ;
  valarray<double> x1 = g1.solve() ;
  g1.print() ;

  cout << "----------------------------------------------------------------------------------" << endl;    
  
//  g1.printSolution() ;

   
  for (auto i : x1) 
      cout << i << ' ' ;
  cout << endl ;    


  return 0;
}










