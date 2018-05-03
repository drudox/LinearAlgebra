# include <iostream>
# include "Matrix.H"
# include <string>



using namespace std;

using namespace mg::numeric::algebra ;

int main(){

  Matrix<double> m1("input17.dat");    
  
 // cout << m1 ;  //(2,2) << endl;    
  m1.print();

  Matrix<double> m2 = m1.minors(2,2)   ;
  
  cout << "-------------------------------------------------" << endl;

  cout << m2 ;    


  cout << "-------------------------------------------------" << endl;


  Matrix<double> m3("input18.dat");    
      
 cout << m3.det() << endl;    
  cout << "-------------------------------------------------" << endl;
      

  return 0;
}


