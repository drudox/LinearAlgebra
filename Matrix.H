#ifndef ___MATRIX_H
#define ___MATRIX_H


# include <initializer_list>
# include <valarray>
# include <cmath>
# include <iomanip>
# include <string>
# include <fstream>
# include <exception>
# include <sstream>
# include <cassert>
# include <iostream>
# include <vector>



namespace mg {
               namespace numeric {
                                    namespace algebra {


// forward declaration
template<typename T>
class Matrix ;

// external function
template <typename T>
std::ostream& operator<<(std::ostream& os , const Matrix<T>& m) noexcept ; 



//----------------------------------------------------------------------------
//
template <typename Type = double >
class Matrix {

  public:
    
      template <typename T>
      friend std::ostream& operator<<(std::ostream& os , const Matrix<T>& m) noexcept ; 

    public:     
    
     constexpr Matrix() {}
      
     constexpr Matrix(const std::size_t r , const std::size_t c ) ;

     constexpr Matrix( std::initializer_list<std::initializer_list<Type>> row ) ;
      
     constexpr Matrix( const std::string& fname  ) ;  

     constexpr Matrix<Type> minors(const std::size_t r , const std::size_t c ) const ;

     Type& operator()(const std::size_t i, const std::size_t j) noexcept ;

     constexpr const Type& operator()(const std::size_t i, const std::size_t j) const noexcept ;

     constexpr auto print() const noexcept ; 
      
     constexpr auto size() const noexcept { return _rows ;} 
      
     constexpr auto size1() const noexcept { return _rows ;} 
      
     constexpr auto size2() const noexcept { return _columns ;} 
      
     constexpr auto cofactor( Matrix<Type>& temp , std::size_t p , std::size_t q ) noexcept; 
      
     constexpr auto det() -> Type  ;
     
   
   
   private: 
      
      std::valarray<Type>    _data    ;
      std::size_t            _rows    ;
      std::size_t            _columns ;  

};


//----------------------------------------------------------------------------
//
template <typename T>   
constexpr Matrix<T>::Matrix(const std::string& fname ) 
{
      std::ifstream f(fname , std::ios::in);

      if(!f)
      {
         throw std::runtime_error("Error in constructor unable to oper file for reading");   
      }
      
      if( fname.find(".mtx") != std::string::npos )
      {

         std::string line ;
         T elem = 0 ;
         T dummy = 0 ;
         auto r = 0 , c = 0 ;       
         auto i=0 , j=0 ;
         getline(f,line) ; // jump out the header
            
         while(getline(f,line))
         {
            std::istringstream ss(line) ;  
            if(i==0)
            {
               ss >> _rows ; 
               ss >> _columns ;
               ss >> dummy ;
             
             _data.resize(_rows*_columns);
          }
          else
          {
             ss >> r ; 
             ss >> c ;
             ss >> elem ;
             
             this->operator()(r,c) = elem ;
          }      
         }   
      }
      else
      {
         T elem = 0;   
         auto i=0, j=0 , k=0;
         
         std::string tmp, line ;
         
         while(getline(f,line))
         {  
            //j=0 ;
            std::istringstream ss(line);
            while(ss >> elem)
            {
                  if(i==0) j++;
            }
            i++ ;
         }
         this->_rows    = i ;
         this->_columns = j ;
         
       //  std::cout << i << 'x' << j << std::endl;
            
         _data.resize(_rows * _columns );   
         
         f.clear() ;
         f.seekg(0 , std::ios::beg);
         
         k=0;
         while(getline(f,line))
         {  
            std::istringstream ss(line);
            while(ss >> elem)
            {     
                  _data[k] = elem ;
                  k++ ;
            }
            i++ ;
         }

      }
}


//-- 
//

template <typename T>
constexpr Matrix<T>::Matrix(const std::size_t r , const std::size_t c ) 
{
      _rows = r ;
      _columns = c ;
      _data.resize(r*c);
}

//---
//
template <typename T>
T& Matrix<T>::operator()(const std::size_t i, const std::size_t j) noexcept 
{
   return _data[(i-1) * _columns + (j-1)] ;
}
//---
//
template <typename T>
constexpr const T& Matrix<T>::operator()(const std::size_t i, const std::size_t j) const noexcept 
{
   return _data[(i-1) * _columns + (j-1)] ;
}



template <typename T>
constexpr auto Matrix<T>::print() const noexcept 
{
    for(auto i=1; i <= _rows ; i++){
       for(auto j=1; j <= _columns ; j++){
            std::cout << std::setw(6) << this->operator()(i,j) << ' ' ;
       }
       std::cout << std::endl;
    }


}


// Non member function
//

template <typename T>
std::ostream& operator<<(std::ostream& os , const Matrix<T>& m) noexcept 
{
      for(auto i=1 ; i<= m._rows ; i++){
           for(auto j=1; j<= m._columns ; j++){
                  os << std::setw(6) << m(i,j) << ' ' ;
           }
           os<< std::endl; 
      }

}


template <typename T>
constexpr Matrix<T> Matrix<T>::minors(const std::size_t r , const std::size_t c ) const 
{
      //if(r > 0 && r <= row && c > 0 && c <= columns)
            
      if(r > 0 && r <= _rows && c > 0 && c <= _columns ) // 1-based
      {
          Matrix<T> res{_rows-1 , _columns-1} ;

          for(std::size_t ri = 1 ; ri <= ( _rows - (r >= _rows)); ri++ )
          {
             for(std::size_t ci = 1 ; ci <= (_columns - (c >= _columns)); ci++ )
             {
                  res(ri - (ri > r), ci - (ci > c)) = this->_data[(ri-1) * _columns + (ci -1)];
             }
          }  
          return res;  
      }
      else
      {
          throw std::runtime_error("Index for minor out of range");  
      }
       
}


template <typename T>
constexpr auto Matrix<T>::det() -> T  
{
   if( !(_rows == _columns ))
   {
     throw std::runtime_error("Matrix must be SQUARE for compute the DETERMINANT WITH MINORS-method");     
   }
   
   Matrix<T> a{*this};   // matrix used to compute the step 

   T d = 0 ; // value of determinant   
   auto rows = a.size1();
   auto cols = a.size2();    
   
   if( rows == cols) // a is a square Matrix 
   {
      if(rows==1)  // 1x1 matrix
      {
         d = a(1,1);   
      }
      else if(rows==2)
      {
         d = a(1,1)*a(2,2) - a(2,1)*a(1,2) ;    
      }
      else
      {
            // matrix equal or larger than 3x3
            for(auto c=1 ; c <= cols ; c++)
            {
                Matrix<T> M = a.minors(1,c);

                d += pow(-1,1+c) * a(1,c) * M.det();
            }
      
      }
    } 
    return d;
}



//function to get cofactor in a temporary matrix




  }//algebra
 }//numeric
}//mg
# endif
