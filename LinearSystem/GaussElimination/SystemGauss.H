# ifndef __GAUSS_FACTOR_CONTAINER_H__
# define __GAUSS_FACTOR_CONTAINER_H__


# include <iostream>
# include <cmath>
# include <initializer_list>
# include <valarray>
# include <string>
# include <fstream>
# include "Matrix.H" 



namespace mg {
                namespace numeric {       
                                    namespace algebra {


/**-------------------------------------------------------
 *  @brief class container (matrix to be factorized)
 *   using Gauss eliminations and the back-sostitution
 *
 *
 * 
 *
 -------------------------------------------------------*/


template <template<typename> typename container_type, typename type>
class Gauss {

      
    public:

      constexpr Gauss( const std::string& fname ,
                       const std::valarray<type> b ,
                       const std::valarray<type> x)  
                                                      : _a{fname} ,
                                                        _b{b},
                                                        _x{x} 

                              {}
      
      constexpr Gauss( const std::string& fname, const std::string& bname   ) ;

      virtual ~Gauss() = default ;


      constexpr auto getA() const noexcept { return _a ;}            
      constexpr auto getB() const noexcept { return _b ;}            
      constexpr auto getX() const noexcept { return _x ;}            
      
      auto set_b(std::valarray<type> b ) { _b = b ; } 

      auto constexpr print() const noexcept ; 

      auto constexpr gauss() noexcept ; 
      
      auto constexpr printElimination() const noexcept ;
      
      std::valarray<type> constexpr solve() noexcept ; 


    private:

       //
             container_type<type> _a    ;  // generic Matrix container ! 
             std::valarray<type>  _b    ;
             std::valarray<type>  _x    ;
};


//-------------------------------------   IMPLEMENTATION ---------------------------------------------


template <template<typename> typename container_type, typename type>
constexpr Gauss<container_type,type>::Gauss( const std::string& fname,
                                             const std::string& bname ) 
                                                                        : _a{fname} 
                                                                              
{                       
            std::ifstream f(bname, std::ios::in );
            if(!f )
            {
               throw std::runtime_error("Unable to open the known vector in constructor");   
            }

            std::string line ;
            type elem =0 ;
         //--
            _b.resize(_a.size());

            auto i=0;
            while(getline(f,line) )
            {
               std::istringstream ss(line) ;
               ss >> _b[i]  ;
               i++ ;                
                  
            }           
         //--   
            _x.resize(_a.size()) ;
}      
    

template <template<typename> typename container_type, typename type>
auto constexpr Gauss<container_type,type>::gauss() noexcept 
{
      
      auto neq = _a.size( ) ;
      //
      type pivot  ;
      type temp   ;
      type fatt   ;  // moltiplicative factor
      std::size_t  pivotRow ;      // pivot row
      
      for(auto i=1; i <= neq ; i++)
      {
         pivot = _a(i,i) ;    // rersearch of the pivotal element and pivotal row     
         pivotRow = i ;   
         
            // research of a larger pivotal element  
            for(auto j = i+1 ; j <= neq ; j++ )
            {
                if( fabs( _a(i,j) ) > fabs(pivot)  ) 
                {
                     pivot    =  _a( j,i ) ;
                     pivotRow =  j ; 
                }
            }

            // check of larger pivot
            if(pivotRow > i )
            {
               for(auto j=1 ; j <= neq ; j++) 
               {
                             temp = _a(i,j) ;
                          _a(i,j) = _a(pivotRow,j) ;
                   _a(pivotRow,j) = temp ;
               }

               temp           = _b[i-1] ;
               _b[i-1]        = _b[pivotRow-1] ;
               _b[pivotRow-1] = temp ;
            }

            //eliminations of the coefficients from the rows i+1 to neq 
            for(auto j=1+i; j <= neq ; j++)
            {
                  fatt    = _a(j,i)/pivot ;
                  _a(j,i) = 0.0 ;  

                  for(auto k=i+1 ; k < neq ; k++ )
                  {
                      _a(j,k) -= fatt * _a(i,k) ;
                  }
                  _b[j-1] -= fatt* _b[i-1] ;
            }


      
      
      }
}

// solve the problem by mean of back replacement..  

template <template<typename>  typename container_type , typename type >
std::valarray<type> constexpr Gauss<container_type,type>::solve() noexcept          
{
      std::size_t neq = _a.size() ;
             type sost ;
             
      // computing of i-th x
      _x[neq-1] = _b[neq-1] / _a(neq,neq) ;
      //
      for(auto i=neq-1 ; i >= 1 ; i--)
      {
          sost = _b[i-1];
          
          for(auto j=i+1 ; j<= neq ; j++ ){
                  sost -= _a(i,j) * _x[j-1] ;
          }
          _x[i-1] = sost / _a(i,i) ;


      }
      return _x ;
}


//---------------------------- Prinst the whole system solved
//
template <template<typename>  typename container_type , typename type >
auto constexpr Gauss<container_type,type>::print() const noexcept         
{
      for(auto i = 1 ; i <= _a.size() ; i++ ){
            for(auto j = 1 ; j <= _a.size() ; j++ ){
                 std::cout  << std::setw(6) << _a(i,j)  << ' ' ;
            }
            std::cout << "   * " << _x[i-1]  << " =   "   << _b[i-1] << std::endl ; 
      }
}


/** ------------------------------------------------------------------------
        print the system before back replacement (Gauss elimination)

---------------------------------------------------------------------------*/
template <template<typename>  typename container_type , typename type >
auto constexpr Gauss<container_type,type>::printElimination() const noexcept 
{
 for(auto i = 1 ; i <= _a.size() ; i++ ){
            for(auto j = 1 ; j <= _a.size() ; j++ ){
                 std::cout  << std::setw(6) << _a(i,j)  << ' ' ;
            }
            std::cout << "| *  x" + std::to_string(i)  << " =   "   << _b[i-1] << std::endl ; 
  }


}

  }//algebra
 }//numeric
}//mg

# endif
