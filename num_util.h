#ifndef __NUM_UTIL__
#define __NUM_UTIL__
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <cassert>

namespace br_num_util{


template< typename T> bool is_nan(T& x){
  return x !=x ;
};


class exp_eps_t{
  public:
   exp_eps_t() { 
     ex0 = 1./std::numeric_limits<double>::epsilon();
     x0  = log( ex0 ); 
   };   
   double operator()(double x) const { 
      return  x<x0 ? exp(x) : ex0; 
   };

  private:
    double x0; 
    double ex0;
};




//////////////////////// pretabulated factorials ///////////////////////////
template <typename T> class factorial_t{
  public:
    factorial_t(const size_t upto=5) {
      _fact.resize( 2 );
      _fact[0]=1; _fact[1]=1;
      grow(upto);
    }

   T operator()(const size_t j) {
      if( j>=_fact.size() ){ grow(j); }
      return _fact[j];
   }
  private:
    std::vector<T> _fact;

    bool grow(const size_t upto);
}; 

template<typename T> bool factorial_t<T>::grow(const size_t upto){
  bool ret=true;
  assert(upto>=_fact.size());
  for( size_t j=_fact.size(); j<upto+1; ++j ){
    // sanity check:
    if( _fact[j-1] > std::numeric_limits<T>::max()/static_cast<T>(j)  ){ 
       std::cerr<<"\n\n******** Factorial: overflow @ j="<<j<<"\n\n"; 
       ret=false;
    }
    _fact.push_back( _fact[j-1]*static_cast<T>(j) );
  }
  return ret;
}

////////////////////  a vector_based buffer: VERY CRUDE, handle w/care
template<typename T> class buffer_t{
  public:
    buffer_t(size_t max){
      if( max < 0 || max>_buf.max_size() ){
         std::cerr<<" buffer_t ctor: max="<<max<<"  max_size()="<<_buf.max_size()<<"\n";
         exit(-111);
      }
      _buf.reserve(max);
      _size=0;
      std::cout<<"buffer_t ctor!\n";
    }
    T& operator[](size_t j){
       assert(j>=0 && j<_size);
       return _buf[j];
    }
    void set_size(const size_t new_size) { 
      assert(new_size>=0 && new_size<_buf.max_size());
      _size=new_size;
    }
    size_t size() const { return _size; }
  private:
    std::vector<T> _buf;
    size_t _size;
};


/////////////////////////////////
template<typename vecT> void print_vector(const vecT& vec, std::ostream& os=std::cout){
  os<<"["<<vec.size()<<"]: ";
  for(size_t j=0; j<vec.size(); ++j){ os<<vec[j]<<" "; }
  os<<"\n";
};


inline std::string bool_to_text(const bool b){ return b ? "true" : "false";};

};


#endif
