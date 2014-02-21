#ifndef __HILBERT_SPACE_HELPER__
#define __HILBERT_SPACE_HELPER__
#include <iostream>
#include <valarray>
#include <algorithm>

/* 
given N *ISING* spins the Hilb space has the dimension 2**N, and
is spanned by the direct product of the spins. 

Each spin s={+1,-1}
Indexing is idx = S_0 + S_1*2 + ... + S_j*2^j + ... S_N*2^{N-1}, 
where S(s=+1)=1, S(s=-1)=0 

*/

/////////////////// spins
typedef short spin_t;
typedef std::vector<spin_t> spins_t;

typedef std::vector<spins_t> basis_t; 

class idx2vec_t{
  public:
    idx2vec_t( const size_t N );
    size_t idx_max() const  {return _max; };
    spins_t idx2vec( const size_t idx ) const ; 
  private:
    std::vector<size_t> _f;
    size_t _max;
    size_t _N;
};



inline idx2vec_t::idx2vec_t( const size_t N ){
  size_t ttt=1;
  _f.resize(N);
  for( size_t j=0; j<N ; ++j ){
    _f[j]=( ttt );
    ttt *=2;
  }
  std::reverse(&_f[0], &_f[N] ) ; //   _f.begin(), _f.end());
  _max=_f[0] *2;
  _N=N;
}



inline spins_t idx2vec_t::idx2vec(const size_t idx) const {
  spins_t vec;
  vec.resize(_N);
  fill(vec.begin(), vec.end(), -101);
  size_t iii=idx;
  for( size_t j=0; j<_N; ++j ){
    vec[j]=iii/_f[j];
    iii -= vec[j]*_f[j];
  }
  // {0,1}  -> {-1,+1}
  for( size_t j=0; j<vec.size(); ++j){ 
   if(vec[j]==0){vec[j]=-1 ;}
  }
  return vec;
}


inline basis_t tab_basis(const size_t N){
  basis_t basis;
  idx2vec_t idx2vec(N);
  for(size_t j=0; j<idx2vec.idx_max(); ++j){
    basis.push_back( idx2vec.idx2vec(j) );
  }
  return basis; 
}



#endif

