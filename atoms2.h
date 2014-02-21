#ifndef __HIGHT_MOMENTS__
#define __HIGHT_MOMENTS__
#include <iostream>
#include <cassert>
#include <omp.h>
#include <valarray>
#include <boost/random.hpp>
#include "num_util.h"
#include "hilbert.h"


/****************
TODO:

7. W2 function: decay faster than linear at x<1 : better?
8. "global update": move'em all but the 1st one : helps? 
9. ALPS scheduler & embarassing parallelism : need one?

10. is all the <bool tab> business really needed? CAREFUL w/ Heisenberg-type models + fast-updates maybe
11. move_one-type update: only recalc order N distances --- fast-update a la Rubtsov
12. namespaces 

*****************/


/////////////////////////////////////////// coordinates
struct coords_t{
  double x;
  double y;
};

inline std::ostream& operator<<(std::ostream& stream, const coords_t& c){
  stream<<"("<<c.x<<","<<c.y<<") ";
  return stream;
};


///////////////////////////////////////////// atoms
typedef std::valarray<coords_t> atoms_t;

inline double atoms_distance(const coords_t& one, const coords_t& two){
  double dx = one.x-two.x;
  double dy = one.y-two.y;
  return sqrt( dx*dx + dy*dy );
};


inline double comp_dist3(const coords_t& one, const coords_t& two){
  double ttt = atoms_distance(one,two);
  return ttt*ttt*ttt;
};



///////////////////////////////////// proxies & subsets
typedef std::valarray<size_t> proxy_t;
typedef std::vector< size_t > vect_index_t;
typedef std::vector< std::pair< double, vect_index_t > > subsets_t;


/////////////////// base functionality for R-s //////////////
class baseR{
  public:
    baseR(const double offset=1., const bool verbose=1) : _offset(offset), 
                                                          _verbose(verbose),
                                                          _mexp() {}
    double magn(const spins_t& spins) const;

    double offset() const {return _offset;}
    double mexp(const double x) const {return _mexp(x);} 
  protected:
    double _offset;
    bool _verbose;
    br_num_util::exp_eps_t _mexp;
};


inline double baseR::magn(const spins_t& spins) const{
  return 1.*std::accumulate( spins.begin(), spins.end(),0 );
};




///////// objects of this one are not supposed to exist on their own
template<bool tab> class rtab_stor{
  public:
    rtab_stor(const size_t Nat) {}
}; // yes, empty!


template<>
class rtab_stor<true>{
   public:
     rtab_stor(const size_t Nat) {
       tab_dist3_.resize( Nat*(Nat-1)/2 );
     } 
     void tabulate_dist3(const atoms_t& a) const;   // will receive atoms_
  //protected:
     mutable std::valarray<double> tab_dist3_;
};


inline void rtab_stor<true>::tabulate_dist3(const atoms_t& a) const{
   size_t n=a.size();
   size_t cnt=0;
   for( size_t i=0; i<n; ++i){
     for( size_t j=i+1; j<n; ++j ){
       cnt=i*n-i*(i+3)/2 +j -1; // 1) contigious,  2) spans [0,n(n-1)/2)
       tab_dist3_[cnt]=comp_dist3( a[i], a[j] );
     }
   }
};




/*** super-tabulated versions of R & F *************/
template<bool tab> class safeR : public baseR  /*@#$ , public rtab_stor<tab> */{
  public:
     safeR(const vect_index_t& subset, const double offset=1, const bool verbose=false) :  
                                                   baseR(offset, verbose),
                                                   /* @#$ rtab_stor<tab>( subset.size() ){  */
                                                   rtab_( subset.size() ) {
       basis_=tab_basis( subset.size());
       proxy_.resize(subset.size());
       atoms_.resize(subset.size());
       for(size_t j=0; j<proxy_.size(); ++j){ proxy_[j] = subset[j]; }

       enes_.resize(basis_.size());
       magns_.resize(basis_.size());
     }
     proxy_t proxy() const {return proxy_; }

     double operator()(const atoms_t& atoms) const;  // need specializing
     double energy(const spins_t& spins) const ;
     double dist3(const size_t& i, const size_t& j) const; 
 

   private:
     basis_t basis_;
     proxy_t proxy_; 
     mutable atoms_t atoms_;
     rtab_stor<tab> rtab_;
     double R_impl() const;

     mutable std::vector<double> enes_;   // working arrays 
     mutable std::vector<double> magns_;
};


/// now specialize dist3:
template<>  // a non-tab version
inline double safeR<false>::dist3(const size_t& i, const size_t& j) const{ 
  return comp_dist3( atoms_[i], atoms_[j] ); 
};

/// a tab version
template<>
inline double safeR<true>::dist3(const size_t& i, const size_t& j) const{
  size_t cnt= i*proxy_.size()-i*(i+3)/2 +j -1; // cf tabulate_dist3()
  return rtab_.tab_dist3_[cnt];
};


// operator()
template<>
inline double safeR<false>::operator()(const atoms_t& atoms) const{
   atoms_=atoms[proxy_];
   return R_impl();   // operates on atoms_
};

template<>
inline double safeR<true>::operator()(const atoms_t& atoms) const{
   atoms_=atoms[proxy_];
   rtab_.tabulate_dist3(atoms_); 
   return R_impl();   // operates on atoms_
};




/// energy & R_impl are template parameter independent really (only propagate the templ param)
template<bool tab> double safeR<tab>::energy(const spins_t& spins) const{
    assert( proxy_.size()==spins.size());
    double ene = 0.;
    for( size_t i=0; i<proxy_.size(); ++i ){
      for( size_t j=i+1 ; j<proxy_.size(); ++j ){
        ene += (spins[i]*spins[j]-_offset) / dist3(i,j); 
      }
    }
    return ene;  
};


template<bool tab>
double safeR<tab>::R_impl() const{
   assert( atoms_.size()>=proxy_.size() );
   size_t spin_space_dim=basis_.size();

   //double ene_max = std::numeric_limits<double>::min(); 
   for( size_t idx=0; idx<spin_space_dim; ++idx ){
     magns_[idx] = magn( basis_[idx] );
     enes_[idx] = energy( basis_[idx] );  // uses member atoms_
   }

   double m2sum=0., Zsum=0.;
   for( size_t idx=0; idx<spin_space_dim; ++idx ){
     double m = magns_[idx];
     double b = mexp(enes_[idx]);
     m2sum += m*m*b;
     Zsum += b;
   }

   return m2sum/Zsum; 
};




/////////////////// F: tabulates i) the subsets, ii) factorials @ ctor, 3) holds its R-s
template<bool tab> class objF3{
  public:
     objF3(const size_t Nat, const double offset=1, const bool verbose=false) : _Nat(Nat), _fact(Nat){ 
       subsets_t S=tabulate_subsets(Nat);
      
////// uncomment these pragmas and it segfaults [gcc @ymir: 2nd only is OK; icpc @elysium: either fails]
       //#pragma omp parallel for 
       for( size_t j=0; j<S.size(); ++j ){
         coefs_.push_back( S[j].first );
       }

      // #pragma omp parallel for 
       for(size_t j=0; j<S.size(); ++j){
         safeR<tab> thisR( S[j].second, offset, verbose ); 
         Rs_.push_back( thisR );         
      }
      
      offset_=offset;
     }
     std::pair<bool, double> operator()(const atoms_t&, const bool verbose=false) const ;
     double fact(const size_t j) const { return _fact(j); }
     double local_Nat() const { return _Nat; }
     double offset() const {return offset_; }
  private:
    size_t _Nat;  
    double offset_;
    mutable br_num_util::factorial_t<double> _fact;
    std::vector< safeR<tab> > Rs_; 
    std::vector< double> coefs_;

    subsets_t tabulate_subsets(const size_t NNN){
      subsets_t outp;
      for(size_t j=0; j<NNN-1 ; ++j){
        vect_index_t idx;
        idx.push_back(j);
        //std::cout<<j<<"  ";
        generate_subsets(idx,NNN, outp);
      }
      return outp;
    };

    void generate_subsets(vect_index_t idx, const size_t N, subsets_t& outp){
      if( idx.size()<N ){
        for( size_t j=idx[idx.size()-1]+1; j<N ; ++j ){
          vect_index_t a1(idx);
          a1.push_back(j);
          double sign=pow(-1.,1.*(N-a1.size()));
          outp.push_back( std::make_pair(sign, a1) );
          //std::cout<<"\n-->"<<j;
          //print_vector(a1);
          generate_subsets(a1,N, outp);
        }
      }
    };
};


template<bool tab>
std::pair<bool, double> objF3<tab>::operator()(const atoms_t& atoms, const bool verbose) const{
    assert(atoms.size()==_Nat);
    assert(atoms.size()>2);
    double sum=0.;
    #pragma omp parallel for reduction(+:sum) 
    for( size_t i=0; i<coefs_.size(); ++i ){
      sum += coefs_[i] * Rs_[i](atoms) ;
      if(verbose){
        std::cout.precision(18);
        std::cerr<<i<<"  sumsofar="<<sum<<"  coef="<<coefs_[i]<<" R="<<Rs_[i](atoms)<<"  "; //<<"\n";
        br_num_util::print_vector(Rs_[i].proxy(), std::cerr);
      }
    }
    int sign = atoms.size()%2==0 ? -1 : +1; 
    sum += 1.*sign*atoms.size(); 

    if(verbose){
      std::cout<<"@F: sum="<<sum<<"  lim="<<std::numeric_limits<double>::epsilon()*10*coefs_.size()<<"\n";
    }

    double lost_it=false;
    if( sum < std::numeric_limits<double>::epsilon()*10*coefs_.size() ){ lost_it=true; sum=0;  }

    if(verbose){
      std::cout<<"@F: sum="<<sum<<"  lost_it="<< (lost_it ? "true" : "false") <<"\n";
    }

    return std::make_pair(lost_it, sum/ _fact(atoms.size()) );
}; 


#endif


