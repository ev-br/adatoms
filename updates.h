#ifndef __UPDATES__
#define __UPDATES__
#include "addr_accpt.h"
#include "num_util.h"
#include "atoms2.h"
#include "matr_ops2.h"

/////////////////////////// power-law distribution: W(x) = \alpha / (1+x)^{1+\alpha}

// can't partially specialize members, hence use the trick from 
// http://stackoverflow.com/questions/1501357/template-specialization-of-particular-members
template<int a> struct i2t {};

template <int a, class RNG> class pwr_weight_t{
  public:
    pwr_weight_t(RNG& rndm) : _rndm(rndm) {}
    int alpha() const {return a;}

    double play_new() const { return play_new_impl( i2t<a>() );  /* return 1./pow(_rndm(), 1./a)-1; */ }
    double probab( const double x) const { return probab_impl( i2t<a>(), x );  /*return 1.*a/pow(1.+x, 1+a);*/ }

  private:
    RNG& _rndm;

    // pow() is excruciatingly slow, hence specialize for \propto(1+x)**3
    double play_new_impl(i2t<2>) const {
      return 1./sqrt(_rndm())-1.; 
    };

    double probab_impl(i2t<2>, double x){
      double xxx=1.+x; 
      xxx = xxx*xxx*xxx;
      return 2./xxx;  
    };

};



// move one atom, uniform distribution
template <typename RNG> addr_accpt move_one_uniform( RNG& rndm, atoms_t& atoms, det_weight_t& www, double dx=12){
  addr_accpt thisone;
  thisone.address();

  atoms_t atoms_new=atoms; // NB: a temporary
  size_t which=(size_t)(1+(atoms.size()-1)*rndm()); 
  assert( which>0 && which<atoms.size() );

  atoms_new[which].x += dx*(rndm() -0.5);
  atoms_new[which].y += dx*(rndm() -0.5);

  int info; 
  double ratio ;
  boost::tie(info, ratio)= www.Wrat(atoms_new, atoms);

  if(info!=0){ 
     std::cerr<<"Houston, we have a problem; info="<<info<<" w/ cnf being...";
     br_num_util::print_vector(atoms, std::cerr);
     exit(-111);
  }

  //metropolis
  if( ratio >1 || rndm() <ratio ){ atoms=atoms_new;  thisone.accept(); }

  return thisone;

};


// move one atom, P(dR) = (1/2\pi)*1/(1+dR)**3
template <typename RNG> addr_accpt move_one_r3( RNG& rndm, atoms_t& atoms, det_weight_t& www){
  addr_accpt thisone;
  thisone.address();

  atoms_t atoms_new=atoms; // NB: a temporary
  size_t which=(size_t)(1+(atoms.size()-1)*rndm()); 
  assert( which>0 && which<atoms.size() );


  double phi=2.*M_PI*rndm();
  double dR=1./sqrt(rndm()) -1.; // p(dR)\propto 1/(1+dR)**3
  atoms_new[which].x += dR*cos(phi);
  atoms_new[which].y += dR*sin(phi);

  int info; 
  double ratio ;
  boost::tie(info, ratio)= www.Wrat(atoms_new, atoms);

  if(info!=0){ 
     std::cerr<<"Houston, we have a problem; info="<<info<<" w/ cnf being...";
     br_num_util::print_vector(atoms,std::cerr);
     exit(-111);
  }

  //metropolis
  if( ratio >1 || rndm() <ratio ){ atoms=atoms_new;  thisone.accept(); }

  return thisone;
};


// move all atoms, uniform distribution
template <typename RNG> addr_accpt move_all_uniform( RNG& rndm, atoms_t& atoms, det_weight_t& www, double dx=12){
  addr_accpt thisone;
  thisone.address();

  atoms_t atoms_new=atoms; // NB: a temporary
  for( size_t j=1; j<atoms.size(); ++j ){
    atoms_new[j].x += dx*(rndm() -0.5);
    atoms_new[j].y += dx*(rndm() -0.5);
  }
 
  int info; 
  double ratio ;
  boost::tie(info, ratio)= www.Wrat(atoms_new, atoms);

  if(info!=0){ 
     std::cerr<<"Houston, we have a problem; info="<<info<<" w/ cnf being...";
     br_num_util::print_vector(atoms, std::cerr);
     exit(-111);
  }

  //metropolis
  if( ratio >1 || rndm() <ratio ){ atoms=atoms_new;  thisone.accept(); }

  return thisone;

};





#endif
