#ifndef __MATR_OPS__
#define __MATR_OPS__
#include <iostream>
#include <numeric>
#include <boost/tuple/tuple.hpp>
#include "michele_matrix.h"
#include "num_util.h"

extern "C" void dsyev_( char *jobz, char *uplo, int &n, double *a, int &lda,
					   double *w, double *work, int &lwork, int &info );


class objW2{
  public:
    objW2(int alpha=2, int beta=3) { 
       if( 2!=alpha || 3!=beta){
         std::cout<<"objW2 ctor: alpha="<<alpha<<"  beta="<<beta<<"\n";
         exit(-111);
       }

       alpha_=alpha; beta_=beta;
       norm_ = 1./(2.+(double)alpha_) - 1./(2.-(double)beta_);  // $\int d^2 \mathbf{r} W_2(|\mathbf{r}|)$
       norm_ *= 2.*M_PI; 
        //std::cout<<"norm="<<norm_<<"\n";
    }
    double alpha() const {return alpha_;}
    double beta() const {return beta_;}
    double operator()(const coords_t& a1, const coords_t& a2) const{
       double rrr = atoms_distance(a1,a2);
       return rrr>1. ? 1./(rrr*rrr*rrr) : rrr*rrr ;
       //return rrr>1. ? 1./std::pow(rrr,beta_) : std::pow(rrr,alpha_) ;  // this is excruciatingly slow
    }
    double norm() const { return norm_; }
  private:
    int alpha_, beta_; 
    double norm_;
};




class det_weight_t{
  public:
    det_weight_t(int N, int alpha=1.,int beta=3.) : W2_(alpha, beta), 
                                                      A_(N,N,0), 
                                                      eigg_(N), eigg2_(N) {}
    double alpha() const { return W2_.alpha(); }
    double beta() const { return W2_.beta(); }
    double norm() const { return W2_.norm(); }

    std::pair<int,double> W(const atoms_t& atoms) const;
    std::pair<int,double> Wrat(const atoms_t& atoms_new, const atoms_t& atoms_old) const;

  private:
    objW2 W2_;
    mutable Matrix A_;
    mutable std::vector<double> eigg_, eigg2_;

    void build_matrix(const atoms_t& atoms) const;
    std::pair<int, std::vector<double> > eigenvalues(const atoms_t& atoms) const;
    std::pair<int, double > det_full(const atoms_t& atoms) const;
};


inline void det_weight_t::build_matrix(const atoms_t& atoms) const{
  assert( A_.NRows() == atoms.size() );
  int N=atoms.size();

  // laplacian 
  for (int i=0; i<N; ++i) {
    double sum=0.;
    for (int j=0; j<N; ++j) {
      if(i!=j){
       A_(i,j) = -W2_( atoms[i],atoms[j]);
       sum -= A_(i,j);
      }
    }
    A_(i,i)=sum; 
  }
 
  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {
      A_(i,j) += 1.;
    }
  }
  //std::cout<<A;
};


inline std::pair<int, std::vector<double> > det_weight_t::eigenvalues(const atoms_t& atoms) const{
   assert( A_.NRows()==atoms.size() );
   build_matrix(atoms);

   int N=atoms.size();
   int lda=N;
   int lwork = 3*N-1;
   Matrix WORK(lwork, 1);    // TODO: move to ctor?
   std::vector<double> W(N,1.);
   int info;

   char a1[]="N"; 
   char a2[]="L";

   dsyev_(a1, a2, N, &A_[0], lda, &W[0], &WORK[0], lwork, info);

   return std::make_pair(info,W);
};



inline std::pair<int, double > det_weight_t::det_full(const atoms_t& atoms) const{
  assert( A_.NRows() == atoms.size() );
  int info;
  boost::tie(info, eigg_) = eigenvalues(atoms);

  assert(eigg_.size()==atoms.size());
  if( info != 0 ){ std::cerr<<"!!!!!!!!!!! def_full: dsyev_'s info="<<info<<"\n"; }

  double det=1.;
  for(size_t j=0; j<eigg_.size(); ++j){
    det *= eigg_[j];
  }

/*
  //@#$
  if( det<1e-10){ 
   std::cout<<" *************** @det_full \n";
   std::cout<<"det = "<<det<<"\n";
   std::cout<<"eigenvalues = ";
   br_num_util::print_vector(eigg_);
  }
*/


  return std::make_pair(info,det);
};


inline std::pair<int,double> det_weight_t::W(const atoms_t& atoms) const{
  return det_full(atoms);
};



inline std::pair<int,double> det_weight_t::Wrat(const atoms_t& atoms_new, const atoms_t& atoms_old) const{
  // eig o'new
  int info;
  boost::tie(info, eigg_) = eigenvalues(atoms_new);

  assert(eigg_.size()==atoms_new.size());
  if( info != 0 ){ std::cerr<<"!!!!!!!!!!! Wrat(new): dsyev_'s info="<<info<<"\n"; }

  // eig o'ol'
  boost::tie(info, eigg2_) = eigenvalues(atoms_old);

  assert(eigg2_.size()==atoms_old.size());
  if( info != 0 ){ std::cerr<<"!!!!!!!!!!! Wrat(old): dsyev_'s info="<<info<<"\n"; }

  // ratio of det-s
  double rat=1.;
  for(size_t j=0; j<eigg_.size(); ++j){
    rat *= eigg_[j]/eigg2_[j];
  }
  return std::make_pair(info, rat);
};



#endif
