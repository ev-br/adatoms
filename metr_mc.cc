#include <iostream>
#include <fstream>
#include <alps/alea.h>
#include "atoms2.h"
#include "num_util.h"
#include "matr_ops2.h"
#include "updates.h"
#include "rndm.h"

using namespace std;
using namespace br_num_util;


void save_atoms(const atoms_t& atoms, std::string fname);
void load_atoms(atoms_t& atoms, std::string fname);


engine_type rng(14u);
generator_type rndm(rng,boost::uniform_01<>());


int main(int argc, char *argv[]){

  if(argc<2){
    cerr<<"Expecting the Nat at the command line arg.\n";
    exit(-1);
  }

  size_t Nat; 
  Nat=(size_t)strtoul(argv[1], 0, 0);
 
  // R & F 
  const bool tabulate=true;
  const bool verbose=false;

  objF3<tabulate>  F3(Nat,1.,verbose);
  det_weight_t www(Nat,2,3);
  pwr_weight_t<2,generator_type> w_r3(rndm);

  size_t nsweeps=(size_t)2 ; //1e4;
  size_t step_measure=(size_t)5e2;
  size_t steps_per_sweep=(size_t)1e6; 
  size_t therm_sweeps=30;

  // initial positions of the atoms
  atoms_t atoms(Nat);
  for(size_t j=0; j<Nat; ++j){
    atoms[j].x=0.1*j; atoms[j].y=0.3*j;
  }

  //load(atoms,"tr.cnf.cut.txt");
  cout<<"# W="<<www.W(atoms).second<<"\n";
  cout<<"# F="<<F3(atoms, verbose).second<<"\n";

  cout<<" F uses tab="<< (tabulate ? "true" : "false")<<"  & offset="<<F3.offset()<<"\n";
  cout<<" probab uses alpha="<<www.alpha()<<"  & beta="<<www.beta()<<"\n";
  cout<<"#Nat count    mean   error\n";

  //for the update(s)
  double dx=12;

  addr_accpt mv_unif, mv_r3, pl_one, mv_all;

  // observables
  size_t neglected=0;

  double sum=0., Z=0.;
  size_t acpt=0;
  double nrm=pow( www.norm()*atoms.size(), (double)(atoms.size()-1.) ) * atoms.size() ;
  alps::RealObservable integr;
  integr.reset(true);      

  // sliding cutoff
  short meas_max=-4, meas_min=-12;
  alps::RealVectorObservable integr_cut;
  integr_cut.reset(true);  

  alps::RealObservable dist2;
  dist2.reset(true);


  // metropolis MC
  for(size_t sweeps=0; sweeps<nsweeps; ++sweeps){
    for(size_t step=0; step<steps_per_sweep ; ++step ){    

     // updates
     mv_unif+=move_one_uniform(rndm,atoms,www,dx);
     mv_r3+=move_one_r3(rndm,atoms,www);
     //mv_all+=move_all_uniform(rndm,atoms,www,dx);
     //pl_one += place_one(rndm, atoms,www,w_r3);

      if( 0==step%step_measure){
        // measure
        bool lost_it; 
        int info;
        double term, denom, term0 ;
        boost::tie(lost_it,term0) = F3(atoms);
        boost::tie(info,denom)=www.W(atoms);

        //cout.precision(12);

        if( info !=0) {std::cerr<<" W() panic: info="<<info<<"  ";
           std::cout<<"sweeps="<<sweeps<<"  step="<<step<<"\n";
           std::cout<<"term= "<<term<<"  denom="<<denom<<"\n";
           print_vector(atoms);
           save_atoms(atoms,"troblesome_cnf.txt");
           exit(-111);
        } 


        if( lost_it ) { ++neglected;  term=0;  }   // throw away either of those!
        else{ term = term0/denom;  }

        double det_CUTOFF=1e-10;
        if(denom<det_CUTOFF){ save_atoms(atoms,"cut_cnf.txt"); }

       // sliding cutoff
       std::valarray<double> integr_slide(0.,meas_max-meas_min) ;
       for( short j=meas_max; j>meas_min; --j ){
         if( denom>pow(10.,j) ) { integr_slide[meas_max -j] = term*nrm; }
       }
       integr_cut << integr_slide;

       sum += term*nrm;
       Z +=1.;
       integr << term*nrm;

       // av distance
       double ddd=0;
       for(size_t j=0; j<atoms.size(); ++j){
         double d2=atoms_distance(atoms[0], atoms[j]);
         d2 *= d2;
         ddd += d2;
       }
       ddd /= atoms.size();
       dist2 << ddd;
      }  // if(0==step%step_measure)
    } // for(int step=...

    // printout
    cout<<"************** sweeps="<<sweeps<<"\n";
    cout<<"probab: alpha="<<www.alpha()<<"  & beta="<<www.beta()<<" F.offset()= "<<F3.offset()<<"\n";
    cout<<"N = "<<Nat<<"   naive=" << sum/Z<< "\n" ; 

    //cout<<"integr="<<integr<<"\n\n";

    cout<<"integr_cut: ( count="<<integr_cut[0].count()<<" )\n";
    for(short j=meas_max; j>meas_min; --j){
       cout<<"   "<<integr_cut[meas_max-j].mean()<<" +/- "<<integr_cut[meas_max-j].error()<<
               "   # (w/det > "<<pow(10.,j)<<" ) "<<
               "  conv: "<<alps::convergence_to_text(integr_cut[meas_max-j].converged_errors())<<
               "  tau="<<integr_cut[meas_max-j].tau()<<"\n";
    }
    print_vector(atoms);
    alps::RealObsevaluator dist(dist2);
    dist = sqrt(dist);
    cout<<"<r>="<<dist.mean()<<" "<<dist.error()<<"\n";
    cout<<"neglected : "<<1.*neglected/Z<<"  acpt(unif)="<<mv_unif()<<"  acpt(r3)="<<mv_r3()
                        <<"  acpt(place)="<<pl_one()<<" acpt(mv_all)="<<mv_all()<<"\n";


    cout<<"\n";
    save_atoms(atoms,"cnf.txt");

    mv_unif.clear();
    mv_r3.clear();
    pl_one.clear();

    if(sweeps==therm_sweeps){
      std::cout<<"#************ therm done\n";
      sum=0.;
      Z=0;
      neglected=0;
      integr.reset(true);
      dist2.reset(true);
      integr_cut.reset(true); 
    }

  }  //while

}; // main



void save_atoms(const atoms_t& atoms, std::string fname){
  ofstream outf;
  outf.open (fname.c_str());
  for(size_t j=0; j<atoms.size(); ++j){
    outf<<atoms[j].x<<"  "<< atoms[j].y<<"\n";
  }
  outf.close();
}


void load_atoms(atoms_t& atoms, std::string fname){
  ifstream inf;
  inf.open (fname.c_str());
  for(size_t j=0; j<atoms.size(); ++j){
    inf>>atoms[j].x>>atoms[j].y;
  }
  inf.close();
}


