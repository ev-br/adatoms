#ifndef __ADDR_ACPT__
#define __ADDR_ACPT__

// assumed usage:
//   addr_accpt counter ;
//   counter += do_an_update() ;            // keep track of
//   measurements["counter"] << counter() ; // checkpoint 
//   


class addr_accpt{

public:
  addr_accpt(){ _addr=0.; _accpt=0.; }
  double operator()() const { return _accpt/(_addr+1e-10); }
  void address() { ++_addr; }
  void accept() { ++(this->_accpt); }
  void clear() { _addr=0.; _accpt=0; }


  addr_accpt& operator=(const addr_accpt& rhs){ 
        _addr=rhs._addr; 
        _accpt=rhs._accpt;
        return *this ; 
  } 
  addr_accpt& operator+=(const addr_accpt& rhs){
        _addr+= rhs._addr;
        _accpt+=rhs._accpt;
        return *this ; 
  }  

private:
  double _addr, _accpt ; 

};


#endif



/******************************** a driver for checking it

#include<iostream>
#include"addr_accpt.h"

addr_accpt update_no(){
  addr_accpt tmp;
  tmp.address(); 
  return tmp ; 
}

addr_accpt update_yes(){
  addr_accpt tmp;
  tmp.address(); 
  tmp.accept();
  return tmp ; 
}



int main(){

addr_accpt X;

X+=update_no();
std::cout<<"1: "<<X()<<"\n";

X+=update_yes();
std::cout<<"2: "<<X()<<"\n";

X+=update_no();
X+=update_no();
std::cout<<"4: "<<X()<<"\n";


}


**************************************/
