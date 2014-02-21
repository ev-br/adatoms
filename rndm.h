#ifndef __RNDM__
#define __RNDM__
#include <boost/random.hpp>
#include <iostream>

typedef boost::mt19937 engine_type; 
typedef boost::variate_generator< engine_type&, boost::uniform_01<> > generator_type ;



#endif
