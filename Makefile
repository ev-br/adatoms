include /home/br/ALPS/share/alps/include.mk       # ratatoskr
#include /usr/shared_apps/packages/alps-1.3.5-goto/share/alps/include.mk  #elysium

CXX=g++
#CXXFLAGS_METR= -I/home/br/alps-install/boost_1_41_0  # -O3  -pg     # @ ymir
#CXXFLAGS_METR= -I/home/br/alps-install/boost_1_43_0   -pg      # @ ratatoskr
#CPPFLAGS_METR= -O3 -DNDEBUG
MCXXFLAGS= -Wno-deprecated -DNDEBUG  #-O3 #-fopenmp #-UNDEBUG # -Wall #-pg


metr: metr_mc.cc 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(MCXXFLAGS) metr_mc.cc -o metr $(LDFLAGS) $(LIBS)


clean   :
	rm a.out metr ising evaluate



