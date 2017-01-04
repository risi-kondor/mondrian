#include "Mondrian_base.hpp"
#include "ThreadManager.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "LazyVectorv.hpp"


std::default_random_engine randomNumberGenerator;

SCALAR Mondrian::Vectorv::dummyZero = 0;
SCALAR Mondrian::Vectorl::dummyZero = 0;
SCALAR Mondrian::Vectorh::dummyZero = 0;
SCALAR Mondrian::LazyVectorv::dummyZero = 0;

namespace Mondrian{  
  bool multithreading=true;
  ThreadManager threadManager(4);
}

mutex cout_mutex;
mutex CoutLock::mx;


//using namespace Mondrian;
//bool Mondrian::multithreading = true;
