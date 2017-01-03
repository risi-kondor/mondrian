#ifndef _MultiLoop
#define _MultiLoop

#include "ThreadBank.hpp"

namespace Mondrian{

extern bool multithreading;


class MultiLoop{
public:

  MultiLoop(const int n, std::function<void(int)> lambda){
    if(multithreading){
      ThreadBank threads(1000);
      for(int i=0; i<n; i++) threads.add(lambda,i);
    }else for(int i=0; i<n; i++) lambda(i);
  }
  
  MultiLoop(const int n1, const int n2, std::function<void(int,int)> lambda){
    if(multithreading){
      for(int i1=0; i1<n1; i1++){
	ThreadBank threads(1000);
	for(int i2=0; i2<n2; i2++)
	  threads.add(lambda,i1,i2);
      }
    }else 
      for(int i1=0; i1<n1; i1++) 
	for(int i2=0; i2<n2; i2++)
	  lambda(i1,i2);
  }

};


template<class OBJ>
class MultiForeach{
public:

  template<class CONTAINER>
  MultiForeach(CONTAINER& V, std::function<void(OBJ)> lambda){
    if(multithreading){
      vector<thread> workers;
      for(auto& x:V) workers.push_back(thread(lambda,x));
      for(auto& w:workers) w.join();
    }else for(auto& x:V) lambda(x);
  }

};


template<class ARG1>
class MultiLoop1{ // should be able to do it with context
public:
  MultiLoop1(const int n, std::function<void(int,ARG1)> lambda, ARG1 arg1){
    if(multithreading){
      ThreadBank threads(1000);
      for(int i=0; i<n; i++) threads.add(lambda,i,arg1);
    }else for(int i=0; i<n; i++) lambda(i,arg1);
  }
};



template<class ARG1, class ARG2>
class MultiLoop2{ // should be able to do it with context
public:
  MultiLoop2(const int n, std::function<void(int,ARG1,ARG2)> lambda, ARG1 arg1, ARG2 arg2){
    if(multithreading){
      ThreadBank threads(1000);
      for(int i=0; i<n; i++) threads.add(lambda,i,arg1,arg2);
    }else for(int i=0; i<n; i++) lambda(i,arg1,arg2);
  }
};

} // namespace Mondrian

#endif




