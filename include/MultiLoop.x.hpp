#ifndef _MultiLoop
#define _MultiLoop

//#include "pMMFbase.hpp"
//#include "ThreadBank.hpp"

//extern bool multithreading;


class MultiLoop{
public:

  MultiLoop(const int n, std::function<void(int)> lambda){
    for(int i=0; i<n; i++) lambda(i);
  }
  
  MultiLoop(const int n1, const int n2, std::function<void(int,int)> lambda){
    for(int i1=0; i1<n1; i1++) 
      for(int i2=0; i2<n2; i2++)
	lambda(i1,i2);
  }

};




#endif






/*
  for(int i=0; i<n; i++) workers.push_back(thread([lambda](int _i){
							nrunningthreads++;
                                                        #ifdef _MULTILOOP_VERBOSE
							  coutlock.lock();
							  cout<<"nrunningthreads="<<nrunningthreads<<endl;
							  coutlock.unlock();
                                                        #endif
							lambda(_i); 
							nrunningthreads--;
						      },i));
*/
