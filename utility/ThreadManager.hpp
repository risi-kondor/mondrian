#ifndef _ThreadManager
#define _ThreadManager

#include <list>
#include <queue>
#include <mutex>
#include <thread>

#include "Mondrian_Base.hpp"


namespace Mondrian{

class ThreadBank;

class ThreadManager{
public:

  ThreadManager(const int _maxthreads):maxthreads(_maxthreads),nthreads(0){}
  ~ThreadManager(){}

public:
  
  void enqueue(ThreadBank* bank);
  void release(ThreadBank* bank);

  int get_nthreads(){lock_guard<mutex> lock(mx); return nthreads;}

private:

  bool is_runnable(ThreadBank* bank);
  void launch(ThreadBank* bank);

public:

  int maxthreads;

private:

  mutex mx;
  int nthreads;
  list<ThreadBank*> queue;

};

}


#include "ThreadBank.hpp"


using namespace Mondrian;


inline void ThreadManager::enqueue(ThreadBank* bank){
  lock_guard<mutex> lock(mx);
  if(is_runnable(bank)) launch(bank);
  else queue.push_back(bank);
}


inline void ThreadManager::release(ThreadBank* bank){
  lock_guard<mutex> lock(mx);
  if(bank->nprivileged>0) bank->nprivileged--;
  else nthreads--;
  for(auto it=queue.begin(); it!=queue.end(); it++)
    if(is_runnable(*it)){
      launch(*it);
      it=queue.erase(it);
    }
}


inline bool ThreadManager::is_runnable(ThreadBank* bank){
  return bank->is_ready() && (bank->nprivileged<bank->maxprivileged || nthreads<maxthreads) ;
}


inline void ThreadManager::launch(ThreadBank* bank){
  if(bank->nprivileged<bank->maxprivileged) bank->nprivileged++;
  else nthreads++;
  bank->gate.unlock();
}



#endif
