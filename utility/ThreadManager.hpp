/* ---------------------------------------------------------------------------
 
  Mondrian: open source multithreaded matrix library 

  Copyright (C) 2017 Imre Risi Kondor
  Copyright (C) 2015 Imre Risi Kondor, Nedelina Teneva, Pramod K Mudrakarta

  Parts of the following code are derived from the pMMF library 
  (https://github.com/risi-kondor/pMMF) which is licensed under the 
  GNU Public License, version 3. This code therefore is also licensed 
  under the terms of the GNU Public License, version 3. 
 
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 3
  of the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

 --------------------------------------------------------------------------- */


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
