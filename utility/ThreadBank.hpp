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


#ifndef _ThreadBank
#define _ThreadBank

#include <thread>
#include "Mondrian_base.hpp"


namespace Mondrian{

  class ThreadManager;
  extern ThreadManager threadManager;

  class ThreadBank{
  public:

    ThreadBank()=delete;
  
    ThreadBank(const int _maxthreads=1000, const int _maxprivileged=1): 
      maxthreads(_maxthreads), maxprivileged(_maxprivileged), nthreads(0), nprivileged(0) {gate.lock();}; 

    ~ThreadBank(){for(auto& th:threads) th.join();}


  public:

    template<class FUNCTION, class OBJ>
    void add(FUNCTION lambda, const OBJ x);
  
    template<class FUNCTION, class OBJ1, class OBJ2>
    void add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2);

    template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3>
    void add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3);

    template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3, class OBJ4>
    void add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3, const OBJ4 x4);

    template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3, class OBJ4, class OBJ5>
    void add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3, const OBJ4 x4, const OBJ5 x5);

    bool is_ready(){return nthreads<maxthreads;}

    void printinfo();

  
  public:

    mutex mx;
    mutex gate;
    atomic<int> nthreads;
    int nprivileged=0; //                                               only to be touched by threadManager
    int maxthreads=4;
    int maxprivileged=1;

    vector<thread> threads;

  };

}


#include "ThreadManager.hpp"

using namespace Mondrian;


template<class FUNCTION, class OBJ>
inline void ThreadBank::add(FUNCTION lambda, const OBJ x){
  lock_guard<mutex> lock(mx); //                                   unnecessary if called from a single thread
  threadManager.enqueue(this);
  gate.lock(); //                                                  gate can only be unlocked by threadManager
  nthreads++;
  threads.push_back(thread([this,lambda](OBJ _x){lambda(_x); nthreads--; threadManager.release(this);},x));
#ifdef _THREADBANKVERBOSE
  printinfo();
#endif
}


template<class FUNCTION, class OBJ1, class OBJ2>
inline void ThreadBank::add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2){
  lock_guard<mutex> lock(mx);
  threadManager.enqueue(this);
  gate.lock();
  nthreads++;
  threads.push_back(thread([this,lambda](OBJ1 _x1, OBJ2 _x2){
	lambda(_x1,_x2); nthreads--; threadManager.release(this);},x1,x2));
#ifdef _THREADBANKVERBOSE
  printinfo();
#endif
}


template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3>
inline void ThreadBank::add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3){
  lock_guard<mutex> lock(mx);
  threadManager.enqueue(this);
  gate.lock();
  nthreads++;
  threads.push_back(thread([this,lambda](OBJ1 _x1, OBJ2 _x2, OBJ3 _x3){
	lambda(_x1,_x2,_x3); nthreads--; threadManager.release(this);},x1,x2,x3));
#ifdef _THREADBANKVERBOSE
  printinfo();
#endif
}

template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3, class OBJ4>
inline void ThreadBank::add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3, const OBJ4 x4){
  lock_guard<mutex> lock(mx);
  threadManager.enqueue(this);
  gate.lock();
  nthreads++;
  threads.push_back(thread([this,lambda](OBJ1 _x1, OBJ2 _x2, OBJ3 _x3, OBJ4 _x4){
	lambda(_x1,_x2,_x3,_x4); nthreads--; threadManager.release(this);},x1,x2,x3,x4));
#ifdef _THREADBANKVERBOSE
  printinfo();
#endif
}

template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3, class OBJ4, class OBJ5>
inline void ThreadBank::add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3, const OBJ4 x4, const OBJ5 x5){
  lock_guard<mutex> lock(mx);
  threadManager.enqueue(this);
  gate.lock();
  nthreads++;
  threads.push_back(thread([this,lambda](OBJ1 _x1, OBJ2 _x2, OBJ3 _x3, OBJ4 _x4, OBJ5 _x5){
	lambda(_x1,_x2,_x3,_x4,_x5); nthreads--; threadManager.release(this);},x1,x2,x3,x4,x5));
#ifdef _THREADBANKVERBOSE
  printinfo();
#endif
}

inline void ThreadBank::printinfo(){
  CoutLock lock;
  cout<<"    (threads: "<<nthreads-nprivileged<<"+"<<nprivileged<<" local, ";
  cout<<threadManager.get_nthreads()<<" global)"<<endl;
}



#endif
