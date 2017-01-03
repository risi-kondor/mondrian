#ifndef _ResourceManager
#define _ResourceManager

#include <thread>
#include <unordered_map>
#include <set>
#include <unordered_set>

#include "Mondrian_base.hpp"
#include "ThreadManager.hpp"
#include "FlushedGuard.hpp"


namespace Mondrian{

  extern ThreadManager threadManager;

  class ManagedJob;
  class ResourceManager;


  class ResourceDependency{
  public:
    ManagedJob* owner;
    int resource;
    bool fulfilled=false;
    ResourceDependency* next_in_line=nullptr;
  public:
    ResourceDependency(ManagedJob* _owner, const int r, const bool _fulfilled=false): 
      owner(_owner), resource(r), fulfilled(_fulfilled) {next_in_line=nullptr;}
  };



  class ManagedJob{
  public:

    int id;
    ResourceManager* owner;
    vector<ResourceDependency*> dependencies;
    function<void()> lambda;

  public:

    ManagedJob(ResourceManager* _owner, int r, function<void()> _lambda, int _id): 
      id(_id),
      owner(_owner),
      lambda(_lambda){
      {CoutLock lock; cout<<"Adding job number "<<id<<" [ "<<r<<" ]"<<endl;}
      registerDependency(r);
      cout<<endl;
    }

    ManagedJob(ResourceManager* _owner, const initializer_list<int>& lst, function<void()> _lambda, int _id): 
      id(_id),
      owner(_owner),
      lambda(_lambda){
      {CoutLock lock; cout<<"Adding job "<<id<<"  [ "; for(auto r:lst) cout<<r<<" "; cout<<"]"<<endl;}
      for(auto r:lst) registerDependency(r);
    }

   ManagedJob(ResourceManager* _owner, const vector<int>& lst, function<void()> _lambda, int _id): 
      id(_id),
      owner(_owner),
      lambda(_lambda){
      {CoutLock lock; cout<<"Adding job "<<id<<"  [ "; for(auto r:lst) cout<<r<<" "; cout<<"]"<<endl;}
      for(auto r:lst) registerDependency(r);
    }

    ~ManagedJob(){for(auto p: dependencies) delete p;}


  public:

    void registerDependency(const int r);
    void launch();

    bool isRunnable() const{
      for(auto r: dependencies) 
	if(!r->fulfilled) return false;
      return true;
    }

  };


    
  class ResourceManager{
  public:

    unordered_set<ManagedJob*> jobs;
    vector<thread> threads;
    unordered_map<int,ResourceDependency*> end_of_queue;

    mutex guardmx;
    mutex cvmx;
    condition_variable working_cv;
    int jobc=0;
    mutex flushedmx;

  public:

    //class FlushedGuard{
    //public:
    //  lock_guard<mutex> lock;
    //  FlushedGuard(ResourceManager& _owner): lock(_owner.guardmx){_owner.flush();}
      //FlushedGuard(FlushedGuard&& x): lock(std::move(lock)){}
    //};


  public:

    ~ResourceManager(){
      unique_lock<mutex> lk(cvmx);
      while(jobs.size()>0){
	working_cv.wait(lk);
	if(jobs.size()>0) cout<<"Spurious wakeup"<<endl;
      }
      for(auto& th:threads) th.join(); //{th->join(); delete th;}
    }


  public:

    void addJob(const int r, function<void()> lambda){
      lock_guard<mutex> flushed_guard(flushedmx);
      lock_guard<mutex> guard(guardmx);
      ManagedJob* job=new ManagedJob(this,r,lambda,jobc++);
      jobs.insert(job);
      if(job->isRunnable()) job->launch();
    }

    void addJob(const initializer_list<int>& lst, function<void()> lambda){
      lock_guard<mutex> flushed_guard(flushedmx);
      lock_guard<mutex> guard(guardmx);
      ManagedJob* job=new ManagedJob(this,lst,lambda,jobc++);
      jobs.insert(job);
      if(job->isRunnable()) job->launch();
    }

    template<class TYPE>
    TYPE addTask(const int r, function<TYPE()> lambda){
      TYPE ret;
      mutex donecv_guard;
      condition_variable donecv;
      lock_guard<mutex> flushed_guard(flushedmx);
      {lock_guard<mutex> guard(guardmx);
	ManagedJob* job=new ManagedJob(this,r,[&donecv, &donecv_guard, &ret, lambda](){
	    ret=lambda();
	    lock_guard<mutex> lk(donecv_guard);
	    donecv.notify_one();
	  },jobc++);
	jobs.insert(job);
	if(job->isRunnable()) job->launch();
      }
      unique_lock<mutex> lk(donecv_guard);
      donecv.wait(lk);
      return ret;
    }

    template<class TYPE>
    TYPE addTask(const initializer_list<int>& lst, function<TYPE()> lambda){
      TYPE ret;
      mutex donecv_guard;
      condition_variable donecv;
      lock_guard<mutex> flushed_guard(flushedmx);
      {lock_guard<mutex> guard(guardmx);
	ManagedJob* job=new ManagedJob(this,lst,[&donecv, &donecv_guard, &ret, lambda](){
	    ret=lambda();
	    lock_guard<mutex> lk(donecv_guard);
	    donecv.notify_one();
	  },jobc++);
	jobs.insert(job);
	if(job->isRunnable()) job->launch();
      }
      unique_lock<mutex> lk(donecv_guard);
      donecv.wait(lk);
      return ret;
    }

    void jobCompleteCallback(ManagedJob* job){
      lock_guard<mutex> guard(guardmx);
      //{CoutLock lock; cout<<"callback"<<endl;}
      unordered_set<ManagedJob*> dependents;
      for(auto r:job->dependencies){
	if(r->next_in_line==nullptr){
	  end_of_queue.erase(r->resource);
	}else{
	  r->next_in_line->fulfilled=true;
	  dependents.insert(r->next_in_line->owner);
	}
      }
      for(auto d: dependents) 
	if(d->isRunnable()) d->launch();
      jobs.erase(job);
      delete job;
      if(jobs.size()==0){lock_guard<mutex> lk(cvmx); working_cv.notify_one();}
    }

    void flush(){
      unique_lock<mutex> lk(cvmx);
      while(jobs.size()>0){
	working_cv.wait(lk);
	if(jobs.size()>0) cout<<"Spurious wakeup"<<endl;
      }
    }

};



inline void ManagedJob::registerDependency(const int r){
  auto it=owner->end_of_queue.find(r);
  if(it==owner->end_of_queue.end()){
    auto req=new ResourceDependency(this,r,true);
    owner->end_of_queue[r]=req;
    dependencies.push_back(req);
  }else{
    auto req=new ResourceDependency(this,r,false);
    owner->end_of_queue[r]->next_in_line=req;
    owner->end_of_queue[r]=req;
    dependencies.push_back(req);
  }
}


inline void ManagedJob::launch(){
  owner->threads.push_back(thread([this](){
				    this->lambda();
				    owner->jobCompleteCallback(this);}));
}



} // namespace Mondrian 


#endif
    /*
    mutex gate;
    atomic<int> nthreads;
    int nprivileged=0;
    int maxthreads=4;
    int maxprivileged=1;
    */
   //ResourceManager()=delete;
  
    //ResourceManager(const int _maxthreads=1000, const int _maxprivileged=1): 
    //  maxthreads(_maxthreads), maxprivileged(_maxprivileged), nthreads(0), nprivileged(0) {gate.lock();}; 

    //    bool is_ready(){return nthreads<maxthreads;}
    /*
    void printinfo(){
      CoutLock lock;
      cout<<"    (threads: "<<nthreads-nprivileged<<"+"<<nprivileged<<" local, ";
      cout<<threadManager.get_nthreads()<<" global)"<<endl;
    }
    */
  
   //template<class FUNCTION>
    //void addJob(function<void(ManagedJob&)>preamble, FUNCTION lambda){
    //  ManagedJob* th=new ManagedJob(preamble);
    //  //preamble(*this);
    //}
