#ifndef _Vectorh
#define _Vectorh

#include <unordered_map>

#include "Cvector.hpp"
#include "SparseVector.hpp"
#include "VectorView.hpp"
#include "GivensRotation.hpp"

extern default_random_engine randomNumberGenerator;



namespace Mondrian{
  
  class Vectorh: public SparseVector{
  public:

    unordered_map<INDEX,SCALAR> tbl;
    static SCALAR dummyZero;


  private:

    mutex insrtmx;


  public:

    Vectorh(): SparseVector(0) {}
    
    Vectorh(const int _n): SparseVector(_n) {}

    Vectorh(const initializer_list<SCALAR> list): Vectorh(list.size()){
      int i=0; for(SCALAR v:list) tbl[i++]=v;}

    Vectorh(const int _n, const initializer_list<ivpair> list): Vectorh(n){
      for(const ivpair& v:list) tbl[v.first]=v.second;}

    ~Vectorh(){}

    
  public: // copying
  
    Vectorh(const Vectorh& x): SparseVector(x.n){
      COPY_WARNING("Vectorh");
      tbl=x.tbl;
    } 

    Vectorh(const Vectorh& x, const _NoWarn dummy): SparseVector(x.n){
      tbl=x.tbl;
    } 

    Vectorh(Vectorh&& x): SparseVector(x.n){
      MOVE_WARNING("Vectorh");
      tbl=std::move(x.tbl);
    }

    Vectorh& operator=(const Vectorh& x){
      ASSIGN_WARNING("Vectorh");
      n=x.n; 
      tbl=x.tbl;
      return *this;
    }
    
    Vectorh& operator=(Vectorh&& x){
      MOVEASSIGN_WARNING("Vectorh");
      n=x.n; 
      tbl=std::move(x.tbl); 
      return *this;
    }

    Vectorh copy() const{
      Vectorh v(n); 
      v.tbl=tbl; 
      return v;
    }

    
  public: // conversions 

    Vectorh(const Vector& v): SparseVector(v.n){
      CONVERT_WARNING("SparseVector","Vectorh");
      v.for_each_filled([this](const int i, const SCALAR v){insert(i,v);});
    }

    Vectorh(const Vector& v, _NoWarn dummy): 
      SparseVector(v.n){
      v.for_each_filled([this](const int i, const SCALAR v){insert(i,v);});
    }

    operator Cvector(){
      CONVERT_WARNING("Vectorh","Cvector");
      Cvector v=Cvector::Zero(n);
      for(auto& p:tbl) v(p.first)=p.second;
      return v;
    }


public: // named constructors
    
    static Vectorh Zero(const int _n) {return Vectorh(_n);}

    static Vectorh Filled(const int _n, const SCALAR t){
      Vectorh v(_n); for(int i=0; i<_n; i++) v.tbl[i]=t; return v;}

    static Vectorh Uniform(const int _n){
      Vectorh v(_n); 
      uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.tbl[i]=distr(randomNumberGenerator); 
      return v;}

    static Vectorh Gaussian(const int _n){
      Vectorh v(_n); 
      normal_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.tbl[i]=distr(randomNumberGenerator); 
      return v;}

    static Vectorh Bernoulli(const int _n, const double p=0.5){
      Vectorh v(_n); 
      bernoulli_distribution distr(p);
      for(int i=0; i<_n; i++) if(distr(randomNumberGenerator)==1) v.tbl[i]=1.0; 
      return v;}


  public: // polymorphism

    virtual string classname() const {return "Vectorh";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer<Vectorh>();}


  public: // comparators

    bool operator==(const Vectorh& x) const {
      for(auto& p:tbl) if(p.second!=const_cast<unordered_map<INDEX,SCALAR>&>(x.tbl)[p.first]) return false; 
      for(auto& p:x.tbl) if(p.second!=const_cast<unordered_map<INDEX,SCALAR>&>(tbl)[p.first]) return false; 
      return true;
    }

    
  public: // element access 

    SCALAR operator()(const int i) const {assert(i<n);
      auto it=tbl.find(i); if(it==tbl.end()) return 0; else return it->second;}
    
    SCALAR read(const int i) const {assert(i<n);
      auto it=tbl.find(i); if(it==tbl.end()) return 0; else return it->second;}
    

    SCALAR& operator()(const int i) {assert(i<n); 
      return tbl[i];}
    
    void set(const int i, const SCALAR v) {assert(i<n); 
      tbl[i]=v;}

    void set_msafe(const int i, const SCALAR v) {assert(i<n);
      lock_guard<mutex> lock(insrtmx);
      tbl[i]=v;}

    SCALAR* ptr(const int i){assert(i<n);
      auto it=tbl.find(i);
      if(it!=tbl.end()) return &it->second;
      else return &tbl.insert({i,0}).first->second;
    }

    bool isFilled(const int i) const {assert(i<n); 
      auto it=tbl.find(i); return (it!=tbl.end());}
    
    int nFilled() const {return tbl.size();}
    
    
  public: //iterators

    void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      for(auto& p:tbl) lambda(p.first,p.second);}

    void for_each_filled(std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(auto& p:tbl) lambda(p.first,p.second);}
    
    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const TYPE, const TYPE)> accumulator, const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(auto& p:tbl)
	t=accumulator(t,p.second);
      return t;
    }

    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const SCALAR&)> lambda, 
      std::function<TYPE(const TYPE, const TYPE)> accumulator=std::plus<TYPE>(), 
      const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(auto& p:tbl)
	t=accumulator(t,lambda(p.second));
      return t;
    }

    INDEX find_best(std::function<bool(const SCALAR&, const SCALAR&)> selector) const{
      assert(tbl.size()>0);
      INDEX besti=tbl.begin()->first;
      SCALAR bestv=tbl.begin()->second;
      for(auto& p:tbl)
	if(selector(bestv,p.second)){besti=p.first; bestv=p.second;}
      return besti;
    }

    
  public: // sparse vector methods

    SCALAR* findptr(const int i){
      auto it=tbl.find(i); 
      if(it!=tbl.end()) return &it->second;
      return &dummyZero;
    }

    const SCALAR* ptr_if_filled(const int i) const{
      auto it=tbl.find(i); 
      if(it!=tbl.end()) return &it->second;
      return nullptr;
    }

    SCALAR* ptr_if_filled(const int i){
      auto it=tbl.find(i); 
      if(it!=tbl.end()) return &it->second;
      return nullptr;
    }

    void insert(const int i, const SCALAR value) {assert(i<n); 
      tbl[i]=value;}
    
    void append(const int i, const SCALAR value) {assert(i<n); 
      tbl[i]=value;}
    
    Vectorh& merge(vector<ivpair>& x){
      for(auto& p:x) tbl[p.first]=p.second;
      return *this;
    }

    void zero(const int i){assert(i<n); 
      tbl.erase(i);}

    void tidy(){
      for(auto it=tbl.begin(); it!=tbl.end(); it++) 
	if(it->second==0) it=tbl.erase(it);}


  public: // scalar-valued operations 

    int nnz() const{
      return tbl.size();}
    
    SCALAR max() const{
      if(tbl.size()==0) return 0;
      SCALAR max=tbl.begin()->second;
      for(auto p:tbl) if(p.second>max) max=p.second;
      return max;
    }
    
    SCALAR max_abs() const{
      if(tbl.size()==0) return 0;
      SCALAR max=fabs(tbl.begin()->second);
      for(auto p:tbl) if(fabs(p.second)>max) max=fabs(p.second);
      return max;
    }
    
    int argmax() const{
      if(tbl.size()==0) return 0;
      int best=tbl.begin()->first; SCALAR max=tbl.begin()->second;
      for(auto p:tbl) if(p.second>max) {best=p.first; max=p.second;}
      return best;
    }
    
    int argmax_abs() const{
      if(tbl.size()==0) return 0;
      int best=tbl.begin()->first; SCALAR max=fabs(tbl.begin()->second);
      for(auto p:tbl) if(fabs(p.second)>max) {best=p.first; max=fabs(p.second);}
      return best;
    }
    
    SCALAR min() const{
      if(tbl.size()==0) return 0;
      SCALAR min=tbl.begin()->second;
      for(auto p:tbl) if(p.second<min) min=p.second;
      return min;
    }
    
    SCALAR min_abs() const{
      if(tbl.size()==0) return 0;
      SCALAR min=fabs(tbl.begin()->second);
      for(auto p:tbl) if(fabs(p.second)<min) min=fabs(p.second);
      return min;
    }
    
    int argmin() const{
      if(tbl.size()==0) return 0;
      int best=tbl.begin()->first; SCALAR min=tbl.begin()->second;
      for(auto p:tbl) if(p.second<min) {best=p.first; min=p.second;}
      return best;
    }
    
    int argmin_abs() const{
      if(tbl.size()==0) return 0;
      int best=tbl.begin()->first; SCALAR min=fabs(tbl.begin()->second);
      for(auto p:tbl) if(fabs(p.second)<min) {best=p.first; min=fabs(p.second);}
      return best;
    }
    
    SCALAR sum() const{
      SCALAR t=0; for(auto& p:tbl) t+=p.second; return t;}

    SCALAR norm1() const{
      SCALAR t=0; for(auto& p:tbl) t+=fabs(p.second); return t;}
    
    SCALAR norm2() const{
      SCALAR t=0; for(auto& p:tbl) t+=p.second*p.second; return t;}
    
    SCALAR diff2(const Vectorh& x) const{ // check!
      assert(x.n==n); SCALAR t=0;
      for(auto& xp:x.tbl) if(xp.second!=0) t+=(xp.second-read(xp.first))*(xp.second-read(xp.first));
      for(auto& p:tbl) if(x.read(p.first)==0) t+=p.second*p.second;
      return t;
    }
    
    SCALAR diff2(const Cvector& x) const{
      assert(x.n==n); SCALAR t=0;
      for(auto& p:tbl) t+=(p.second-x(p.first))*(p.second-x(p.first));
      for(int i=0; i<n; i++) 
	if(!isFilled(i)) t+=x(i)*x(i);
      return t;
    }
    

  public: // scalar valued arithmetic

    SCALAR dot(const Vectorh& x) const{
      SCALAR t=0; for(auto& p:x.tbl) t+=p.second*this->read(p.first); return t;}

    
  public: // in-place operations 
    
    Vectorh& operator*=(const SCALAR c){
      for(auto& p:tbl) p.second*=c; return *this;}
    Vectorh& operator/=(const SCALAR c){
      for(auto& p:tbl) p.second/=c; return *this;}

    template<class VECTOR>
    Vectorh& operator*=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:tbl) p.second*=x(p.first); return *this;}
    template<class VECTOR>
    Vectorh& operator/=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:tbl) p.second/=x(p.first); return *this;
    }
    
    Vectorh& add(const Vectorh& x, const SCALAR c=1){
      for(auto& xp:x.tbl) tbl[xp.first]+=c*xp.second; return *this;}

    Vectorh& operator+=(const Vectorh& x){return add(x,1);}
    Vectorh& operator-=(const Vectorh& x){return add(x,-1);}

    Vectorh& operator*=(const Vectorh& x){
      auto it=x.tbl.begin();
      for(auto& p:tbl){
	while(it!=x.tbl.end() && it->first<p.first) it++;
	if(it!=x.tbl.end() && it->first==p.first) p.second*=it->second;
	else p.second=0;
      }
      return *this; 
    }
    
    Vectorh& operator/=(const Vectorh& x){
      auto it=x.tbl.begin();
      for(auto& p:tbl){
	while(it!=x.tbl.end() && it->first<p.first) it++;
	if(it!=x.tbl.end() && it->first==p.first) p.second/=it->second;
	else p.second=1.0/0.0;
      }
      return *this; 
    }
    

  public: // Givens rotations 
  
  Vectorh& apply(const GivensRotation& Q){
    assert(Q.i1<n); assert(Q.i2<n);
    SCALAR* p1=findptr(Q.i1); 
    SCALAR* p2=findptr(Q.i2);
    if(p1==&dummyZero){
      if(p2==&dummyZero) return *this; 
      SCALAR x2=*p2;
      *p2=Q.cos*x2;
      insert(Q.i1,-Q.sin*x2);
    }else{
      SCALAR x1=*p1;
      if(p2==&dummyZero){
	*p1=Q.cos*x1;
	insert(Q.i2,Q.sin*x1);
      }else{
	*p1=Q.cos*x1-Q.sin*(*p2);
	*p2=Q.sin*x1+Q.cos*(*p2);
     }
    }
    return *this;
  }

  Vectorh& applyT(const GivensRotation& Q){
    assert(Q.i1<n); assert(Q.i2<n);
    SCALAR* p1=findptr(Q.i1); 
    SCALAR* p2=findptr(Q.i2);
    if(p1==&dummyZero){
      if(p2==&dummyZero) return *this; 
      SCALAR x2=*p2;
      *p2=Q.cos*x2;
      insert(Q.i1,Q.sin*x2);
    }else{
      SCALAR x1=*p1;
      if(p2==&dummyZero){
	*p1=Q.cos*x1;
	insert(Q.i2,-Q.sin*x1);
      }else{
	*p1=Q.cos*x1+Q.sin*(*p2);
	*p2=-Q.sin*x1+Q.cos*(*p2);
     }
    }
    return *this;
  }



  public: // KpointOp
    
    template<int k>
    Vectorh& apply(const KpointOp<k>& Q){ 
      int i=0; for(;i<k;i++) if(tbl.find(Q.map(i))!=tbl.end()) break;
      if(i==k) return *this;
      SCALAR* vptr[k];
      SCALAR temp[k];
      for(int i=0; i<k; i++){
	auto it=tbl.find(Q.map(i));
	if(it==tbl.end())
	  vptr[i]=&tbl.insert({Q.map(i),0}).first->second;
	else 
	  vptr[i]=&it->second;;
	temp[i]=*vptr[i];
      }
      for(int i=0; i<k; i++){
	SCALAR s=0; 
	for(int j=0; j<k; j++) 
	  s+=Q.q[j*k+i]*temp[j]; 
	*vptr[i]=s;
      }
      return *this;
    }
     

    template<int k>
    Vectorh& applyT(const KpointOp<k>& Q){ 
      int i=0; for(;i<k;i++) if(tbl.find(Q.map(i))!=tbl.end()) break;
      if(i==k) return *this;
      SCALAR* vptr[k];
      SCALAR temp[k];
      for(int i=0; i<k; i++){
	auto it=tbl.find(Q.map(i));
	if(it==tbl.end())
	  vptr[i]=&tbl.insert({Q.map(i),0}).first->second;
	else 
	  vptr[i]=&it->second;;
	temp[i]=*vptr[i];
      }
      for(int i=0; i<k; i++){
	SCALAR s=0; 
	for(int j=0; j<k; j++) 
	  s+=Q.q[i*k+j]*temp[j]; 
	*vptr[i]=s;
      }
      return *this;
    }
     

  public: // I/O 
    
    //string str(const Sparse dummy) const;
    //string str(const Dense dummy) const {return Vector::str(Dense());}
    //string str() const{return Vector::str();}
    
    
  public: // Python interface
    
    Vectorh(double* numpyInDblArray, int numpyInSize): Vectorh(numpyInSize){
      for(int i=0; i<n; i++) {
	double v=numpyInDblArray[i];
	if(v!=0) tbl[i]=v;
      }
    }
    
    void np(double** numpyOutArray1, int* numpyOutLen1){
      *numpyOutLen1=n;
      *numpyOutArray1=new double[n];
      std::fill(*numpyOutArray1,*numpyOutArray1+n,0);
      for(auto& p:tbl) (*numpyOutArray1)[p.first]=p.second; 
    }
  

  public:  
    

    //bool equals(const Vector& x) const {return (*this)==(static_cast<const Vectorl&>(x));} // needs work!

    
  };
  
}




#endif



// ----------------------------------------------------------------------------------------------------------


    //bool isDense() const {return false;}
    //public: // attributes

    //bool isSparse() const {return true;}


    /*
  public: // KpointOperators

  Vectorh& apply(const KpointOperator& Q){ // it is not assumed that the indices in Q are sorted
    vector<int> missing;
    for(int i=0; i<Q.k; i++) if(!isFilled(Q.map(i))) missing.push_back(Q.map(i)); 
    if(missing.size()==Q.k) return *this;
    for(int i=0; i<missing.size(); i++) insert(missing[i],0);
    for(int i=0; i<Q.k; i++){
      auto it=tbl.find(Q.map(i));
      Q.tempp[i]=&it->second; 
      Q.temp[i]=it->second;
    }
    for(int i=0; i<Q.k; i++){
      double s=0; 
      for(int j=0; j<Q.k; j++) s+=Q.q[j*Q.k+i]*Q.temp[j]; 
      *Q.tempp[i]=s;
    }
    return *this;
  }


  Vectorh& applyT(const KpointOperator& Q){ // it is not assumed that the indices in Q are sorted
    vector<int> missing;
    for(int i=0; i<Q.k; i++) if(!isFilled(Q.map(i))) missing.push_back(Q.map(i)); 
    if(missing.size()==Q.k) return *this;
    for(int i=0; i<missing.size(); i++) insert(missing[i],0);
    for(int i=0; i<Q.k; i++){
      auto it=tbl.find(Q.map(i));
      Q.tempp[i]=&it->second; 
      Q.temp[i]=it->second;
    }
    for(int i=0; i<Q.k; i++){
      double s=0; 
      for(int j=0; j<Q.k; j++) s+=Q.q[i*Q.k+j]*Q.temp[j]; 
      *Q.tempp[i]=s;
    }
    return *this;
  }
    */
