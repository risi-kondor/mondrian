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


#ifndef _Vectorl
#define _Vectorl

#include <list>

#include "Cvector.hpp"
#include "SparseVector.hpp"
#include "VectorView.hpp"
#include "GivensRotation.hpp"

extern default_random_engine randomNumberGenerator;



namespace Mondrian{
  
  class Vectorl: public SparseVector{
  public:

    list<ivpair> lst;
    static SCALAR dummyZero;


  private:

    mutex insrtmx;


  public:

    Vectorl(): SparseVector(0) {}

    Vectorl(const int _n): SparseVector(_n) {}

    Vectorl(const initializer_list<SCALAR> list): Vectorl(list.size()){
      int i=0; for(SCALAR v:list) lst.push_back(ivpair(i++,v));}

    Vectorl(const int _n, const initializer_list<ivpair> list): Vectorl(n){
      for(const ivpair& v:list) lst.push_back(v);}

    ~Vectorl(){}


  public: // copying

    Vectorl(const Vectorl& x): SparseVector(x.n){
      COPY_WARNING("Vectorl");
      lst=x.lst;
    } 

    Vectorl(const Vectorl& x, const _NoWarn dummy): SparseVector(x.n){
      lst=x.lst;
    } 

    Vectorl(Vectorl&& x): SparseVector(x.n){
      MOVE_WARNING("Vectorl");
      lst=std::move(x.lst);
    }

    Vectorl& operator=(const Vectorl& x){
      ASSIGN_WARNING("Vectorl");
      n=x.n; 
      lst=x.lst;
      return *this;
    }
    
    Vectorl& operator=(Vectorl&& x){
      MOVEASSIGN_WARNING("Vectorl");
      n=x.n; 
      lst=std::move(x.lst); 
      return *this;
    }

    Vectorl copy() const{
      Vectorl v(n); 
      v.lst=lst; 
      return v;
    }

    
  public: // named constructors
    
    static Vectorl Zero(const int _n) {return Vectorl(_n);}

    static Vectorl Filled(const int _n, const SCALAR t){
      Vectorl v(_n); for(int i=0; i<_n; i++) v.lst.push_back(ivpair(i,t)); return v;}

    static Vectorl Uniform(const int _n){
      Vectorl v(_n); 
      //default_random_engine randomNumberGenerator;
      uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.lst.push_back(ivpair(i,distr(randomNumberGenerator))); 
      return v;}

    static Vectorl Gaussian(const int _n){
      Vectorl v(_n); 
      //default_random_engine randomNumberGenerator;
      normal_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.lst.push_back(ivpair(i,distr(randomNumberGenerator))); 
      return v;}

    static Vectorl Bernoulli(const int _n, const double p=0.5){
      Vectorl v(_n); 
      //default_random_engine randomNumberGenerator;
      bernoulli_distribution distr(p);
      for(int i=0; i<_n; i++) if(distr(randomNumberGenerator)==1) v.lst.push_back(ivpair(i,1.0)); 
      return v;}


  public: // polymorphism

    virtual string classname() const {return "Vectorl";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer<Vectorl>();}


  public: // comparators

    bool operator==(const Vectorl& x) const {
      auto xit=x.lst.begin();
      for(auto& v:lst){
	if(v.second==0) continue;
	for(;xit!=x.lst.end() && xit->second==0; xit++){}
	if(xit->first!=v.first) return false;
	if(xit->second!=v.second) return false;
      }
      for(;xit!=x.lst.end(); xit++)
	if(xit->second!=0) return false; 
      return true;
    }

    
  public: // element access 

    SCALAR operator()(const int i) const{assert(i<n);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}
      if(it==lst.end() || it->first>i) return 0; else return it->second;
    }

    SCALAR read(const int i) const{assert(i<n);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}
      if(it==lst.end() || it->first>i) return 0; else return it->second;
    }


    SCALAR& operator()(const int i){assert(i<n);
     auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;} 
     if(it==lst.end() || it->first>i) it=lst.insert(it,ivpair(i,0));
     return it->second;
    }
    
    void set(const int i, const SCALAR v) {assert(i<n);
      //auto lb=lower_bound(lst.begin(),lst.end(),ivpair(i,v),ivpairComparator);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;} 
      if(it==lst.end()||it->first>i) lst.insert(it,ivpair(i,v));
      else it->second=v;
    }

    void set_msafe(const int i, const SCALAR v) {assert(i<n);
      lock_guard<mutex> lock(insrtmx);
      set(i,v);}

    SCALAR* ptr(const int i){assert(i<n);
      //auto lb=lower_bound(lst.begin(),lst.end(),ivpair(i,0),ivpairComparator);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}    
      if(it!=lst.end()&&it->first==i) return &it->second;
      else return &lst.insert(it,ivpair(i,0))->second;
    }

    const SCALAR* ptr_if_filled(const int i) const {assert(i<n);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}    
      if(it!=lst.end()&&it->first==i) return &it->second;
      else return nullptr;
    }

    SCALAR* ptr_if_filled(const int i){assert(i<n);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}    
      if(it!=lst.end()&&it->first==i) return &it->second;
      else return nullptr;
    }


    bool isFilled(const int i)const {assert(i<n);
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}
      return (it!=lst.end() && it->first==i);}
    
    int nFilled() const {return lst.size();}


  public: // function mappings

    void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      for(auto& p:lst) lambda(p.first,p.second);}

    void for_each_filled(std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(auto& p:lst) lambda(p.first,p.second);}
    
    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const TYPE, const TYPE)> accumulator, const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(auto& p:lst)
	t=accumulator(t,p.second);
      return t;
    }

    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const SCALAR&)> lambda, 
      std::function<TYPE(const TYPE, const TYPE)> accumulator=std::plus<TYPE>(), 
      const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(auto& p:lst)
	t=accumulator(t,lambda(p.second));
      return t;
    }

    INDEX find_best(std::function<bool(const SCALAR&, const SCALAR&)> selector) const{
      assert(lst.size()>0);
      INDEX besti=lst.begin()->first;
      SCALAR bestv=lst.begin()->second;
      for(auto& p:lst)
	if(selector(bestv,p.second)){besti=p.first; bestv=p.second;}
      return besti;
    }


  public: // sparse vector methods

    SCALAR* findptr(const int i){
      auto it=lst.begin(); while(it!=lst.end() && it->first<i) {it++;} 
      if(it!=lst.end() && it->first==i) return &it->second;
      return &dummyZero;
    }

    void insert(const int i, const SCALAR value) {assert(i<n); (*this)(i)=value;}

    void append(const int i, const SCALAR value) {assert(i<n); lst.push_back(ivpair(i,value));}

    Vectorl& merge(vector<ivpair>& x){
      std::sort(x.begin(),x.end(),ivpairComparator);
      auto it=lst.begin();
      for(auto& p:x){
	auto lb=lower_bound(it,lst.end(),p,ivpairComparator);
	if(lb==lst.end()) {lst.push_back(p); continue;}
	if(lb->first>p.first) lst.insert(lb,p);
	else lb->second=p.second;
	it=lb;
      }
      return *this;
    }
    
    void zero(const int i){assert(i<n); 
      auto it=lst.begin(); while(it!=lst.end() && it->first<i){it++;}
      if(it!=lst.end() && it->first==i) lst.erase(it);}
    
    void tidy() {int i=-1;
      for(auto it=lst.begin(); it!=lst.end(); it++){
	while(it->first==i || it->second==0) it=lst.erase(it); 
	if(it!=lst.end()) i=it->first;}
    }


  public: // conversions 

    Vectorl(const SparseVector& v): SparseVector(v.n){
      CONVERT_WARNING("SparseVector","Vectorl");
      v.sort(); v.for_each_filled([this](const int i, const SCALAR v){append(i,v);});
    }

    Vectorl(const SparseVector& v, const _NoWarn dummy): 
      SparseVector(v.n){
      v.sort(); v.for_each_filled([this](const int i, const SCALAR v){append(i,v);});
    }

    operator Cvector(){
      CONVERT_WARNING("Vectorl","Cvector");
      Cvector v=Cvector::Zero(n);
      for(auto& p:lst) v(p.first)=p.second;
      return v;
    }

    
  public: // scalar-valued operations 

    int nnz() const {const_cast<Vectorl&>(*this).tidy(); return lst.size();}
    
    SCALAR max() const{
      if(lst.size()==0) return 0;
      SCALAR max=lst.begin()->second;
      for(auto p:lst) if(p.second>max) max=p.second;
      return max;
    }
    
    SCALAR max_abs() const{
      if(lst.size()==0) return 0;
      SCALAR max=fabs(lst.begin()->second);
      for(auto p:lst) if(fabs(p.second)>max) max=fabs(p.second);
      return max;
    }
    
    int argmax() const{
      if(lst.size()==0) return 0;
      int best=lst.begin()->first; SCALAR max=lst.begin()->second;
      for(auto p:lst) if(p.second>max) {best=p.first; max=p.second;}
      return best;
    }
    
    int argmax_abs() const{
      if(lst.size()==0) return 0;
      int best=lst.begin()->first; SCALAR max=fabs(lst.begin()->second);
      for(auto p:lst) if(fabs(p.second)>max) {best=p.first; max=fabs(p.second);}
      return best;
    }
    
    SCALAR min() const{
      if(lst.size()==0) return 0;
      SCALAR min=lst.begin()->second;
      for(auto p:lst) if(p.second<min) min=p.second;
      return min;
    }
    
    SCALAR min_abs() const{
      if(lst.size()==0) return 0;
      SCALAR min=fabs(lst.begin()->second);
      for(auto p:lst) if(fabs(p.second)<min) min=fabs(p.second);
      return min;
    }
    
    int argmin() const{
      if(lst.size()==0) return 0;
      int best=lst.begin()->first; SCALAR min=lst.begin()->second;
      for(auto p:lst) if(p.second<min) {best=p.first; min=p.second;}
      return best;
    }
    
    int argmin_abs() const{
      if(lst.size()==0) return 0;
      int best=lst.begin()->first; SCALAR min=fabs(lst.begin()->second);
      for(auto p:lst) if(fabs(p.second)<min) {best=p.first; min=fabs(p.second);}
      return best;
    }
    
    SCALAR sum() const{
      SCALAR t=0; for(auto& p:lst) t+=p.second; return t;}

    SCALAR norm1() const{
      SCALAR t=0; for(auto& p:lst) t+=fabs(p.second); return t;}
    
    SCALAR norm2() const{
      SCALAR t=0; for(auto& p:lst) t+=p.second*p.second; return t;}
    
    SCALAR diff2(const Vectorl& x) const {
      assert(x.n==n); SCALAR t=0; 
      auto it=lst.begin();
      for(auto& xv:x.lst){
	while(it!=lst.end() && it->first<xv.first) {t+=(it->second)*(it->second); it++;}
	if(it!=lst.end() && it->first==xv.first) {t+=(it->second-xv.second)*(it->second-xv.second);it++;}
	else {t+=xv.second*xv.second;}
      }
      return t; 
    }

    SCALAR diff2(const Cvector& x) const {
      assert(x.n==n); SCALAR t=0; 
      int ix=0;
      for(auto& p:lst){
	for(int i=ix; i<p.first; i++) t+=x.array[i]*x.array[i];
	t+=(p.second-x.array[p.first])*(p.second-x.array[p.first]);
	ix=p.first+1;	
      }
      for(int i=ix; i<n; i++) t+=x.array[i]*x.array[i];
      return t; 
    }


  public: // scalar valued arithmetic

    SCALAR dot(const Vectorl& x) const{
      SCALAR t=0; 
      auto xit=x.lst.begin();
      for(auto& v:lst){
	int i=v.first;
	while(xit!=x.lst.end() && xit->first<i) xit++;
	if(xit==x.lst.end()) break;
	if(xit->first==i) t+=v.second*xit->second;
      }
      return t;
    }
    
    
  public: // in-place operations 
    
    Vectorl& operator*=(const SCALAR c){
      for(auto& p:lst) p.second*=c; return *this;}
    Vectorl& operator/=(const SCALAR c){
      for(auto& p:lst) p.second/=c; return *this;}

    template<class VECTOR>
    Vectorl& operator*=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:lst) p.second*=x(p.first); return *this;}
    template<class VECTOR>
    Vectorl& operator/=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:lst) p.second/=x(p.first); return *this;
    }
    
    Vectorl& add(const Vectorl& x, const SCALAR c=1){
      auto it=lst.begin();
      for(auto& xv:x.lst){
	while(it!=lst.end() && it->first<xv.first) it++;
	if(it!=lst.end() && it->first==xv.first) it->second+=c*xv.second;
	else {lst.insert(it,ivpair(xv.first,c*xv.second));}
      }
      return *this;
    }

    Vectorl& operator+=(const Vectorl& x){return add(x,1);}
    Vectorl& operator-=(const Vectorl& x){return add(x,-1);}

    Vectorl& operator*=(const Vectorl& x){
      auto it=x.lst.begin();
      for(auto& p:lst){
	while(it!=x.lst.end() && it->first<p.first) it++;
	if(it!=x.lst.end() && it->first==p.first) p.second*=it->second;
	else p.second=0;
      }
      return *this; 
    }
    
    Vectorl& operator/=(const Vectorl& x){
      auto it=x.lst.begin();
      for(auto& p:lst){
	while(it!=x.lst.end() && it->first<p.first) it++;
	if(it!=x.lst.end() && it->first==p.first) p.second/=it->second;
	else p.second=1.0/0.0;
      }
      return *this; 
    }
    

  public: // Givens rotations 

    Vectorl& apply(const GivensRotation& Q){
      assert(Q.i1<n); assert(Q.i2<n);
      if(Q.i1>Q.i2) return applyT(Q);
      auto it1=lst.begin(); 
      while(it1!=lst.end() && it1->first<Q.i1) {it1++;} 
      if(it1->first==Q.i1){
	SCALAR x1=it1->second;
	auto it2=it1; it2++; 
	while(it2!=lst.end() && it2->first<Q.i2) {it2++;} 
	if(it2->first==Q.i2){
	  it1->second=Q.cos*x1-Q.sin*(it2->second);
	  it2->second=Q.sin*x1+Q.cos*(it2->second);
	}else{
	  it1->second=Q.cos*x1;
	  lst.insert(it2,ivpair(Q.i2,Q.sin*x1));
	}
      }else{
	auto it2=it1; it2++; 
	while(it2!=lst.end() && it2->first<Q.i2) {it2++;} 
	if(it2->first==Q.i2){
	  SCALAR x2=it2->second;
	  it2->second=Q.cos*x2;
	  lst.insert(it1,ivpair(Q.i1,-Q.sin*x2));
	}
      }
      return *this;
    }
    
    Vectorl& applyT(const GivensRotation& Q){
      assert(Q.i1<n); assert(Q.i2<n);
      if(Q.i1>Q.i2) return apply(Q);
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

     
  public: // KpointOp<k>
    
    template<int k>
    Vectorl& apply(const KpointOp<k>& Q){ 
      // it is assumed that Q.map is sorted
      auto it=lst.begin();
      int i=0;
      for(;i<k;i++){
	it=lower_bound(it,lst.end(),ivpair(Q.map(i),0),ivpairComparator);
	if(it==lst.end()) i=k;
	else if(it->first==Q.map(i)) break; 
      }
      if(i==k) return *this;
      SCALAR* vptr[k];
      SCALAR temp[k];
      it=lst.begin();
      for(int i=0; i<k; i++){
	int ix=Q.map(i);
	while(it!=lst.end()&&it->first<ix) it++;
	if(it==lst.end()||it->first>ix){
	  vptr[i]=&lst.insert(it,ivpair(ix,0))->second;
	  temp[i]=*vptr[i];
	}else{
	  vptr[i]=&it->second; 
	  temp[i]=it->second;
	}
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
    Vectorl& applyT(const KpointOp<k>& Q){ 
      // it is assumed that Q.map is sorted
      auto it=lst.begin();
      int i=0;
      for(;i<k;i++){
	it=lower_bound(it,lst.end(),ivpair(Q.map(i),0),ivpairComparator);
	if(it==lst.end()) i=k;
	else if(it->first==Q.map(i)) break; 
      }
      if(i==k) return *this;
      SCALAR* vptr[k];
      SCALAR temp[k];
      it=lst.begin();
      for(int i=0; i<k; i++){
	int ix=Q.map(i);
	while(it!=lst.end()&&it->first<ix) it++;
	if(it==lst.end()||it->first>ix){
	  vptr[i]=&lst.insert(it,ivpair(ix,0))->second;
	  temp[i]=*vptr[i];
	}else{
	  vptr[i]=&it->second; 
	  temp[i]=it->second;
	}
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
    
    Vectorl(double* numpyInDblArray, int numpyInSize): Vectorl(numpyInSize){
      for(int i=0; i<n; i++) {
	double v=numpyInDblArray[i];
	if(v!=0) lst.push_back(ivpair(i,v));
      }
    }
    
    void np(double** numpyOutArray1, int* numpyOutLen1){
      *numpyOutLen1=n;
      *numpyOutArray1=new double[n];
      std::fill(*numpyOutArray1,*numpyOutArray1+n,0);
      for(auto& p:lst) (*numpyOutArray1)[p.first]=p.second; 
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
    
    Vectorl& apply(const KpointOperator& Q){ // it is assumed that the indices in Q are sorted
      if(Q.k==0) return *this;
      vector<int> missing; auto it=lst.begin();
      for(int i=0; i<Q.k; i++) {int ix=Q.map(i);
	auto it2=find_if(it,lst.end(),[&Q,ix](ivpair& p)->bool{return p.first==ix;});
	if(it2==lst.end()) missing.push_back(ix); else it=it2++;
      }
      if(missing.size()==Q.k) return *this;
      for(int i=0; i<missing.size(); i++) insert(missing[i],0);
      int coi=0; int ix=Q.map(coi);
      for(auto it=lst.begin(); coi<Q.k; it++){
	assert(it!=lst.end());
	if(it->first==ix){Q.tempp[coi]=&it->second; Q.temp[coi]=it->second; ix=Q.map(++coi);}
      }
      for(int i=0; i<Q.k; i++){
	double s=0; 
	for(int j=0; j<Q.k; j++) s+=Q.q[j*Q.k+i]*Q.temp[j]; 
	*Q.tempp[i]=s;
      }
      return *this;
    }
    
    
    Vectorl& applyT(const KpointOperator& Q){ // it is assumed that the indices in Q are sorted
      if(Q.k==0) return *this;
      vector<int> missing; auto it=lst.begin();
      for(int i=0; i<Q.k; i++) {int ix=Q.map(i);
	auto it2=find_if(it,lst.end(),[&Q,ix](ivpair& p)->bool{return p.first==ix;});
	if(it2==lst.end()) missing.push_back(ix); else it=it2++;
      }
      if(missing.size()==Q.k) return *this;
      for(int i=0; i<missing.size(); i++) insert(missing[i],0);
      int coi=0; int ix=Q.map(coi);
      for(auto it=lst.begin(); coi<Q.k; it++){
	assert(it!=lst.end());
	if(it->first==ix){Q.tempp[coi]=&it->second; Q.temp[coi]=it->second; ix=Q.map(++coi);}
      }
      for(int i=0; i<Q.k; i++){
	double s=0; 
	for(int j=0; j<Q.k; j++) s+=Q.q[i*Q.k+j]*Q.temp[j]; 
	*Q.tempp[i]=s;
      }
      return *this;
    }
    */    
