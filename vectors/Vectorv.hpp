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


#ifndef _Vectorv
#define _Vectorv

#include "Cvector.hpp"
#include "SparseVector.hpp"
#include "GivensRotation.hpp"


extern default_random_engine randomNumberGenerator;


namespace Mondrian{
  
  class Vectorv: public SparseVector{
  public:  
    
    vector<ivpair> vec;
    bool sorted;
    static SCALAR dummyZero;


  private:

    mutex insrtmx;


  public:

    Vectorv(): SparseVector(0), sorted(false){}

    Vectorv(const int _n, const bool _sorted=false): SparseVector(_n), sorted(_sorted){}

    Vectorv(const initializer_list<SCALAR> list): Vectorv(list.size()){
      int i=0; for(SCALAR v:list) vec.push_back(ivpair(i++,v));}

    Vectorv(const int _n, const initializer_list<ivpair> list): Vectorv(n){
      for(const ivpair& v:list) vec.push_back(v);}

    ~Vectorv(){}


  public:

    Vectorv(const Vectorv& x): Vectorv(x.n){
      COPY_WARNING("Vectorv");
      vec=x.vec; 
      sorted=x.sorted;
    } 

    Vectorv(const Vectorv& x, const _NoWarn dummy): Vectorv(x.n){
      vec=x.vec; 
      sorted=x.sorted;
    } 

    Vectorv(Vectorv&& x): Vectorv(x.n){
      MOVE_WARNING("Vectorv");
      vec=std::move(x.vec);
      sorted=x.sorted;
    }

    Vectorv& operator=(const Vectorv& x){
      ASSIGN_WARNING("Vectorv");
      n=x.n; 
      vec=x.vec;
      sorted=x.sorted;
      return *this;
    }
    
    Vectorv& operator=(Vectorv&& x){
      MOVEASSIGN_WARNING("Vectorv");
      n=x.n; 
      vec=std::move(x.vec);
      sorted=x.sorted;
      return *this;
    }

    Vectorv copy(){
      Vectorv v(n); 
      v.vec=vec; 
      v.sorted=sorted;
      return v;
    }

    
  public: // conversions 

    Vectorv(const Cvector& x): SparseVector(x.n){
      for(int i=0; i<n; i++) if(x(i)!=0) vec.push_back(ivpair(i,x(i))); sorted=0;
    }

    Vectorv(const SparseVector& v): SparseVector(v.n){
      CONVERT_WARNING("SparseVector","Vectorv");
      v.for_each_filled([this](const int i, const SCALAR v){insert(i,v);});
    }
    
    operator Cvector(){
      CONVERT_WARNING("Vectorv","Cvector");
      Cvector v=Cvector::Zero(n);
      for(auto& p:vec) v(p.first)=p.second;
      return v;
    }
    

  public: // named constructors

    static Vectorv Zero(const int _n) {return Vectorv(_n,true);}

    static Vectorv Filled(const int _n, const SCALAR t){
      Vectorv v(_n,true); v.vec.resize(_n); 
      for(int i=0; i<_n; i++) v.vec[i]=ivpair(i,t); return v;
    }

    static Vectorv Uniform(const int _n){
      Vectorv v(_n,true); v.vec.resize(_n); 
      uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.vec[i]=ivpair(i,distr(randomNumberGenerator)); 
      return v;
    }

    static Vectorv Gaussian(const int _n){
      Vectorv v(_n,true); v.vec.resize(_n); 
      normal_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.vec[i]=ivpair(i,distr(randomNumberGenerator)); 
      return v;
    }

    static Vectorv Bernoulli(const int _n, const double p=0.5){
      Vectorv v(_n,true); 
      bernoulli_distribution distr(p);
      for(int i=0; i<_n; i++) if(distr(randomNumberGenerator)==1) v.vec.push_back(ivpair(i,1.0)); 
      return v;
    }


  public: // polymorphism

    virtual string classname() const {return "Vectorv";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer<Vectorv>();}


  public: // comparators

    bool operator==(const Vectorv& x) const { // check this
      const_cast<Vectorv&>(*this).sort(); 
      const_cast<Vectorv&>(x).sort(); 
      auto xit=x.vec.begin();
      for(auto& v:vec){
	if(v.second==0) continue;
	for(;xit!=x.vec.end() && xit->second==0; xit++){}
	if(xit->first!=v.first) return false;
	if(xit->second!=v.second) return false;
      }
      for(;xit!=x.vec.end(); xit++)
	if(xit->second!=0) return false; 
      return true;
    }

    
  public: // element access 
    
    SCALAR operator()(const int i) const {assert(i<n);
      auto it=vec.begin(); while(it!=vec.end() && it->first!=i){it++;}
      if(it==vec.end()) return 0; else return it->second;}

    SCALAR read(const int i) const {assert(i<n);
      auto it=vec.begin(); while(it!=vec.end() && it->first!=i){it++;}
      if(it==vec.end()) return 0; else return it->second;}


    SCALAR& operator()(const int i) {assert(i<n);
      auto it=vec.begin(); while(it!=vec.end() && it->first!=i){it++;}
      if(it==vec.end()) {vec.push_back(ivpair(i,0)); sorted=0; return vec.back().second;}
      return it->second;}

    void set(const int i, const SCALAR v) {assert(i<n);
      auto it=vec.begin(); while(it!=vec.end() && it->first!=i){it++;}
      if(it==vec.end()) {vec.push_back(ivpair(i,v)); sorted=0;}
      else it->second=v;
    }

    void set_msafe(const int i, const SCALAR v) {assert(i<n);
      lock_guard<mutex> lock(insrtmx);
      set(i,v);}

    SCALAR* ptr(const int i){assert(i<n);
      auto it=vec.begin(); 
      if(sorted) it=lower_bound(vec.begin(),vec.end(),ivpair(i,0),ivpairComparator);
      else while(it!=vec.end()&&it->first!=i) it++;
      if(it!=vec.end()&&it->first==i) return &it->second;
      sorted=0;
      return &vec.insert(vec.end(),ivpair(i,0))->second;
    }


    bool isFilled(const int i) const {assert(i<n);
      auto it=vec.begin(); while(it!=vec.end() && it->first!=i){it++;}
      return (it!=vec.end());}

    int nFilled() const {return vec.size();}
    

  public: // iterators

    void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      for(auto& p:vec) lambda(p.first,p.second);}
    
    void for_each_filled(std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(auto& p:vec) lambda(p.first,p.second);}

    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const TYPE, const TYPE)> accumulator, const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(auto& p:vec)
	t=accumulator(t,p.second);
      return t;
    }

    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const SCALAR&)> lambda, 
      std::function<TYPE(const TYPE, const TYPE)> accumulator=std::plus<TYPE>(), 
      const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(auto& p:vec)
	t=accumulator(t,lambda(p.second));
      return t;
    }

    INDEX find_best(std::function<bool(const SCALAR&, const SCALAR&)> selector) const{
      assert(vec.size()>0);
      INDEX besti=vec[0].first;
      SCALAR bestv=vec[0].second;
      for(auto& p:vec)
	if(selector(bestv,p.second)){besti=p.first; bestv=p.second;}
      return besti;
    }


  public: // sparse vector methods

    SCALAR* findptr(const int i){
      auto it=vec.begin(); 
      if(sorted){
	while(it!=vec.end() && it->first<i) {it++;} 
	if(it!=vec.end() && it->first==i) return &it->second;}
      else 
	for(;it!=vec.end();it++) 
	  if(it->first==i) return &it->second;
      return &dummyZero;
    }

    SCALAR* ptr_if_filled(const int i){
      auto it=vec.begin(); 
      if(sorted){
	while(it!=vec.end() && it->first<i) {it++;} 
	if(it!=vec.end() && it->first==i) return &it->second;}
      else 
	for(;it!=vec.end();it++) 
	  if(it->first==i) return &it->second;
      return nullptr;
    }

    const SCALAR* ptr_if_filled(const int i) const{
      auto it=vec.begin(); 
      if(sorted){
	while(it!=vec.end() && it->first<i) {it++;} 
	if(it!=vec.end() && it->first==i) return &it->second;}
      else 
	for(;it!=vec.end();it++) 
	  if(it->first==i) return &it->second;
      return nullptr;
    }

    vector<ivpair>::iterator find(const int i){
      auto it=vec.begin(); while(it!=vec.end() && it->first!=i) it++;
      return it;
    }

    void insert(const int i, const SCALAR value){
      assert(i<n); vec.push_back(ivpair(i,value)); sorted=0;}
    
    Vectorv& insert(const ivpair& p){
      vec.push_back(p); sorted=0; return *this;}
    
    void append(const int i, const SCALAR value){
      assert(i<n); vec.push_back(ivpair(i,value));}

    Vectorv& merge(vector<ivpair>& x){
      if(!sorted){
	for(auto p:x){
	  auto it=find(p.first);
	  if(it==vec.end()) vec.push_back(p);
	  else it->second=p.second;}
	return *this;
      }
      std::sort(x.begin(),x.end(),ivpairComparator);
      vector<ivpair> newv;
      newv.reserve(vec.size()+x.size());
      auto it=vec.begin();
      for(auto& p:x){
	auto lb=lower_bound(it,vec.end(),p,ivpairComparator);
	newv.insert(newv.end(),it,lb);
	if(lb==vec.end()) {newv.push_back(p); it=lb; continue;}
	if(lb->first>p.first) newv.push_back(p);
	else lb->second=p.second;
	it=lb;
      }
      newv.insert(newv.end(),it,vec.end());
      vec=std::move(newv);
      return *this;
    }
    
    void zero(const int i){assert(i<n);
      auto it=vec.begin(); 
      while(it!=vec.end() && it->first!=i){it++;}
      if(it!=vec.end()) it->second=0;}

    void sort() const {
      if(!sorted) std::sort(const_cast<vector<ivpair>&>(vec).begin(), const_cast<vector<ivpair>&>(vec).end(), ivpairComparator); 
      const_cast<Vectorv&>(*this).sorted=true;
    }
    
    void tidy(){
      sort(); auto v(vec); vec.clear(); int i=-1;
      for(auto& p:v) if(p.first!=i && p.second!=0) {vec.push_back(p); i=p.first;}
    }
  

  
  public: // scalar-valued operations 

    int nnz() const {const_cast<Vectorv&>(*this).tidy(); return vec.size();}
    
    SCALAR max() const{
      if(vec.size()==0) return 0;
      SCALAR max=vec.begin()->second;
      for(auto p:vec) if(p.second>max) max=p.second;
      return max;
    }
    
    SCALAR max_abs() const{
      if(vec.size()==0) return 0;
      SCALAR max=fabs(vec.begin()->second);
      for(auto p:vec) if(fabs(p.second)>max) max=fabs(p.second);
      return max;
    }
    
    int argmax() const{
      if(vec.size()==0) return 0;
      int best=vec.begin()->first; SCALAR max=vec.begin()->second;
      for(auto p:vec) if(p.second>max) {best=p.first; max=p.second;}
      return best;
    }
    
    int argmax_abs() const{
      if(vec.size()==0) return 0;
      int best=vec.begin()->first; SCALAR max=fabs(vec.begin()->second);
      for(auto p:vec) if(fabs(p.second)>max) {best=p.first; max=fabs(p.second);}
      return best;
    }
    
    SCALAR min() const{
      if(vec.size()==0) return 0;
      SCALAR min=vec.begin()->second;
      for(auto p:vec) if(p.second<min) min=p.second;
      return min;
    }
    
    SCALAR min_abs() const{
      if(vec.size()==0) return 0;
      SCALAR min=fabs(vec.begin()->second);
      for(auto p:vec) if(fabs(p.second)<min) min=fabs(p.second);
      return min;
    }
    
    int argmin() const{
      if(vec.size()==0) return 0;
      int best=vec.begin()->first; SCALAR min=vec.begin()->second;
      for(auto p:vec) if(p.second<min) {best=p.first; min=p.second;}
      return best;
    }
    
    int argmin_abs() const{
      if(vec.size()==0) return 0;
      int best=vec.begin()->first; SCALAR min=fabs(vec.begin()->second);
      for(auto p:vec) if(fabs(p.second)<min) {best=p.first; min=fabs(p.second);}
      return best;
    }
    
    SCALAR sum() const{
      SCALAR t=0; for(auto& p:vec) t+=p.second; return t;}

    SCALAR norm1() const{
      SCALAR t=0; for(auto& p:vec) t+=fabs(p.second); return t;}
    
    SCALAR norm2() const{
      SCALAR t=0; for(auto& p:vec) t+=p.second*p.second; return t;}
    
    SCALAR diff2(const Vectorv& x) const {
      assert(x.n==n); SCALAR t=0; 
      const_cast<Vectorv&>(*this).sort(); 
      const_cast<Vectorv&>(x).sort(); 
      auto it=vec.begin();
      for(auto& xv:x.vec){
	while(it!=vec.end() && it->first<xv.first) {t+=(it->second)*(it->second); it++;}
	if(it!=vec.end() && it->first==xv.first) {t+=(it->second-xv.second)*(it->second-xv.second);it++;}
	else {t+=xv.second*xv.second;}
      }
      return t; 
    }

    SCALAR diff2(const Cvector& x) const {
      assert(x.n==n); SCALAR t=0; 
      const_cast<Vectorv&>(*this).tidy(); 
      int ix=0;
      for(auto& p:vec){
	for(int i=ix; i<p.first; i++) t+=x.array[i]*x.array[i];
	t+=(p.second-x.array[p.first])*(p.second-x.array[p.first]);
	ix=p.first+1;
      }
      for(int i=ix; i<n; i++) t+=x.array[i]*x.array[i];
      return t; 
    }

    
  public: // scalar valued arithmetic

    SCALAR dot(const Vectorv& x) const{
      SCALAR t=0; 
      const_cast<Vectorv&>(*this).sort(); 
      const_cast<Vectorv&>(x).sort(); 
      auto xit=x.vec.begin();
      for(auto& v:vec){
	int i=v.first;
	//while(xit!=x.vec.end() && xit->first<i) xit++;
	xit=lower_bound(xit,x.vec.end(),ivpair(i,0),ivpairComparator);
	if(xit==x.vec.end()) break;
	if(xit->first==i) {t+=v.second*xit->second; xit++;}
       }
      return t;
    }
    
    
  public: // in-place operations 
    
    Vectorv& operator*=(const SCALAR c){
      for(auto& p:vec) p.second*=c; return *this;}
    Vectorv& operator/=(const SCALAR c){
      for(auto& p:vec) p.second/=c; return *this;}
    
    template<class VECTOR>
    Vectorv& operator*=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:vec) p.second*=x(p.first); return *this;}
    template<class VECTOR>
    Vectorv& operator/=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:vec) p.second/=x(p.first); return *this;}
    
    Vectorv& add(const Vectorv& x, const SCALAR c=1){
      sort();
      const_cast<Vectorv&>(x).sort(); 
      vector<ivpair> newbies;
      auto it=vec.begin();
      for(auto& xv:x.vec){
	while(it!=vec.end() && it->first<xv.first) it++;
	if(it!=vec.end() && it->first==xv.first) it->second+=c*xv.second;
	else newbies.push_back(ivpair(xv.first,c*xv.second));
      }
      if(newbies.size()!=0) {vec.insert(vec.end(),newbies.begin(),newbies.end()); sorted=0;}
      return *this; 
    }

    Vectorv& operator+=(const Vectorv& x){return add(x,1);}
    Vectorv& operator-=(const Vectorv& x){return add(x,-1);}

    Vectorv& operator*=(const Vectorv& x){
      sort();
      const_cast<Vectorv&>(x).sort(); 
      auto it=x.vec.begin();
      for(auto& p:vec){
	while(it!=x.vec.end() && it->first<p.first) it++;
	if(it!=x.vec.end() && it->first==p.first) p.second*=it->second;
	else p.second=0;
      }
      return *this; 
    }
    
    Vectorv& operator/=(const Vectorv& x){
      sort();
      const_cast<Vectorv&>(x).sort(); 
      auto it=x.vec.begin();
      for(auto& p:vec){
	while(it!=x.vec.end() && it->first<p.first) it++;
	if(it!=x.vec.end() && it->first==p.first) p.second/=it->second;
	else p.second=1.0/0.0;
      }
      return *this; 
    }
    

public: // Givens rotations 
    
    Vectorv& apply(const GivensRotation& Q){
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
    
    Vectorv& applyT(const GivensRotation& Q){
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

     
  public: // KpointOp<k>
    
    template<int k>
    Vectorv& apply(const KpointOp<k>& Q){ // improve!
      // it is assumed that Q.map is sorted
      SCALAR* vptr[k];
      SCALAR temp[k];
      for(int i=0; i<k; i++) vec.push_back(ivpair(-3,0));
      int nfound=0;
      if(sorted){
	auto it=vec.begin();
	for(int i=0; i<k; i++){
	  int ix=Q.map(i);
	  //it=lower_bound(it,vec.end(),ivpair(ix,0),ivpairComparator); // the -3 entries mess this up 
	  while(it!=vec.end()&&it->first<ix) it++;
	  if(it==vec.end()||it->first>ix) vptr[i]=nullptr;
	  else{vptr[i]=&it->second; temp[i]=it->second; nfound++;}
	}
      }else{
	for(int i=0; i<k; i++){
	  int ix=Q.map(i);
	  auto it=vec.begin();
     	  while(it!=vec.end()&&it->first!=ix) it++;
	  if(it==vec.end()) vptr[i]=nullptr;
	  else{vptr[i]=&it->second; temp[i]=it->second; nfound++;}
	}
      }
      if(nfound==0){
	for(int i=0; i<k; i++) 
	  vec.pop_back();
	return *this;
      }
      int c=vec.size()-k;
      for(int i=0; i<k; i++)
	if(vptr[i]==nullptr){
	  vec[c].first=Q.map(i);
	  vptr[i]=&vec[c].second;
	  temp[i]=0;
	  c++;
	}
      for(int i=0; i<k; i++){
	SCALAR s=0; 
	for(int j=0; j<k; j++) 
	  s+=Q.q[j*k+i]*temp[j]; 
	*vptr[i]=s;
      }
      for(int i=0; i<nfound; i++) 
	vec.pop_back();
      if(nfound<k) sorted=false;
      return *this;
    }
     
    template<int k>
    Vectorv& applyT(const KpointOp<k>& Q){ 
      // it is assumed that Qmap is sorted
      SCALAR* vptr[k];
      SCALAR temp[k];
      for(int i=0; i<k; i++) vec.push_back(ivpair(-1,0));
      int nfound=0;
      if(sorted){
	auto it=vec.begin();
	for(int i=0; i<k; i++){
	  int ix=Q.map(i);
	  //it=lower_bound(it,vec.end(),ivpair(ix,0),ivpairComparator);
	  while(it!=vec.end()&&it->first<ix) it++;
	  if(it==vec.end()||it->first>ix) vptr[i]=nullptr;
	  else{vptr[i]=&it->second; temp[i]=it->second; nfound++;}
	}
      }else{
	for(int i=0; i<k; i++){
	  int ix=Q.map(i);
	  auto it=vec.begin();
     	  while(it!=vec.end()&&it->first!=ix) it++;
	  if(it==vec.end()) vptr[i]=nullptr;
	  else{vptr[i]=&it->second; temp[i]=it->second; nfound++;}
	}
      }
      if(nfound==0){
	for(int i=0; i<k; i++) 
	  vec.pop_back();
	return *this;
      }
      int c=vec.size()-k;
      for(int i=0; i<k; i++)
	if(vptr[i]==nullptr){
	  vec[c].first=Q.map(i);
	  vptr[i]=&vec[c].second;
	  temp[i]=0;
	  c++;
	}
      for(int i=0; i<k; i++){
	SCALAR s=0; 
	for(int j=0; j<k; j++) 
	  s+=Q.q[i*k+j]*temp[j]; 
	*vptr[i]=s;
      }
      for(int i=0; i<nfound; i++) 
	vec.pop_back();
      if(nfound<k) sorted=false;
      return *this;
    }

    
  public: // I/O 
    
    //string str(const Sparse dummy) const;
    //string str(const Dense dummy) const {return Vector::str(Dense());}

    string str() const {
      ostringstream result;
      for(auto& p:vec) cout<<"("<<p.first<<","<<p.second<<")"<<endl;
      //result<<"(";
      //cout<<vec.size()<<endl;
      //for(int i=0; i<vec.size()-1; i++) cout<<i<<endl;
      //result<<vec[i].first<<",";
      //result<<vec[vec.size()-1].first<<")";
      return result.str();
    }
    
  

  //bool equals(const Vector& x) const {return (*this)==(static_cast<const Vectorv&>(x));} // needs work!

  
public: // Python interface

  Vectorv(double* numpyInDblArray, int numpyInSize): Vectorv(numpyInSize){
    for(int i=0; i<n; i++) {
      double v=numpyInDblArray[i];
      if(v!=0) vec.push_back(ivpair(i,v));
    }
  }

  void np(double** numpyOutArray1, int* numpyOutLen1){
    *numpyOutLen1=n;
    *numpyOutArray1=new double[n];
    std::fill(*numpyOutArray1,*numpyOutArray1+n,0);
    for(auto& p:vec) (*numpyOutArray1)[p.first]=p.second; 
  }
  

};
  
}




#endif



// ----------------------------------------------------------------------------------------------------------


    //Vectorv* clone() const {return new Vectorv(*this);}
    //bool isDense() const {return false;}
   /*
  public: // KpointOperators
    
    Vectorv& apply(const KpointOperator& Q){ // it is assumed that the indices in Q are sorted
      if(Q.k==0) return *this;
      sort(); vector<int> missing; auto it=vec.begin(); // improve! does not take sorted nature into account
      for(int i=0; i<Q.k; i++) {int ix=Q.map(i);
	auto it2=find_if(it,vec.end(),[&Q,ix](ivpair& p)->bool{return p.first==ix;});
	if(it2==vec.end()) missing.push_back(ix);
	else it=it2++;
      }
      if(missing.size()==Q.k) return *this;
      for(int i=0; i<missing.size(); i++) insert(missing[i],0);
      int coi=0; int ix=Q.map(coi);
      for(auto it=vec.begin(); coi<Q.k; it++){
	assert(it!=vec.end());
	if(it->first==ix){Q.tempp[coi]=&it->second; Q.temp[coi]=it->second; ix=Q.map(++coi);}
      }
      for(int i=0; i<Q.k; i++){
	double s=0; 
	for(int j=0; j<Q.k; j++) s+=Q.q[j*Q.k+i]*Q.temp[j]; 
	*Q.tempp[i]=s;
      }
      return *this;
    }
    
    Vectorv& applyT(const KpointOperator& Q){ // it is assumed that the indices in Q are sorted
      if(Q.k==0) return *this; 
      sort(); vector<int> missing; auto it=vec.begin();
      for(int i=0; i<Q.k; i++) {int ix=Q.map(i);
	auto it2=find_if(it,vec.end(),[&Q,ix](ivpair& p)->bool{return p.first==ix;});
	if(it2==vec.end()) missing.push_back(ix);
	else it=it2++;
      }
      if(missing.size()==Q.k) return *this;
      for(int i=0; i<missing.size(); i++) insert(missing[i],0);
      int coi=0; int ix=Q.map(coi);
      for(auto it=vec.begin(); coi<Q.k; it++){
	assert(it!=vec.end());
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
