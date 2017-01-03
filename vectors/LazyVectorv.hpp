#ifndef _LazyVectorv
#define _LazyVectorv


#include "Vectorv.hpp"


extern default_random_engine randomNumberGenerator;


namespace Mondrian{
  
  class LazyVectorv: public SparseVector{
  public:  
    
    vector<ivpair> vec;
    vector<ivpair> chg;

    //bool sorted=true;
    int chg_max=100;
    static SCALAR dummyZero;


  public:

    LazyVectorv(): SparseVector(0){}

    LazyVectorv(const int _n): SparseVector(_n){chg.reserve(chg_max);}

    LazyVectorv(const initializer_list<SCALAR> list): LazyVectorv(list.size()){
      int i=0; for(SCALAR v:list) vec.push_back(ivpair(i++,v)); sort();}

    LazyVectorv(const int _n, const initializer_list<ivpair> list): LazyVectorv(n){
      for(const ivpair& v:list) vec.push_back(v); sort();}

    ~LazyVectorv(){}


  public:

    LazyVectorv(const LazyVectorv& x): LazyVectorv(x.n){
      COPY_WARNING("LazyVectorv");
      vec=x.vec; 
      chg=x.chg;
    } 

    LazyVectorv(const LazyVectorv& x, const _NoWarn dummy): LazyVectorv(x.n){
      vec=x.vec;
      chg=x.chg;
    } 

    LazyVectorv(LazyVectorv&& x): LazyVectorv(x.n){
      MOVE_WARNING("LazyVectorv");
      vec=std::move(x.vec);
      chg=std::move(x.chg);
    }

    LazyVectorv& operator=(const LazyVectorv& x){
      ASSIGN_WARNING("LazyVectorv");
      n=x.n; 
      vec=x.vec;
      chg=x.chg;
      return *this;
    }
    
    LazyVectorv& operator=(LazyVectorv&& x){
      MOVEASSIGN_WARNING("LazyVectorv");
      n=x.n; 
      vec=std::move(x.vec);
      chg=std::move(x.chg);
      return *this;
    }

    LazyVectorv copy(){
      LazyVectorv v(n); 
      v.vec=vec; 
      v.chg=chg;
      return v;
    }

    
public: // conversions 
  
    LazyVectorv(const Cvector& x): LazyVectorv(x.n){
      for(int i=0; i<n; i++) if(x(i)!=0) vec.push_back(ivpair(i,x(i)));
    }

    LazyVectorv(const SparseVector& v): SparseVector(v.n){
      CONVERT_WARNING("SparseVector","LazyVectorv");
      v.for_each_filled([this](const int i, const SCALAR v){insert(i,v);});
    }
  
    operator Cvector(){
      CONVERT_WARNING("LazyVectorv","Cvector");
      Cvector v=Cvector::Zero(n);
      for(auto& p:vec) v(p.first)=p.second;
      return v;
    }


  public: // named constructors

    static LazyVectorv Zero(const int _n) {return LazyVectorv(_n);}

    static LazyVectorv Filled(const int _n, const SCALAR t){
      LazyVectorv v(_n); v.vec.resize(_n); for(int i=0; i<_n; i++) v.vec[i]=ivpair(i,t); return v;}

    static LazyVectorv Uniform(const int _n){
      LazyVectorv v(_n); v.vec.resize(_n); 
      uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.vec[i]=ivpair(i,distr(randomNumberGenerator)); 
      return v;}

    static LazyVectorv Gaussian(const int _n){
      LazyVectorv v(_n); v.vec.resize(_n); 
      normal_distribution<SCALAR> distr;
      for(int i=0; i<_n; i++) v.vec[i]=ivpair(i,distr(randomNumberGenerator)); 
      return v;}

    static LazyVectorv Bernoulli(const int _n, const double p=0.5){
      LazyVectorv v(_n); 
      bernoulli_distribution distr(p);
      for(int i=0; i<_n; i++) if(distr(randomNumberGenerator)==1) v.vec.push_back(ivpair(i,1.0)); 
      return v;}


  public: // polymorphism

    virtual string classname() const {return "LazyVectorv";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer<LazyVectorv>();}


  public: // comparators

    bool operator==(const LazyVectorv& x) const {
      const_cast<LazyVectorv&>(*this).flush(); 
      const_cast<LazyVectorv&>(x).flush();
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
      auto p=find_ptr(i);
      if(p!=nullptr) return p->second;
      else return 0;
    }

    SCALAR read(const int i) const {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) return p->second;
      else return 0;
    }

    void set(const int i, const SCALAR v) {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) {p->second=v; return;}
      chg.push_back(ivpair(i,v));
      return;
    }

    void set_if_filled(const int i, const SCALAR v) {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) p->second=v;
      return;
    }

    SCALAR& operator()(const int i) {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) return p->second;
      chg.push_back(ivpair(i,0));
      return chg.back().second;
    }

    bool isFilled(const int i) const {assert(i<n);
      return (find_ptr(i)!=nullptr);
    }
    
    int nFilled() const{
      return vec.size()+chg.size();
    }
    

  public: // iterators

    void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      for(auto& p:vec) lambda(p.first,p.second);
      for(auto& p:chg) lambda(p.first,p.second);
    }
    
    void for_each_filled(std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(auto& p:vec) lambda(p.first,p.second);
      for(auto& p:chg) lambda(p.first,p.second);
    }


  public: // sparse vector methods

    SCALAR* findptr(const int i){
      auto it=find_ptr(i);
      if(it!=nullptr) return &it->second;
      return &dummyZero;
    }

    const ivpair* find_ptr(const int i) const{
      auto it=chg.begin(); 
      for(;it!=chg.end();it++) 
	if(it->first==i) return &*it;
      it=lower_bound(vec.begin(),vec.end(),ivpair(i,0),ivpairComparator);
      if(it!=vec.end()&&it->first==i) return &*it;
      else return nullptr;
    }

    ivpair* find_ptr(const int i){
      auto it=chg.begin(); 
      for(;it!=chg.end();it++) 
	if(it->first==i) return &*it;
      it=lower_bound(vec.begin(),vec.end(),ivpair(i,0),ivpairComparator);
      if(it!=vec.end()&&it->first==i) return &*it;
      else return nullptr;
    }

    ivpair* get_ptr(const int i){
      auto it=chg.begin(); 
      for(;it!=chg.end();it++) 
	if(it->first==i) return &*it;
      it=lower_bound(vec.begin(),vec.end(),ivpair(i,0),ivpairComparator);
      if(it!=vec.end()&&it->first==i) return &*it;
      chg.push_back(ivpair(i,0));
      return &chg.back();
    }

    pair<ivpair*,ivpair*> get_ptr_unless_all_zero(int i1, int i2){
      assert(i1<n); assert(i2<n);
      bool swapp=false; if(i1>i2) {int t=i1; i1=i2; i2=t; swapp=true;}
      auto it1=chg.begin(); 
      auto it2=chg.begin(); 
      for(;it1!=chg.end();it1++) 
	if((it1)->first==i1) break;
      for(;it2!=chg.end();it2++) 
	if((it2)->first==i2) break;
      if(it1==chg.end()){
	it1=lower_bound(vec.begin(),vec.end(),ivpair(i1,0),ivpairComparator);
	if(it2==chg.end()) it2=lower_bound(it1,vec.end(),ivpair(i2,0),ivpairComparator);	
      }else{
	if(it2==chg.end()) it2=lower_bound(vec.begin(),vec.end(),ivpair(i2,0),ivpairComparator);	
      }
      ivpair* p1;
      if(it1!=vec.end()&&it1->first==i1) p1=&*it1;
      else{
	if(it2==vec.end()||it2->first!=i2) return pair<ivpair*,ivpair*>(nullptr,nullptr);
	else chg.push_back(ivpair(i1,0)); p1=&chg.back();
      }
      ivpair* p2;
      if(it2!=vec.end()&&it2->first==i2) p2=&*it2;
      else {chg.push_back(ivpair(i2,0)); p2=&chg.back();}
      if(swapp) return pair<ivpair*,ivpair*>(p2,p1);
      else return pair<ivpair*,ivpair*>(p1,p2);
    }

    void insert(const int i, const SCALAR value){
      chg.push_back(ivpair(i,value));}
    
    LazyVectorv& insert(const ivpair& p){
      chg.push_back(p); 
      return *this;
    }
    
    void append(const int i, const SCALAR value){assert(i<n); 
      chg.push_back(ivpair(i,value));}

    
    LazyVectorv& merge(vector<ivpair>& x){
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
    
    
    void sort() const {}
  
    
    void tidy(){}


  public: // lazy vector methods

    void flush() const{
      const_cast<LazyVectorv&>(*this).flush();
    }

    void flush(){
      std::sort(chg.begin(),chg.end(),ivpairComparator);
      vector<ivpair> newv;
      newv.reserve(vec.size()+chg.size());
      auto it=vec.begin();
      for(auto& p:chg){
	auto lb=lower_bound(it,vec.end(),p,ivpairComparator);
	newv.insert(newv.end(),it,lb);
	if(lb==vec.end()) {newv.push_back(p); it=lb; continue;}
	if(lb->first>p.first) newv.push_back(p);
	else lb->second=p.second;
	it=lb;
      }
      newv.insert(newv.end(),it,vec.end());
      vec=std::move(newv);
      chg.clear();
    }

    void setCache(const int s){
      chg_max=s; chg.reserve(s);}

  public: // scalar-valued operations 
    
    int nnz() const {
      const_cast<LazyVectorv&>(*this).tidy(); 
      return vec.size();
    }
   
    template<class COMPARATOR>
    const ivpair* find_best(const COMPARATOR& comp) const{ // improve?
      const ivpair* best;
      if(chg.size()==0) best=&*vec.begin();
      else best=&*chg.begin();
      for(auto& p:chg) if(comp(p,*best)) best=&p;
      for(auto& p:vec) if(comp(p,*best)) best=&p;
      return best;
    }
    
    SCALAR max() const{
      if(nFilled()==0) return 0;
      return (find_best(ivlarger))->second;
    }
    
    SCALAR max_abs() const{
      if(nFilled()==0) return 0;
      return (find_best(ivlarger_abs))->second;
    }

    INDEX argmax() const{
      if(nFilled()==0) return 0;
      return (find_best(ivlarger))->first;
    }
    
    INDEX argmax_abs() const{
      if(nFilled()==0) return 0;
      return (find_best(ivlarger_abs))->first;
    }

    SCALAR min() const{
      if(nFilled()==0) return 0;
      return (find_best(ivsmaller))->second;
    }
    
    SCALAR min_abs() const{
      if(nFilled()==0) return 0;
      return (find_best(ivsmaller_abs))->second;
    }

    INDEX argmin() const{
      if(nFilled()==0) return 0;
      return (find_best(ivsmaller))->first;
    }
    
    INDEX argmin_abs() const{
      if(nFilled()==0) return 0;
      return (find_best(ivsmaller_abs))->first;
    }
    
    SCALAR sum() const{
      SCALAR t=0; 
      for(auto& p:chg) t+=p.second; 
      for(auto& p:vec) t+=p.second; 
      return t;
    }

    SCALAR norm1() const{
      SCALAR t=0; 
      for(auto& p:chg) t+=fabs(p.second); 
      for(auto& p:vec) t+=fabs(p.second); 
      return t;
    }
    
    SCALAR norm2() const{
      SCALAR t=0; 
      for(auto& p:chg) t+=p.second*p.second; 
      for(auto& p:vec) t+=p.second*p.second; 
      return t;
    }
    
    SCALAR diff2(const LazyVectorv& x) const {
      cout<<"Unimplemented"<<endl;
      return 0;
    }
    
    SCALAR diff2(const Cvector& x) const {
      assert(x.n==n); flush(); SCALAR t=0; 
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

    SCALAR dot(const LazyVectorv& x) const{
      SCALAR t=0; 
      const_cast<LazyVectorv&>(*this).sort(); 
      const_cast<LazyVectorv&>(x).sort(); 
      auto xit=x.vec.begin();
      for(auto& v:vec){
	int i=v.first;
	while(xit!=x.vec.end() && xit->first<i) xit++;
	if(xit==x.vec.end()) break;
	if(xit->first==i) t+=v.second*xit->second;
      }
      return t;
    }
    
    
  public: // in-place operations 
    
    LazyVectorv& operator*=(const SCALAR c){
      for(auto& p:chg) p.second*=c; 
      for(auto& p:vec) p.second*=c; 
      return *this;
    }
    
    LazyVectorv& operator/=(const SCALAR c){
      for(auto& p:chg) p.second/=c; 
      for(auto& p:vec) p.second/=c; 
      return *this;
    }
    
    template<class VECTOR>
    LazyVectorv& operator*=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:chg) p.second*=x(p.first); 
      for(auto& p:vec) p.second*=x(p.first); 
      return *this;
    }
    
    template<class VECTOR>
    LazyVectorv& operator/=(const VECTOR& x){assert(x.n==n); 
      for(auto& p:chg) p.second/=x(p.first); 
      for(auto& p:vec) p.second/=x(p.first); 
      return *this;
    }
    
    LazyVectorv& add(const LazyVectorv& x, const SCALAR c=1){
      sort();
      const_cast<LazyVectorv&>(x).sort(); 
      vector<ivpair> newbies;
      auto it=vec.begin();
      for(auto& xv:x.vec){
	while(it!=vec.end() && it->first<xv.first) it++;
	if(it!=vec.end() && it->first==xv.first) it->second+=c*xv.second;
	else newbies.push_back(ivpair(xv.first,c*xv.second));
      }
      if(newbies.size()!=0) {vec.insert(vec.end(),newbies.begin(),newbies.end());}
      return *this; 
    }

    LazyVectorv& operator+=(const LazyVectorv& x){return add(x,1);}
    LazyVectorv& operator-=(const LazyVectorv& x){return add(x,-1);}

    LazyVectorv& operator*=(const LazyVectorv& x){
      const_cast<LazyVectorv&>(x).sort(); 
      auto it=x.vec.begin();
      for(auto& p:vec){
	while(it!=x.vec.end() && it->first<p.first) it++;
	if(it!=x.vec.end() && it->first==p.first) p.second*=it->second;
	else p.second=0;
      }
      return *this; 
    }
    
    LazyVectorv& operator/=(const LazyVectorv& x){
      sort();
      const_cast<LazyVectorv&>(x).sort(); 
      auto it=x.vec.begin();
      for(auto& p:vec){
	while(it!=x.vec.end() && it->first<p.first) it++;
	if(it!=x.vec.end() && it->first==p.first) p.second/=it->second;
	else p.second=1.0/0.0;
      }
      return *this; 
    }
    

  public: // Givens rotations 
    
    LazyVectorv& apply(const GivensRotation& Q){
      auto pp=get_ptr_unless_all_zero(Q.i1,Q.i2);
      if(pp.first==nullptr) return *this;
      SCALAR* p1=&pp.first->second; 
      SCALAR* p2=&pp.second->second;
      SCALAR x1=*p1;
      *p1=Q.cos*x1-Q.sin*(*p2);
      *p2=Q.sin*x1+Q.cos*(*p2);
      if(chg.size()>chg_max) flush();
      return *this;
    }
    
    LazyVectorv& applyT(const GivensRotation& Q){
      auto pp=get_ptr_unless_all_zero(Q.i1,Q.i2);
      if(pp.first==nullptr) return *this;
      SCALAR* p1=&pp.first->second; 
      SCALAR* p2=&pp.second->second;
      SCALAR x1=*p1;
      *p1=Q.cos*x1+Q.sin*(*p2);
      *p2=-Q.sin*x1+Q.cos*(*p2);
      if(chg.size()>chg_max) flush();
      return *this;
    }
    


     
  public: // I/O 
    
    //string str(const Sparse dummy) const;
    //string str(const Dense dummy) const {return Vector::str(Dense());}

    string str() const {
      ostringstream result;
      for(auto& p:vec) result<<"("<<p.first<<","<<p.second<<")"<<endl;
      result<<" - - - - "<<endl;
      for(auto& p:chg) result<<"("<<p.first<<","<<p.second<<")"<<endl;
      result<<endl;
      //result<<"(";
      //cout<<vec.size()<<endl;
      //for(int i=0; i<vec.size()-1; i++) cout<<i<<endl;
      //result<<vec[i].first<<",";
      //result<<vec[vec.size()-1].first<<")";
      return result.str();
    }
    
  

  //bool equals(const Vector& x) const {return (*this)==(static_cast<const LazyVectorv&>(x));} // needs work!

  
public: // Python interface

  LazyVectorv(double* numpyInDblArray, int numpyInSize): LazyVectorv(numpyInSize){
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


    //LazyVectorv* clone() const {return new Vectorv(*this);}
    //bool isDense() const {return false;}
   /*
  public: // KpointOperators
    
    LazyVectorv& apply(const KpointOperator& Q){ // it is assumed that the indices in Q are sorted
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
    
    LazyVectorv& applyT(const KpointOperator& Q){ // it is assumed that the indices in Q are sorted
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
