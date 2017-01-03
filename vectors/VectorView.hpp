#ifndef _VectorView
#define _VectorView

#include "IndexMap.hpp"


namespace Mondrian{

template<class VECTOR>
class VectorView: public Vector{
public:

  VECTOR& v;
  IndexMap map;


public:

  VectorView(VECTOR& _v, const IndexMap& _map): 
    Vector(_map.nsource), v(_v), map(_map){}

  VectorView(VECTOR& _v, IndexMap&& _map): 
    Vector(_map.nsource), v(_v), map(std::move(_map)){}


public: // copying and assignment

  VectorView(const VectorView& x): Vector(x.n), v(x.v), map(x.map){}

  VectorView(VectorView&& x): Vector(x.n), v(x.v), map(std::move(x.map)){}
  
  // VectorView& operator=...

  ~VectorView(){}


public: // polymorphism
  
  static string classname() {return "VectorView<"+VECTOR::classname()+">";}
  virtual VectorSpecializerBase* specializer() const{
    return new VectorSpecializer<VectorView<VECTOR> >();}


public: // attributes

  bool isSparse() const {return v.isSparse();}


public: // conversions

  // operator Cvector() const; inherited from Vector


public: // element access

  SCALAR read(const int i) const {assert(i<n); return v(map(i));} 
  SCALAR operator()(const int i) const {assert(i<n); return v(map(i));} 

  void set(const int i, const SCALAR x){assert(i<n); v.set(map(i),x);}
  SCALAR& operator()(const int i) {assert(i<n); return v(map(i));} 

  bool isFilled(const int i) const {return v.isFilled(map(i));}
  int nFilled() const {int t=0; for(int i=0; i<n; i++) t+=v.isFilled(map(i)); return t;} 


public: //assignments

  template<class VECTOR2>
  VectorView<VECTOR>& operator=(const VECTOR2& x){
    assert(n==x.n);
    for(int i=0; i<n; i++) v(map(i))=x(i);
    return *this;
  }

  
public: // iterators

  void for_each(std::function<void(const INDEX, SCALAR&)> lambda){
    for(int i=0; i<n; i++) lambda(i,v(map(i)));}
  void for_each(std::function<void(const INDEX, SCALAR)> lambda) const{
    for(int i=0; i<n; i++) lambda(i,v.read(map(i)));}

  void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
    for(int i=0; i<n; i++) if(v.isFilled(map(i))) lambda(i,v(map(i)));}
  void for_each_filled(std::function<void(const INDEX, SCALAR)> lambda) const{
    for(int i=0; i<n; i++) if(v.isFilled(map(i))) lambda(i,v.read(map(i)));}


public: // views

  VectorView<VectorView<VECTOR> > operator()(const IndexMap& phi) {
    return VectorView<VectorView<VECTOR> >(*this,phi);}
  // should this be done with map composition?


public: // vector valued arithmetic

  VECTOR mult(const SCALAR c) const{
    VECTOR r(n); for(int i=0; i<n; i++) v(i)=c*v(map(i)); return r;}

  template<class VECTOR2>
  VECTOR plus(const VECTOR2& w) const{
    VECTOR r(n); for(int i=0; i<n; i++) v(i)=w(i)+v(map(i)); return r;}
  template<class VECTOR2>
  VECTOR minus(const VECTOR2& w) const{
    VECTOR r(n); for(int i=0; i<n; i++) v(i)=w(i)+v(map(i)); return r;}

  VECTOR operator*(const SCALAR c) const {return this->mult(c);}
  template<class VECTOR2> VECTOR operator+(const VECTOR2& w) const {return this->plus(w);}
  template<class VECTOR2> VECTOR operator-(const VECTOR2& w) const {return this->minus(w);}


public: // scalar valued methods
    
  int nnz() const{int t=0; for(int i=0; i<n; i++) if (v(map(i))!=0) t++; return t;}
  SCALAR sum() const{SCALAR t=0; for(int i=0; i<n; i++) t+=v(map(i)); return t;}
  SCALAR norm1() const {SCALAR t=0; for(int i=0; i<n; i++) t+=fabs(v(map(i))); return t;}
  SCALAR norm2() const {SCALAR t=0; for(int i=0; i<n; i++) t+=v(map(i))*v(map(i)); return t;}

  template<class VECTOR2>
  SCALAR diff2(const VECTOR2& x) const{
    assert(x.n==n); SCALAR t=0; 
    for(int i=0; i<n; i++) t+=(v(map(i))-x(i))*(v(map(i))-x(i)); 
    return t;}
  
  SCALAR max() const{
    assert(n>0); SCALAR t=v(map(0)); 
    for(int i=1; i<n; i++) if(v(map(i))>t) t=v(map(i)); 
    return t;}
  SCALAR max_abs() const{
    assert(n>0); SCALAR t=fabs(v(map(0))); 
    for(int i=1; i<n; i++) if(fabs(v(map(i)))>t) t=v(map(i)); 
    return t;}
  INDEX argmax() const{ 
    assert(n>0); int best=0; SCALAR max=v(map(0));  
    for(int i=0; i<n; i++) if(v(map(i))>max) {best=i; max=v(map(i));}
    return best;}
  INDEX argmax_abs() const{
    if(n==0) return 0; int best=0; SCALAR max=fabs(v(map(0)));  
    for(int i=0; i<n; i++) if(fabs(v(map(i)))>max) {best=i; max=fabs(v(map(i)));}
    return best;}
  
  SCALAR min() const{
    assert(n>0); SCALAR t=v(map(0)); 
    for(int i=1; i<n; i++) if(v(map(i))<t) t=v(map(i)); 
    return t;}
  SCALAR min_abs() const{
    assert(n>0); SCALAR t=fabs(v(map(0))); 
    for(int i=1; i<n; i++) if(fabs(v(map(i)))<t) t=v(map(i)); 
    return t;}
  INDEX argmin() const{ 
    assert(n>0); int best=0; SCALAR max=v(map(0));  
    for(int i=0; i<n; i++) if(v(map(i))<max) {best=i; max=v(map(i));}
    return best;}
  INDEX argmin_abs() const{
    if(n==0) return 0; int best=0; SCALAR max=fabs(v(map(0)));  
    for(int i=0; i<n; i++) if(fabs(v(map(i)))<max) {best=i; max=fabs(v(map(i)));}
    return best;}
  
  
 

public:


};


  /* Where to put this? 
template<class VECTOR>
VectorView<VECTOR>::operator Cvector() const{
  Cvector w(n);
  for(int i=0; i<n; i++) w(i)=v(map(i))l
    return w;
}
  */ 

}
#endif
