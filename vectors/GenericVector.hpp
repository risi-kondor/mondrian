#ifndef _GenericVector
#define _GenericVector

#include "Vector.hpp"

namespace Mondrian{


class GenericVector: public Vector{
public:

  GenericVector(){}
  GenericVector(const GenericVector& x){
    obj=x.spec->clone(*x.obj);
    spec=obj->specializer();}
  GenericVector(GenericVector&& x){
    obj=x.obj; x.obj=nullptr;
    spec=obj->specializer();}
    
  ~GenericVector(){delete obj; delete spec;}

  Vector* obj=nullptr;
  VectorSpecializerBase* spec;


public:

  string classname() const {return "GenericVector";}
  virtual VectorSpecializerBase* specializer() const {
    return new VectorSpecializer<GenericVector>();}


public:

  template<class VECTOR>
  GenericVector static create(const VECTOR& v){
    GenericVector r; 
    r.obj=new VECTOR(v);
    return r;
  }


public: 

  bool isDense() const {return obj->isDense();}
  bool isSparse() const {return obj->isSparse();}
    
  bool operator==(const GenericVector& x) const {return spec->equals(*obj,*x.obj);}
  bool operator!=(const GenericVector& x) const {return !(spec->equals(*obj,*x.obj));}

  int nFilled() const {return obj->nFilled();}
  
  SCALAR& operator()(const int i){return (*obj)(i);}
  SCALAR operator()(const int i) const {return (*obj)(i);}
  SCALAR read(const int i) const {return obj->read(i);}
    
  void for_each(std::function<void(const INDEX, SCALAR&)> lambda) {obj->for_each(lambda);}
  void for_each(std::function<void(const INDEX, const SCALAR)> lambda) const {
    const_cast<const Vector*>(obj)->for_each(lambda);}

  GenericVector& operator+=(const SCALAR& x) {spec->increment(*obj,x); return *this;}
  GenericVector& operator-=(const SCALAR& x) {spec->decrement(*obj,x); return *this;}
  GenericVector& operator*=(const SCALAR& x) {spec->multiply_by(*obj,x); return *this;}
  GenericVector& operator/=(const SCALAR& x) {spec->divide_by(*obj,x); return *this;}

  GenericVector& operator+=(const GenericVector& x) {spec->increment(*obj,*x.obj); return *this;}
  GenericVector& operator-=(const GenericVector& x) {spec->decrement(*obj,*x.obj); return *this;}
  GenericVector& operator*=(const GenericVector& x) {spec->multiply_by(*obj,*x.obj); return *this;}
  GenericVector& operator/=(const GenericVector& x) {spec->divide_by(*obj,*x.obj); return *this;}

  SCALAR dot(const GenericVector& x) const {return spec->dot(*obj,*x.obj);}

  //GenericVector& operator*(const SCALAR x) const {return GenericVector(spec->times(*obj,x));}
  //GenericVector& operator/(const SCALAR x) const {return GenericVector(spec->by(*obj,x));}
  //GenericVector& operator+(const GenericVector& x) const {return GenericVector(spec->plus(*obj,*x.obj));}
  //GenericVector& operator-(const GenericVector& x) const {return GenericVector(spec->minus(*obj,*x.obj));}

  int nnz() const {return obj->nnz();}

  SCALAR max() const {return obj->max();}
  SCALAR max_abs() const {return obj->max();}
  int argmax() const {return obj->argmax();}
  int argmax_abs() const {return obj->argmax_abs();}
  
  SCALAR min() const {return obj->min();}
  SCALAR min_abs() const {return obj->min();}
  int argmin() const {return obj->argmin();}
  int argmin_abs() const {return obj->argmin_abs();}
  
  SCALAR sum() const {return obj->sum();}
  SCALAR norm1() const {return obj->norm1();}
  SCALAR norm2() const {return obj->norm2();}

  string str() const{return obj->str();}


};


inline ostream& operator<<(ostream& stream, const GenericVector& x){stream<<x.str(); return stream;}
  

}

#endif
