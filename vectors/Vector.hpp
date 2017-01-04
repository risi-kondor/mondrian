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


#ifndef _Vector
#define _Vector

#include "Mondrian_base.hpp"
#include <math.h>

namespace Mondrian{

  class VectorSpecializerBase;

  class GivensRotation;
  class KpointOperator;

  class Vector{
  public:

    int n;


  public:
    
    Vector(): n(0){}
    Vector(const int _n=0): n(_n){}
    virtual ~Vector(){};



  public: // polymorphism

    //static string classname() const=0;
    virtual VectorSpecializerBase* specializer() const=0;

    
  public: // attributes

    virtual bool isSparse() const {return false;}
    virtual bool isMultithreaded() const {return false;}

    
  public: // comparisons
    
    template<class VECTOR2>
    bool operator==(const VECTOR2& x) const{
      if(x.n!=n) return false;
      for(int i=0; i<n; i++) if((*this)(i)!=x(i)) return false;
      return true;}

    template<class VECTOR2>
    bool operator!=(const VECTOR2& x) const {return !((*this)==x);}


  public: // element access
        
    virtual SCALAR operator()(const int i) const=0; 
    virtual SCALAR read(const int i) const=0; 

    virtual void set(const int i, const SCALAR x)=0;
    virtual SCALAR& operator()(const int i)=0; 

    virtual bool isFilled(const int i) const=0;
    virtual int nFilled() const=0;


  public: // iterators 
    
    virtual void for_each(std::function<void(const INDEX, SCALAR&)> lambda){
      for(int i=0; i<n; i++) lambda(i,(*this)(i));}
    virtual void for_each(std::function<void(const INDEX, SCALAR)> lambda) const{
      for(int i=0; i<n; i++) lambda(i,this->read(i));}

    virtual void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      for(int i=0; i<n; i++) if(this->isFilled(i)) lambda(i,(*this)(i));}
    virtual void for_each_filled(std::function<void(const INDEX, SCALAR)> lambda) const{
      for(int i=0; i<n; i++) if(this->isFilled(i)) lambda(i,this->read(i));}


  public: // conversions 
  
    //virtual operator Cvector() const;

    
  public: // in-place operators
    
    virtual void operator+=(const SCALAR& x) {for(int i=0; i<n; i++) (*this)(i)+=x;}
    virtual void operator-=(const SCALAR& x) {for(int i=0; i<n; i++) (*this)(i)-=x;}
    virtual void operator*=(const SCALAR& x) {for(int i=0; i<n; i++) (*this)(i)*=x;}
    virtual void operator/=(const SCALAR& x) {for(int i=0; i<n; i++) (*this)(i)/=x;}
    
    template<class VECTOR2>
    void operator+=(const VECTOR2& x) {assert(x.n==n); for(int i=0; i<n; i++) (*this)(i)+=x(i);}
    template<class VECTOR2>
    void operator-=(const VECTOR2& x) {assert(x.n==n); for(int i=0; i<n; i++) (*this)(i)-=x(i);}
    template<class VECTOR2>
    void operator*=(const VECTOR2& x) {assert(x.n==n); for(int i=0; i<n; i++) (*this)(i)*=x(i);}
    template<class VECTOR2>
    void operator/=(const VECTOR2& x) {assert(x.n==n); for(int i=0; i<n; i++) (*this)(i)/=x(i);}

    virtual void exp(const SCALAR x){
      for(int i=0; i<n; i++) (*this)(i)=::exp(x*(*this)(i));
    }

  public: // scalar valued arithmetic

    template<class VECTOR2>
    SCALAR dot(const VECTOR2& x) const{assert(x.n==n);
      SCALAR t=0; for(int i=0; i<n; i++) t+=(*this)(i)*x(i); return t;}


  public: // matrix valued methods

    //Cmatrix outer(const Vector& x) const;
    

  public: // scalar valued methods
    
    virtual SCALAR max() const=0;
    virtual SCALAR max_abs() const=0;
    virtual INDEX argmax() const=0;
    virtual INDEX argmax_abs() const=0;

    virtual SCALAR min() const=0;
    virtual SCALAR min_abs() const=0;
    virtual INDEX argmin() const=0;
    virtual INDEX argmin_abs() const=0;

    virtual SCALAR sum() const=0;
    virtual SCALAR norm1() const=0;
    virtual SCALAR norm2() const=0;

    virtual int nnz() const=0;


  public: // operators 

    // trouble with covariant return types in virtual functions 
    // virtual Vector& apply(const GivensRotation& Q)=0;
    // virtual Vector& apply(const KpointOperator& Q)=0;
    // virtual Vector& applyT(const GivensRotation& Q)=0;
    // virtual Vector& applyT(const KpointOperator& Q)=0;


  public: // I/O
    
    //virtual saveTo(MatrixOF& file) const =0;

    virtual string str() const{
      ostringstream oss; oss<<"(";
      //oss.precision(3); 
    //oss.setf(ios_base::fixed, ios_base::floatfield);
      for(int i=0; i<n-1; i++){oss<<(*this)(i)<<",";}
      if(n>0) oss<<(*this)(n-1); oss<<")";
      return oss.str();
    }

    //virtual void print() const{
    //  cout<<str()<<endl;
    //}
   
  public: // Python 

    const char* __str__() {
      ostringstream ostream; 
      ostream<<str()<<endl; 
      return ostream.str().c_str();
    }


};



  class VectorSpecializerBase{
  public:

    virtual ~VectorSpecializerBase(){};
    virtual Vector* clone(const Vector& x) const=0;
    virtual bool equals(const Vector& x, const Vector& y) const=0;
    virtual void increment(Vector& x, const SCALAR y) const=0;
    virtual void decrement(Vector& x, const SCALAR y) const=0;
    virtual void multiply_by(Vector& x, const SCALAR y) const=0;
    virtual void divide_by(Vector& x, const SCALAR y) const=0;
    virtual void increment(Vector& x, const Vector& y) const=0;
    virtual void decrement(Vector& x, const Vector& y) const=0;
    virtual void multiply_by(Vector& x, const Vector& y) const=0;
    virtual void divide_by(Vector& x, const Vector& y) const=0;
    virtual SCALAR dot(const Vector& x, const Vector& y) const=0;
    //virtual Vector times(const Vector& x, const SCALAR y) const=0;
    //virtual Vector by(const Vector& x, const SCALAR y) const=0;
    //virtual Vector plus(const Vector& x, const Vector& y) const=0;
    //virtual Vector minus(const Vector& x, const Vector& y) const=0;

  };



  template<class VECTOR>
  class VectorSpecializer: public VectorSpecializerBase{
  public:

    VECTOR* clone(const Vector& x) const {
      return new VECTOR(dynamic_cast<const VECTOR&>(x));}
    bool equals(const Vector& x, const Vector& y) const{
      return dynamic_cast<const VECTOR&>(x)==(dynamic_cast<const VECTOR&>(y));}
    void increment(Vector& x, const SCALAR y) const {dynamic_cast<VECTOR&>(x)+=y;}
    void decrement(Vector& x, const SCALAR y) const {dynamic_cast<VECTOR&>(x)-=y;}
    void multiply_by(Vector& x, const SCALAR y) const {dynamic_cast<VECTOR&>(x)*=y;}
    void divide_by(Vector& x, const SCALAR y) const {dynamic_cast<VECTOR&>(x)/=y;}
    void increment(Vector& x, const Vector& y) const {dynamic_cast<VECTOR&>(x)+=dynamic_cast<const VECTOR&>(y);}
    void decrement(Vector& x, const Vector& y) const {dynamic_cast<VECTOR&>(x)-=dynamic_cast<const VECTOR&>(y);}
    void multiply_by(Vector& x, const Vector& y) const {dynamic_cast<VECTOR&>(x)*=dynamic_cast<const VECTOR&>(y);}
    void divide_by(Vector& x, const Vector& y) const {dynamic_cast<VECTOR&>(x)/=dynamic_cast<const VECTOR&>(y);}
    SCALAR dot(const Vector& x, const Vector& y) const {return dynamic_cast<const VECTOR&>(x).dot(dynamic_cast<const VECTOR&>(y));}
    //VECTOR times(const Vector& x, const SCALAR y) const {return dynamic_cast<const VECTOR&>(x)*y;}
    //Vector by(const Vector& x, const SCALAR y) const {return dynamic_cast<const VECTOR&>(x)/y;}
    //Vector plus(const Vector& x, const Vector& y) const {return dynamic_cast<const VECTOR&>(x)+dynamic_cast<const VECTOR&>(y);}
    //Vector minus(const Vector& x, const Vector& y) const {return dynamic_cast<const VECTOR&>(x)-dynamic_cast<const VECTOR&>(y);}
    
  };

  
  inline ostream& operator<<(ostream& stream, const Vector& x){stream<<x.str(); return stream;}


} // namespace Mondrian

#endif


