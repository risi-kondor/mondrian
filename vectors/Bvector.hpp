#ifndef _Bvector
#define _Bvector

#include "Cvector.hpp"

extern default_random_engine randomNumberGenerator;

namespace Mondrian{

  class Bvector: public Cvector{
  public: 
    
    Bvector(const Bvector& x): Bvector(x.n){
      COPY_WARNING("Bvector");
      std::copy(x.array,x.array+n,array);
    } 

    Bvector(Bvector&& x): Cvector(x.n){
      MOVE_WARNING("Bvector");
      array=x.array; x.array=nullptr;
    }

    Bvector& operator=(const Bvector& x){
      ASSIGN_WARNING("Bvector");
      n=x.n; delete[] array; array=new SCALAR[n]; 
      std::copy(x.array,x.array+n,array); 
      return *this;
    }
    
    Bvector& operator=(Bvector&& x){
      MOVEASSIGN_WARNING("Bvector");
      n=x.n; delete[] array; array=x.array; x.array=nullptr; 
      return *this;
    }

    Bvector copy() const {
      Bvector v(n); 
      std::copy(array,array+n,v.array); 
      return v;
    }
    
    Bvector shallow() const {Bvector v(n); v.array=array; return v;}
    void assign(const Bvector& x){assert(x.n==n); std::copy(x.array,x.array+n,array);}
    void detach(){array=nullptr;}

  public: // constructors
    
    using Cvector::Cvector;

  public: // named constructors
    
    static Bvector Zero(const int n) {
      Bvector v(n); for(int i=0; i<n; i++) v.array[i]=0; return v;}
    
    static Bvector Filled(const int n, const SCALAR t){
      Bvector v(n); for(int i=0; i<n; i++) v.array[i]=t; return v;}
    
    static Bvector Gaussian(const int n){
      Bvector v(n); normal_distribution<SCALAR> distr;
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;
    }

  public: // in-place operators
    
    Bvector& operator*=(const SCALAR& x) {cblas_dscal(n,x,array,incr); return *this;}
    Bvector& operator/=(const SCALAR& x) {cblas_dscal(n,1.0/x,array,incr); return *this;}
    
    Bvector& operator+=(const Bvector& x) {assert(x.n==n); blas_daxpy(n,1.0,x.array,incr,array,incr); return *this;}
    Bvector& operator-=(const Bvector& x) {assert(x.n==n); blas_daxpy(n,-1.0,x.array,incr,array,incr);  return *this;}

    Bvector& add(const Bvector& x) {assert(n==x.n); blas_daxpy(n,1.0,x.array,incr,array,incr); return *this;}
    Bvector& add(const Bvector& x, const SCALAR c) {assert(n==x.n); blas_daxpy(n,c,x.array,incr,array,incr); return *this;}


  public: // scalar-valued operations 

    SCALAR dot(const Bvector& x) const{assert(x.n==n); return blas_ddot(n, array, incr, x.array, incr);}

  public: // vector-valued operations 
    
    Bvector operator*(const SCALAR& x) const {Bvector v=copy(); cblas_dscal(n, x, v.array, incr); return v;}
    Bvector operator/(const SCALAR& x) const {Bvector v=copy(); cblas_dscal(n,1.0/x,v.array,incr); return v;}
    
    Bvector operator+(const Bvector& x) const{assert(x.n==n); Bvector v=copy(); blas_daxpy(n,1.0,x.array,incr,v.array,incr); return v;}
    Bvector operator-(const Bvector& x) const{assert(x.n==n); Bvector v=copy(); blas_daxpy(n,-1.0,x.array,incr,v.array,incr); return v;}


  public: // scalar valued methods
    
    SCALAR norm1() const {return blas_dasum(n, array, incr);}
    SCALAR norm2() const {double res=blas_dnrm2(n, array, incr); return res*res;}

    SCALAR max_abs() const {return array[argmax_abs()];}
    INDEX argmax_abs() const{
      if(n==0) return 0; 
      return blas_idamax(n, array, incr);
    }

    SCALAR min_abs() const {return array[argmin_abs()];}
    INDEX argmin_abs() const{ 
      if(n==0) return 0; 
      return blas_idamin(n, array, incr);
    }
    
public:
  int incr = 1;

};

}

#endif
