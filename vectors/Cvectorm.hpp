#ifndef _Cvectorm
#define _Cvectorm

#include "Mondrian_base.hpp"
#include "ThreadBank.hpp"
#include "Cvector.hpp"
#include "VectorView.hpp"
#include "Detachable.hpp"
#include "GivensRotation.hpp"

extern default_random_engine randomNumberGenerator;



namespace Mondrian{

  class AsCvectorm;


  class Cvectorm: public Cvector{
  public: 
    
    
    int maxthreads=1000;
    int nthreads=4;

    bool isMultithreaded() const{return true;}


  public: // constructors

    Cvectorm(): Cvector() {}

    Cvectorm(const int _n): Cvector(_n) {}

    Cvectorm(const int _n, const SCALAR* _array): 
      Cvector(_n,_array) {}
    
    Cvectorm(const initializer_list<SCALAR> list): 
      Cvector(list){}

    Cvectorm(const int _n, const initializer_list<ivpair> list): 
      Cvector(_n,list){}


  public: // copying

    Cvectorm(const Cvectorm& x): 
      Cvector(x,_NoWarn()){
      COPY_WARNING("Cvectorm");
    } 

    Cvectorm(const Cvectorm& x, const _NoWarn dummy): 
      Cvector(x.n, dummy){}

    Cvectorm(Cvectorm&& x): 
      Cvector(std::move(x),_NoWarn()){
      MOVE_WARNING("Cvectorm");
    }

    Cvectorm& operator=(const Cvectorm& x){
      Cvector::operator=(x);
      return *this;
    }
    
    Cvectorm& operator=(Cvectorm&& x){
      Cvector::operator=(std::move(x));
      return *this;
    }

    Cvectorm copy() const {
      return Cvector::copy();
    }
    
    Cvectorm shallow() const {
      return Cvector::shallow();
    }

    void assign(const Cvectorm& x){
      return Cvector::assign(x);
    }


  public: // downcasting 

    Cvectorm(const Cvector& x):
      Cvector(x,_NoWarn()){
      DOWNCASTCOPY_WARNING("Cvector","Cvectorm");
    }

    Cvectorm(Cvector&& x):
      Cvector(std::move(x)){}


  public: // polymorphism
    
    static string classname() {return "Cvectorm";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer<Cvectorm>();}


  public: // function mappings


    void pfor_each_filled(std::function<void(const INDEX, const SCALAR&)> lambda) const{
      int nblocks=::max(::min(nthreads,n/10),1);
      const int w=n/nblocks;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nblocks; I++){
	int imax=(I==nblocks-1) ? (n):((I+1)*w);
	threads.add([this,lambda](const int imin, const int imax){
	    for(int i=imin; i<imax; i++) lambda(i,array[i]);
	  },I*w,imax);
      }
    }

    
    void pfor_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      int nblocks=::max(::min(nthreads,n/10),1);
      const int w=n/nblocks;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nblocks; I++){
	int imax=(I==nblocks-1) ? (n):((I+1)*w);
	threads.add([this,lambda](const int imin, const int imax){
	    for(int i=imin; i<imax; i++) lambda(i,array[i]);
	  },I*w,imax);
      }
    }

    
    void pfor_each(std::function<void(const INDEX, const SCALAR&)> lambda) const{
      pfor_each_filled(lambda);}


    void pfor_each(std::function<void(const INDEX, SCALAR&)> lambda){
      pfor_each_filled(lambda);}


    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const TYPE, const TYPE)> accumulator, const TYPE& t0=TYPE()) const{
      int nblocks=::max(::min(nthreads,n/10),1);
      const int w=n/nblocks;
      TYPE dump[nblocks];
      {ThreadBank threads(maxthreads);
	for(int I=0; I<nblocks; I++){
	  int imax=(I==nblocks-1) ? (n):((I+1)*w);
	  threads.add([this,accumulator,t0](const int imin, const int imax, TYPE* target){
	      TYPE t(t0);
	      for(int i=imin; i<imax; i++)
		t=accumulator(t,array[i]);
	      *target=t;
	    },I*w,imax,&dump[I]);
	}
      }
      TYPE result(t0);
      for(int I=0; I<nblocks; I++)
	result=accumulator(result,dump[I]);
      return result;
    }

    
    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const SCALAR&)> lambda, 
      std::function<TYPE(const TYPE, const TYPE)> accumulator=std::plus<TYPE>(), 
      const TYPE& t0=TYPE()) const{
      
      int nblocks=::max(::min(nthreads,n/10),1);
      const int w=n/nblocks;
      TYPE dump[nblocks];
      {ThreadBank threads(maxthreads);
	for(int I=0; I<nblocks; I++){
	  int imax=(I==nblocks-1) ? (n):((I+1)*w);
	  threads.add([this,lambda,accumulator,t0](const int imin, const int imax, TYPE* target){
	      TYPE t(t0);
	      for(int i=imin; i<imax; i++)
		t=accumulator(t,lambda(array[i]));
	      *target=t;
	    },I*w,imax,&dump[I]);
	}
      }
      TYPE result(t0);
      for(int I=0; I<nblocks; I++)
	result=accumulator(result,dump[I]);
      return result;
    }

       
    INDEX find_best(std::function<bool(const SCALAR&, const SCALAR&)> selector) const{
      int nblocks=::max(::min(nthreads,n/10),1);
      const int w=n/nblocks;
      SCALAR bestv[nblocks];
      INDEX  besti[nblocks];
      {ThreadBank threads(maxthreads);
	for(int I=0; I<nblocks; I++){
	  int imax=(I==nblocks-1) ? (n):((I+1)*w);
	  threads.add([this,selector](const int imin, const int imax, SCALAR* vtarget, INDEX* itarget){
	      INDEX besti=imin;
	      SCALAR bestv=array[imin];
	      for(int i=imin+1; i<imax; i++)
		if(selector(bestv,array[i])){besti=i; bestv=array[i];}
	      *vtarget=bestv;
	      *itarget=besti;
	    },I*w,imax,&bestv[I],&besti[I]);
	}
      }
      INDEX best_i=besti[0];
      SCALAR best_v=bestv[0];
      for(int I=1; I<nblocks; I++)
	if(selector(best_v,bestv[I])){best_i=besti[I]; best_v=bestv[I];}
      return best_i;
    }
       
  
  public: // in-place arithmetic
    

    void operator+=(const SCALAR& x) {pfor_each([x](INDEX i, SCALAR& v){v+=x;});}
    void operator-=(const SCALAR& x) {pfor_each([x](INDEX i, SCALAR& v){v-=x;});}
    void operator/=(const SCALAR& x) {pfor_each([x](INDEX i, SCALAR& v){v/=x;});}
    void operator*=(const SCALAR& x) {pfor_each([x](INDEX i, SCALAR& v){v*=x;});}
    
    void operator+=(const Cvectorm& x) {pfor_each([&x](INDEX i, SCALAR& v){v+=x.array[i];});}
    void operator-=(const Cvectorm& x) {pfor_each([&x](INDEX i, SCALAR& v){v-=x.array[i];});}
    void operator*=(const Cvectorm& x) {pfor_each([&x](INDEX i, SCALAR& v){v*=x.array[i];});}
    void operator/=(const Cvectorm& x) {pfor_each([&x](INDEX i, SCALAR& v){v/=x.array[i];});}

    Cvectorm& add(const Cvectorm& x) {assert(n==x.n); for(int i=0; i<n; i++) array[i]+=x.array[i]; return *this;}
    Cvectorm& add(const Cvectorm& x, const SCALAR c) {assert(n==x.n); for(int i=0; i<n; i++) array[i]+=c*x.array[i]; return *this;}


  public: // scalar-valued arithmetic

    SCALAR dot(const Cvectorm& x) const{assert(x.n==n);
      SCALAR t=0; for(int i=0; i<n; i++) t+=array[i]*x.array[i]; return t;}


  public: // vector-valued arithmetic 
    
    Cvectorm mult(const SCALAR& x) const {
      Cvectorm v(n); for(int i=0; i<n; i++) v.array[i]=array[i]*x; return v;}

    Cvectorm plus(const Cvectorm& x) const {assert(x.n==n); 
      Cvectorm v(n); for(int i=0; i<n; i++) v.array[i]=array[i]+x.array[i]; return v;}
    Cvectorm minus(const Cvectorm& x) const {assert(x.n==n); 
      Cvectorm v(n); for(int i=0; i<n; i++) v.array[i]=array[i]-x.array[i]; return v;}

    Cvectorm operator*(const SCALAR& x) const {return this->mult(x);}
    Cvectorm operator/(const SCALAR& x) const {return this->mult(1.0/x);}
    Cvectorm operator+(const Cvectorm& x) const {return this->plus(x);}
    Cvectorm operator-(const Cvectorm& x) const {return this->minus(x);}


  public: // scalar valued methods
    
    SCALAR max() const {assert(n>0);
      return accumulate<SCALAR>([](const SCALAR& a, const SCALAR& b){return std::max(a,b);},array[0]);}
    SCALAR max_abs() const {assert(n>0);
      return accumulate<SCALAR>([](const SCALAR& a, const SCALAR& b){return std::max(a,fabs(b));},fabs(array[0]));}
    INDEX argmax() const {assert(n>0); 
      return find_best([](const SCALAR& a, const SCALAR& b){return (b>a);});}
    INDEX argmax_abs() const {assert(n>0);
      return find_best([](const SCALAR& a, const SCALAR& b){return (fabs(b)>fabs(a));});}

    SCALAR min() const {assert(n>0);
      return accumulate<SCALAR>([](const SCALAR& a, const SCALAR& b){return std::min(a,b);},array[0]);}
    SCALAR min_abs() const {assert(n>0);
      return accumulate<SCALAR>([](const SCALAR& a, const SCALAR& b){return std::min(a,fabs(b));},fabs(array[0]));}
    INDEX argmin() const {assert(n>0);
      return find_best([](const SCALAR& a, const SCALAR& b){return (b<a);});}
    INDEX argmin_abs() const {assert(n>0);
      return find_best([](const SCALAR& a, const SCALAR& b){return (fabs(b)<fabs(a));});}

    SCALAR sum() const{
      return accumulate<SCALAR>([](const SCALAR& a, const SCALAR& b){return a+b;});}
    SCALAR norm1() const{
      return accumulate<SCALAR>([](const SCALAR& a, const SCALAR& b){return a+fabs(b);});}
    SCALAR norm2() const{ // TODO!
      return accumulate<SCALAR>([](const SCALAR& a){return a*a;});}
    SCALAR diff2(const Cvector& x) const{
      assert(x.n==n); SCALAR t=0; for(int i=0; i<n; i++) t+=(array[i]-x.array[i])*(array[i]-x.array[i]); return t;}
    int nnz() const{
      return accumulate<int>([](const SCALAR& a){return (a!=0);});}


  public: // Python interface

    Cvectorm(double* numpyInDblArray, int numpyInSize): Cvectorm(numpyInSize){
      std::copy(numpyInDblArray,numpyInDblArray+n,array);
    }  

  };



  class AsCvectorm: public Cvectorm{
  public:
    AsCvectorm()=delete;
    AsCvectorm& operator=(const AsCvectorm& x)=delete;
    AsCvectorm(Cvector& x): Cvectorm(x.shallow()){}
    AsCvectorm(const Cvector& x): Cvectorm(x.shallow()){}
  };

  AsCvectorm as_multithreaded(Cvector& x){
    return AsCvectorm(x);}

  const AsCvectorm as_multithreaded(const Cvector& x){
    return AsCvectorm(x);}

  Cvectorm as_multithreaded(Cvector&& x){
    return Cvectorm(std::move(x));}




  inline ostream& operator<<(ostream& stream, const Cvectorm& v){stream<<v.str(); return stream;}






} // namespace Mondrian







#endif

//bool equals(const Vector& x) const {return (*this)==(static_cast<const Cvector&>(x));}

//  void for_each(std::function<void(const INDEX, SCALAR&)> lambda) {for(int i=0; i<n; i++) lambda(i,array[i]);}
//  void for_each(std::function<void(const INDEX, const SCALAR)> lambda) const {for(int i=0; i<n; i++) lambda(i,array[i]);}

//  int size() const {return n;}


    //Cvector* clone() const {return new Cvector(*this);}
    /*
  public: // KpointOperators

    Cvector& apply(const KpointOperator& Q){
      for(int i=0; i<Q.k; i++) assert(Q.map(i)<n);
      for(int i=0; i<Q.k; i++) {Q.tempp[i]=&array[Q.map(i)]; Q.temp[i]=*Q.tempp[i];}
      for(int i=0; i<Q.k; i++){
	double s=0; for(int j=0; j<Q.k; j++) s+=Q.q[j*Q.k+i]*Q.temp[j]; *Q.tempp[i]=s;}
      return *this;
    }
    
    Cvector& applyT(const KpointOperator& Q){
      for(int i=0; i<Q.k; i++) assert(Q.map(i)<n);
      for(int i=0; i<Q.k; i++) {Q.tempp[i]=&array[Q.map(i)]; Q.temp[i]=*Q.tempp[i];}
      for(int i=0; i<Q.k; i++){
	double s=0; for(int j=0; j<Q.k; j++) s+=Q.q[i*Q.k+j]*Q.temp[j]; *Q.tempp[i]=s;}
      return *this;
    }
    */
   /*
  public: // named constructors
    
    static Cvectorm Zero(const int n) {
      Cvectorm v(n); for(int i=0; i<n; i++) v.array[i]=0; return v;}
    
    static Cvectorm Filled(const int n, const SCALAR t){
      Cvectorm v(n); for(int i=0; i<n; i++) v.array[i]=t; return v;}
        
    static Cvectorm Uniform(const int n){
      Cvectorm v(n); uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;}

    static Cvectorm Gaussian(const int n){
      Cvectorm v(n); normal_distribution<SCALAR> distr;
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;}

    static Cvectorm Bernoulli(const int n, const double p=0.5){
      Cvectorm v(n); bernoulli_distribution distr(p);
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;}
    */
