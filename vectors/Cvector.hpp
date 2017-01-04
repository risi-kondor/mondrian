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


#ifndef _Cvector
#define _Cvector

#include "Mondrian_base.hpp"
#include "Vector.hpp"
#include "VectorView.hpp"
#include "Detachable.hpp"
#include "GivensRotation.hpp"

extern default_random_engine randomNumberGenerator;



namespace Mondrian{

  class Cvector: public Vector, public Detachable{
  public: 
    
    SCALAR* array;


  public:

    Cvector(): Vector(0) {array=NULL;}

    Cvector(const int _n): Vector(_n) {array=new SCALAR[n];}

    Cvector(const int _n, const SCALAR* _array): Cvector(_n) {
      std::copy(_array,_array+n,array);}
    
    Cvector(const initializer_list<SCALAR> list): Cvector(list.size()){
      array=new SCALAR[n]; int i=0; for(SCALAR v:list) array[i++]=v;}

    Cvector(const int _n, const initializer_list<ivpair> list): Cvector(_n){
      for(int i=0; i<n; i++) array[i]=0;
      for(const ivpair& v:list) {assert(v.first<n); array[v.first]=v.second;}
    }

    ~Cvector(){delete[] array;}

    
  public: // copying

    Cvector(const Cvector& x): Cvector(x.n){
      COPY_WARNING("Cvector");
      std::copy(x.array,x.array+n,array);
    } 

    Cvector(const Cvector& x, const _NoWarn dummy): Cvector(x.n){
      std::copy(x.array,x.array+n,array);
    } 

    Cvector(Cvector&& x): Vector(x.n){
      MOVE_WARNING("Cvector");
      array=x.array; x.array=nullptr;
    }

    Cvector& operator=(const Cvector& x){
      ASSIGN_WARNING("Cvector");
      n=x.n; delete[] array; array=new SCALAR[n]; 
      std::copy(x.array,x.array+n,array); 
      return *this;
    }
    
    Cvector& operator=(Cvector&& x){
      MOVEASSIGN_WARNING("Cvector");
      n=x.n; delete[] array; array=x.array; x.array=nullptr; 
      return *this;
    }

    Cvector copy() const {
      Cvector v(n); 
      std::copy(array,array+n,v.array); 
      return v;
    }
    
    Cvector shallow() const {
      Cvector v(n); v.array=array; return v;}

    void assign(const Cvector& x){
      assert(x.n==n); std::copy(x.array,x.array+n,array);}

    void detach(){array=nullptr;}


  public: // named constructors
    
    static Cvector Zero(const int n) {
      Cvector v(n); for(int i=0; i<n; i++) v.array[i]=0; return v;}
    
    static Cvector Filled(const int n, const SCALAR t){
      Cvector v(n); for(int i=0; i<n; i++) v.array[i]=t; return v;}
        
    static Cvector Uniform(const int n){
      Cvector v(n); uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;}

    static Cvector Gaussian(const int n){
      Cvector v(n); normal_distribution<SCALAR> distr;
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;}

    static Cvector Bernoulli(const int n, const double p=0.5){
      Cvector v(n); bernoulli_distribution distr(p);
      for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator); 
      return v;}
    

  public: // polymorphism
    
    static string classname() {return "Cvector";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer<Cvector>();}


  public: // attributes

    bool isSparse() const {return false;}


  public: // comparisons
    
    bool operator==(const Cvector& x) const{
      if(x.n!=n) return false;
      for(int i=0; i<n; i++) if(array[i]!=x.array[i]) return false;
      return true;}


  public: // element access 
    
    SCALAR operator()(const int i) const {assert(i<n); return array[i];}
    SCALAR read(const int i) const {assert(i<n); return array[i];}

    SCALAR& operator()(const int i){assert(i<n); return array[i];}
    void set(const int i, const SCALAR v){assert(i<n); array[i]=v;}
    void set_msafe(const int i, const SCALAR v){assert(i<n); array[i]=v;}
    SCALAR* ptr(const int i){assert(i<n); return &array[i];}
    
    bool isFilled(const int i) const {assert(i<n); return true;}
    int nFilled() const {return n;}


  public: // sparse vector methods

    void insert(const int i, const SCALAR& v){assert(i<n); array[i]=v;}
    void sort(){}
    void tidy(){}


  public: // iterators 
    
    void for_each(std::function<void(const INDEX, SCALAR&)> lambda){
      for(int i=0; i<n; i++) lambda(i,array[i]);}
    void for_each(std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(int i=0; i<n; i++) lambda(i,array[i]);}
    void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      for(int i=0; i<n; i++) lambda(i,array[i]);}
    void for_each_filled(std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(int i=0; i<n; i++) lambda(i,array[i]);}

    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const TYPE, const TYPE)> accumulator, const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(int i=0; i<n; i++)
	t=accumulator(t,array[i]);
      return t;
    }

    template<class TYPE>
    TYPE accumulate(std::function<TYPE(const SCALAR&)> lambda, 
      std::function<TYPE(const TYPE, const TYPE)> accumulator=std::plus<TYPE>(), 
      const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(int i=0; i<n; i++)
	t=accumulator(t,lambda(array[i]));
      return t;
    }

    INDEX find_best(std::function<bool(const SCALAR&, const SCALAR&)> selector) const{
      assert(n>0);
      INDEX besti=0;
      SCALAR bestv=array[0];
      for(int i=1; i<n; i++)
	if(selector(bestv,array[i])){besti=i; bestv=array[i];}
      return besti;
    }


  public: // views

    VectorView<Cvector> operator()(const IndexMap& phi) {
      return VectorView<Cvector>(*this,phi);}
    //const VectorView<Cvector> operator()(const IndexMap& phi) const {return VectorView<Cvector>(*this,phi);}

    Detached<Cvector> subvector_view(const int i, const int _n){
      Detached<Cvector> R; R.n=_n; R.array=array+i; return R;}


  public: // conversions
  
    //Cvector(const BlockedVector<Cvector>& w); -> Blocked Vector
    //Cvector(BlockedVector<Cvector>&& w);

    //#ifdef _withEigen
    //Cvector(const EigenVectorXdAdaptor& X);
    //operator EigenVectorXdAdaptor() const;
    //#endif


  public: // remappings

    Cvector remap(const IndexMap& map) const{assert(map.nsource==n);
      Cvector r(n); for(int i=0; i<n; i++) r(i)=array[map(i)]; return r;}
    Cvector remap(const Inverse<IndexMap>& imap) const{assert(imap.obj.nsource==n);
      Cvector r(n); for(int i=0; i<n; i++) r(imap.obj(i))=array[i]; return r;}
  

  public: // in-place arithmetic
    
    void operator+=(const SCALAR& x) {for(int i=0; i<n; i++) array[i]+=x;}
    void operator-=(const SCALAR& x) {for(int i=0; i<n; i++) array[i]-=x;}
    void operator*=(const SCALAR& x) {for(int i=0; i<n; i++) array[i]*=x;}
    void operator/=(const SCALAR& x) {for(int i=0; i<n; i++) array[i]/=x;}
    
    void operator+=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]+=x.array[i];}
    void operator-=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]-=x.array[i];}
    void operator*=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]*=x.array[i];}
    void operator/=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]/=x.array[i];}

    Cvector& add(const Cvector& x) {assert(n==x.n); for(int i=0; i<n; i++) array[i]+=x.array[i]; return *this;}
    Cvector& add(const Cvector& x, const SCALAR c) {assert(n==x.n); for(int i=0; i<n; i++) array[i]+=c*x.array[i]; return *this;}


  public: // scalar-valued arithmetic

    SCALAR dot(const Cvector& x) const{assert(x.n==n);
      SCALAR t=0; for(int i=0; i<n; i++) t+=array[i]*x.array[i]; return t;}


  public: // vector-valued arithmetic 
    
    Cvector mult(const SCALAR& x) const {
      Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]*x; return v;}

    Cvector plus(const Cvector& x) const {assert(x.n==n); 
      Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]+x.array[i]; return v;}
    Cvector minus(const Cvector& x) const {assert(x.n==n); 
      Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]-x.array[i]; return v;}

    Cvector operator*(const SCALAR& x) const {return this->mult(x);}
    Cvector operator/(const SCALAR& x) const {return this->mult(1.0/x);}
    Cvector operator+(const Cvector& x) const {return this->plus(x);}
    Cvector operator-(const Cvector& x) const {return this->minus(x);}


  public: // scalar valued methods
    
    SCALAR max() const{
      SCALAR t=array[0]; for(int i=1; i<n; i++) if(array[i]>t) t=array[i]; return t;}
    SCALAR max_abs() const{
      SCALAR t=fabs(array[0]); for(int i=1; i<n; i++) if(fabs(array[i])>t) t=array[i]; return t;}
    INDEX argmax() const{ 
      if(n==0) return 0; int best=0; SCALAR max=array[best];  
      for(int i=0; i<n; i++) if(array[i]>max) {best=i; max=array[i];}
      return best;}
    INDEX argmax_abs() const{
      if(n==0) return 0; int best=0; SCALAR max=fabs(array[best]);  
      for(int i=0; i<n; i++) if(fabs(array[i])>max) {best=i; max=fabs(array[i]);}
      return best;}

    SCALAR min() const{
      SCALAR t=array[0]; for(int i=1; i<n; i++) if(array[i]<t) t=array[i]; return t;}
    SCALAR min_abs() const{
      SCALAR t=fabs(array[0]); for(int i=1; i<n; i++) if(fabs(array[i])<t) t=array[i]; return t;}
    INDEX argmin() const{ 
      if(n==0) return 0; int best=0; SCALAR min=array[best];  
      for(int i=0; i<n; i++) if(array[i]<min) {best=i; min=array[i];}
      return best;}
    INDEX argmin_abs() const{ 
      if(n==0) return 0; int best=0; SCALAR min=fabs(array[best]);  
      for(int i=0; i<n; i++) if(fabs(array[i])<min) {best=i; min=array[i];}
      return best;}
    
    SCALAR sum() const{SCALAR t=0; for(int i=0; i<n; i++) t+=array[i]; return t;}
    SCALAR norm1() const {SCALAR t=0; for(int i=0; i<n; i++) t+=fabs(array[i]); return t;}
    SCALAR norm2() const {SCALAR t=0; for(int i=0; i<n; i++) t+=array[i]*array[i]; return t;}
    SCALAR diff2(const Cvector& x) const {assert(x.n==n); 
      SCALAR t=0; for(int i=0; i<n; i++) t+=(array[i]-x.array[i])*(array[i]-x.array[i]); return t;}
    int nnz() const{int t=0; for(int i=0; i<n; i++) if (array[i]!=0) t++; return t;}


  public: // Givens rotations 

    Cvector& apply(const GivensRotation& Q){
      assert(Q.i1<n); assert(Q.i2<n);
      SCALAR x1=array[Q.i1]; 
      SCALAR x2=array[Q.i2]; 
      array[Q.i1]=Q.cos*x1-Q.sin*x2;
      array[Q.i2]=Q.sin*x1+Q.cos*x2; 
      return *this;
    }

    Cvector& applyT(const GivensRotation& Q){
      assert(Q.i1<n); assert(Q.i2<n);
      SCALAR x1=array[Q.i1]; 
      SCALAR x2=array[Q.i2]; 
      array[Q.i1]=Q.cos*x1+Q.sin*x2;
      array[Q.i2]=-Q.sin*x1+Q.cos*x2; 
      return *this;
    }


  public: // KpointOp
    
    template<int k>
    Cvector& apply(const KpointOp<k>& Q){
      SCALAR temp[k];
      SCALAR* tempp[k];
      for(int i=0; i<k; i++) assert(Q.map(i)<n);
      for(int i=0; i<k; i++) {tempp[i]=&array[Q.map(i)]; temp[i]=*tempp[i];}
      for(int i=0; i<k; i++){
	double s=0; 
	for(int j=0; j<k; j++) 
	  s+=Q.q[j*k+i]*temp[j]; 
	*tempp[i]=s;
      }
      return *this;
    }
    
    template<int k>
    Cvector& applyT(const KpointOp<k>& Q){
      SCALAR temp[k];
      SCALAR* tempp[k];
      for(int i=0; i<k; i++) assert(Q.map(i)<n);
      for(int i=0; i<k; i++) {tempp[i]=&array[Q.map(i)]; temp[i]=*tempp[i];}
      for(int i=0; i<k; i++){
	double s=0; 
	for(int j=0; j<k; j++) 
	  s+=Q.q[i*k+j]*temp[j]; 
	*tempp[i]=s;
      }
      return *this;
    }
    

  public:

    //Cvector(MatrixIF& file);
    //virtual saveTo(MatrixOF& file) const=0;

    //void print() const{
    //  cout<<str()<<endl;}

  public: // Python interface

    Cvector(double* numpyInDblArray, int numpyInSize): Cvector(numpyInSize){
      std::copy(numpyInDblArray,numpyInDblArray+n,array);
    }

    void np(double** numpyOutArray1, int* numpyOutLen1){
      *numpyOutLen1=n;
      *numpyOutArray1=new double[n];
      std::copy(array,array+n,*numpyOutArray1);
    }


};





inline ostream& operator<<(ostream& stream, const Cvector& v){stream<<v.str(); return stream;}






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
