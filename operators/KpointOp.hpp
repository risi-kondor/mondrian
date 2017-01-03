#ifndef _KpointOp
#define _KpointOp

#include "Mondrian_base.hpp"
#include "IndexMap.hpp"
#include "Cmatrix.hpp"
#include "CmatrixLA.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

  template<int k>
  class KpointOp{
  public:
    
    //int k;
    double* q;
    //double* temp;
    //SCALAR** tempp;
    IndexMap map;

    
  public:

    //KpointOp()=delete;

    KpointOp(): map(k){
      q=new double[k*k];}
    
    KpointOp(const IndexMap& _map): map(_map){
      assert(map.nsource==k); q=new double[k*k];}
    
    KpointOp(const IndexMap& _map, const Cmatrix& M): map(_map){
      assert(map.nsource==k); assert(M.nrows==k); assert(M.ncols==k);
      q=new double[k*k];  std::copy(M.array,M.array+k*k,q);
    }

    KpointOp(IndexMap&& _map, Cmatrix&& M): map(std::move(_map)){
      assert(map.nsource==k); assert(M.nrows==k); assert(M.ncols==k);
      q=M.array; M.array=nullptr;
    }

    ~KpointOp() {delete[] q;}


  public: // copying 
      
    KpointOp(const KpointOp<k>& x): map(x.map){
      q=new double[k*k]; for(int i=0; i<k*k; i++) q[i]=x.q[i];
    }

    KpointOp(KpointOp<k>&& x): map(x.map){
      q=x.q; x.q=nullptr;
    }

    KpointOp& operator=(KpointOp<k>& x){
      delete[] q; //delete[] temp;
      k=x.k; map=x.map;
      q=new double[k*k]; for(int i=0; i<k*k; i++) q[i]=x.q[i];
      return *this;
    }

    KpointOp& operator=(KpointOp<k>&& x){
      delete[] q; //delete[] temp; 
      map=std::move(map);
      q=x.q; x.q=nullptr;
      return *this;
    }
  

  public: // conversions 

    operator Cmatrix() const{
      return Cmatrix(k,k,q);
    }
    
    
  public: // named constructors

    static KpointOp<k> RandomRotation(const int n){
      KpointOp<k> r(IndexMap::Random(k,n),Cmatrix::Gaussian(k,k));
      for(int j=0; j<k; j++){
	//Cvector vj=Cvector(k,r.q+j*k);
	Detached<Cvector> vj=as_Cmatrix(k,k,r.q).view_of_column(j);
	for(int i=0; i<j; i++){
	  Cvector vi=Cvector(k,r.q+i*k);
	  vj-=vi*(vj.dot(vi));
	}
	vj/=sqrt(vj.norm2());
      }
      return r;
    }


  public: // operations

    void applyTo(SCALAR** vptr) const{
      SCALAR temp[k];
      for(int i=0; i<k; i++) temp[i]=*vptr[i];
      for(int i=0; i<k; i++){
	double s=0; 
	for(int j=0; j<k; j++) 
	  s+=q[j*k+i]*temp[j]; 
	*vptr[i]=s;
      }
    }
    
   void applyToT(SCALAR** vptr) const{
      SCALAR temp[k];
      for(int i=0; i<k; i++) temp[i]=*vptr[i];
      for(int i=0; i<k; i++){
	double s=0; 
	for(int j=0; j<k; j++) 
	  s+=q[i*k+j]*temp[j]; 
	*vptr[i]=s;
      }
    }
    

  public: // I/O
    
    string str() const{
      ostringstream result;
      result<<"k-point operator on "<<map.str()<<":"<<endl<<as_Cmatrix(k,k,q)<<endl;
      return result.str();
    }

  };


  template<int k>
  inline ostream& operator<<(ostream& stream, const KpointOp<k>& x){stream<<x.str(); return stream;}
  

}

#endif
