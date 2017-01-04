#ifndef _Mondrian_base
#define _Mondrian_base

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <random>
#include <functional>
#include <algorithm>

#include <math.h>
#include <assert.h>

using namespace std;


#define COPY_WARNING(name) cout<<"Warning: "<<#name<<" copied."<<endl

//#define ASSIGN_WARNING(name)
#define ASSIGN_WARNING(name) cout<<"Warning: "<<#name<<" assigned."<<endl

#define MOVE_WARNING(name) 
//#define MOVE_WARNING(name) cout<<"Warning: "<<#name<<" moved."<<endl

#define MOVEASSIGN_WARNING(name)
//#define MOVEASSIGN_WARNING(name) cout<<"Warning: "<<#name<<" move-assigned."<<endl

#define TRANSPOSE_WARNING(name)
//#define TRANSPOSE_WARNING(name) cout<<"Warning: "<<#name<<" transposed."<<endl

//#define CONVERT_WARNING(name1,name2)
#define CONVERT_WARNING(name1,name2) cout<<"Warning: "<<#name1<<" copy converted to "<<#name2<<"."<<endl;

#define MOVECONVERT_WARNING(name1,name2) cout<<"Warning: "<<#name1<<" move converted to "<<#name2<<"."<<endl;

//#define COPYINTO_WARNING(name1,name2) cout<<"Warning: "<<#name1<<" copied into "<<#name2<<"."<<endl;
//#define MOVEINTO_WARNING(name1,name2) cout<<"Warning: "<<#name1<<" moved into "<<#name2<<"."<<endl;

#define DOWNCAST_WARNING(name1,name2)
//#define DOWNCAST_WARNING(name1,name2) cout<<"Message: "<<#name1<<" downcast to "<<#name2<<"."<<endl;

#define DOWNCASTCOPY_WARNING(name1,name2) cout<<"Warning: "<<#name1<<" downcast to "<<#name2<<" by copying."<<endl;

#define SYMMETRY_UNSAFE(name1) cout<<"Warning: function "<<#name1<<" is not symmetry-safe."<<endl;
#define SYMMETRY_FORBIDDEN(name1) cout<<"Warning: "<<#name1<<" not implemented for symmetric matrices."<<endl; 

#define UNIMPL_WARNING(name1) cout<<"Warning: "<<#name1<<" not implemented."<<endl;



typedef int    INDEX;
typedef double SCALAR;


#include "Serializable.hpp"
#include "Inverse.hpp"


namespace Mondrian{

  // vectors

  class Cvector;
  class Bvector;
  class Vectorv;
  class Vectorl;
  class Vectorh;
  class LazyVectorv;

  class Cvectorm;

  template<class VECTOR> 
  class VectorView;


  // matrices

  class Cmatrix;
  class Cmatrixm;
  class CmatrixLA;
  class SymmCmatrix;
  class Bmatrix;

  template<class VECTOR> 
  class MatrixX;
  template<class VECTOR> 
  class SymmMatrixX;

  template<class VECTOR> 
  class MatrixXm;
  template<class VECTOR> 
  class SymmMatrixXm;

  template<class MATRIX> 
  class MatrixView;

  
  // index maps

  class IndexMap;
  class IndexBiMap;
  class BindexMap;
  class BtoBindexMap;
  class BtoBindexBiMap;


  // operators 

  class GivensRotation;
  class KpointOperator;

  template<int k>
  class KpointOp;


  // markers 

  class _Zero{};
  class _Uninitialized{};
  class _NoWarn{};
  class _Downcast{};
  //  class _ConstDowncast{};

  class _Sparse{};

}

  class EigenVectorXdAdaptor;
  class EigenMatrixXdAdaptor;



class CoutLock{
public:
  CoutLock(): lock(mx){}
  lock_guard<mutex> lock;
  static mutex mx;
};


template<typename TYPE>
inline void swapp(TYPE& x, TYPE& y) {TYPE t=x; x=y; y=t;}

template<class TYPE>
TYPE snatch(TYPE& x){TYPE t=x; x=0; return t;}

template<class TYPE>
TYPE* snatchptr(TYPE*& x){TYPE* t=x; x=nullptr; return t;}

template<typename TYPE>
inline void move_over(TYPE*& ptr, TYPE*& result) {delete ptr; ptr=result; result=nullptr;}

template<class TYPE>
TYPE* temp2ptr(TYPE&& x){return new TYPE(x);}

template<typename TYPE>
inline void replace(TYPE*& ptr, TYPE* result) {delete ptr; ptr=result;}

template<class TYPE>
void debug(TYPE v){cout<<v<<endl;}

//template<class TYPE>
//copy_to_array(TYPE::iterator _beg, TYPE::iterator:: _end)



typedef pair<int,int> iipair;

struct ivpair{
  ivpair()=default;
  ivpair(const INDEX& _first, const SCALAR& _second):first(_first),second(_second){}
  bool operator==(const ivpair& x) const {return (first==x.first)&&(second==x.second);}
  ivpair operator+(const ivpair& x) const {return ivpair(0,second+x.second);}
  string str() const {ostringstream result; result<<"("<<first<<","<<second<<")"; return result.str();}  
  INDEX first; 
  SCALAR second;
  //struct larger;
  //struct{
  //  bool operator()(const ivpair& a, const ivpair& b) const {return a.second>b.second;}} larger;
  static struct{
    bool operator()(const ivpair& a, const ivpair& b) const {return fabs(a.second)>fabs(b.second);}} larger_abs;
};

struct{
  bool operator()(const ivpair& a, const ivpair& b) const {return a.second>b.second;}
}ivlarger;
struct{
  bool operator()(const ivpair& a, const ivpair& b) const {return fabs(a.second)>fabs(b.second);}
}ivlarger_abs;
struct{
  bool operator()(const ivpair& a, const ivpair& b) const {return a.second>b.second;}
}ivsmaller;
struct{
  bool operator()(const ivpair& a, const ivpair& b) const {return fabs(a.second)>fabs(b.second);}
}ivsmaller_abs;


struct iivtriple{
  iivtriple()=default;
  iivtriple(const INDEX& _first, const INDEX& _second, const SCALAR& _third): 
    first(_first),second(_second), third(_third){}
  bool operator==(const iivtriple& x) const {return (first==x.first)&&(second==x.second)&&(third==x.third);}
  string str() const {ostringstream result; result<<"("<<first<<","<<second<<","<<third<<")"; return result.str();}  
  INDEX first, second; 
  SCALAR third;
};

struct{
  bool operator()(const ivpair& a, const ivpair& b){return a.first<b.first;} 
}ivpairComparator;




#include "package.hpp"


#endif

