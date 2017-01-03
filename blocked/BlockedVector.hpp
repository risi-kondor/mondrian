#ifndef _BlockedVector
#define _BlockedVector

#include "Vector.hpp"
#include "BlockStructure.hpp"


namespace Mondrian{
  
  template<class VECTOR>
  class BlockedVector{ //: public Detachable{ // inherit from Vector?
  public:
    
    BlockedVector(): n(0), nb(0){};
    BlockedVector(const int _nb, const int _n=0);
    BlockedVector(const BlockStructure& b);
    ~BlockedVector();

    int n; // inherit from Vector?
    int nb;
    VECTOR** blocks=nullptr;

    
  public: //copying

    BlockedVector(const BlockedVector<VECTOR>& x);
    BlockedVector(BlockedVector<VECTOR>&& x);
    BlockedVector& operator=(const BlockedVector<VECTOR>& x);
    BlockedVector& operator=(BlockedVector<VECTOR>&& x);
    BlockedVector<VECTOR> copy() const;

        
  public: // named constructors
    
    static BlockedVector<VECTOR> Zero(const BlockStructure& b);
    static BlockedVector<VECTOR> Filled(const BlockStructure& b, const SCALAR t);
    static BlockedVector<VECTOR> Gaussian(const BlockStructure& b);


  public: // block access 
    
    VECTOR& block(const int i) {assert(i<nb); return *blocks[i];}
    const VECTOR& block(const int i) const {assert(i<nb); return *blocks[i];}

    void for_each_block(std::function<void(const INDEX, VECTOR&)> lambda){
      for(int b=0; b<nb; b++) lambda(b,*blocks[b]);}
    void for_each_block(std::function<void(const INDEX, const VECTOR&)> lambda) const{
      for(int b=0; b<nb; b++) lambda(b,*blocks[b]);}

    VECTOR*& block_ptr(const int i) {assert(i<nb); return blocks[i];}
    const VECTOR* block_ptr(const int i) const {assert(i<nb); return blocks[i];}


  public: // element access

    SCALAR& operator()(const int i) {assert(i<n);
      int t=0; int j=0; for(; j<nb && t+blocks[j]->n<i; j++) t+=blocks[j]->n;  
      return (*blocks[j])(i-t);}
    SCALAR operator()(const int i) const {assert(i<n);
      int t=0; int j=0; for(; j<nb && t+blocks[j]->n<i; j++) t+=blocks[j]->n; 
      return (*blocks[j])(i-t);}    

    void for_each(std::function<void(const INDEX, SCALAR&)> lambda) {
      int t=0; for(int b=0; b<nb; b++) for(int i=0; i<blocks[b]->n; i++) lambda(t++,(*blocks[b])(i));}
    void for_each(std::function<void(const INDEX, const SCALAR)> lambda) const {
      int t=0; for(int b=0; b<nb; b++) for(int i=0; i<blocks[b]->n; i++) lambda(t++,(*blocks[b])(i));}
    

  public: // in-place arithmetic
    
    BlockedVector<VECTOR>& operator+=(const SCALAR x);
    BlockedVector<VECTOR>& operator-=(const SCALAR x);
    BlockedVector<VECTOR>& operator*=(const SCALAR x);
    BlockedVector<VECTOR>& operator/=(const SCALAR x);

    template<class VECTOR2>  
    BlockedVector<VECTOR>& operator+=(const BlockedVector<VECTOR2>& x);
    template<class VECTOR2>  
    BlockedVector<VECTOR>& operator-=(const BlockedVector<VECTOR2>& x);


  public: // arithmetic
    
    template<class VECTOR2>
    auto operator+(const BlockedVector<VECTOR2>& x) const -> decltype((*blocks[0])+x);
    template<class VECTOR2>
    auto operator-(const BlockedVector<VECTOR2>& x) const -> decltype((*blocks[0])-x); 

    template<class VECTOR2>
    SCALAR dot(const BlockedVector<VECTOR2>& x) const;  


  public:
    
    string str() const;

  };


  // ---- Constructors ---------------------------------------------------------------------------------------------


  template<class VECTOR>
  BlockedVector<VECTOR>::BlockedVector(const int _nb, const int _n): n(_n), nb(_nb){
    blocks=new VECTOR*[nb];}

  template<class VECTOR>
  BlockedVector<VECTOR>::BlockedVector(const BlockStructure& B): BlockedVector(B.size(),0){
    for(int i=0; i<nb; i++){blocks[i]=new VECTOR(B[i]); n+=B[i];}
  }  

  
  // ---- Copying --------------------------------------------------------------------------------------------------


  template<class VECTOR>
  BlockedVector<VECTOR>::BlockedVector(const BlockedVector<VECTOR>& X): BlockedVector(X.nb,X.n){
    MultiLoop(nb,[this,&X](const int i){blocks[i]=new VECTOR(*X.blocks[i]);});}
  
  template<class VECTOR>
  BlockedVector<VECTOR>::BlockedVector(BlockedVector<VECTOR>&& X): BlockedVector(X.nb,X.n){
    for(int i=0; i<nb; i++){blocks[i]=X.blocks[i]; X.blocks[i]=nullptr;}; X.nb=0; X.n=0;}

  template<class VECTOR>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator=(const BlockedVector<VECTOR>& X){
    for(int i=0; i<nb; i++) delete blocks[i]; delete[] blocks;
    nb=X.nb; blocks=new VECTOR*[nb];
    MultiLoop(nb,[this,X](const int i){blocks[i]=new VECTOR(*X.blocks[i]);});
    return *this;
  }

  template<class VECTOR>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator=(BlockedVector<VECTOR>&& X){
    for(int i=0; i<nb; i++) delete blocks[i]; delete[] blocks;
    nb=X.nb; blocks=new VECTOR*[nb]; X.nb=0; X.n=0;
    for(int i=0; i<nb; i++){blocks[i]=X.blocks[i]; X.blocks[i]=nullptr;};
    return *this;
  }
  
  template<class VECTOR>
  BlockedVector<VECTOR> BlockedVector<VECTOR>::copy() const{
    BlockedVector M(nb,n);
    MultiLoop(nb,[this,&M](const int i){M.blocks[i]=new VECTOR(blocks[i]->copy());});
    return M;
  }

  template<class VECTOR>
  BlockedVector<VECTOR>::~BlockedVector(){
    for(int i=0; i<nb; i++) delete blocks[i]; delete blocks;
  }


  // ---- Named Constructors ---------------------------------------------------------------------------------------
  
  
  template<class VECTOR>
  BlockedVector<VECTOR> BlockedVector<VECTOR>::Zero(const BlockStructure& B){
    BlockedVector<VECTOR> v(B.v.size(),0);
    for(int i=0; i<B.size(); i++) {v.blocks[i]=new VECTOR(VECTOR::Zero(B[i])); v.n+=B[i];}
    return v;
  }

  template<class VECTOR>
    BlockedVector<VECTOR> BlockedVector<VECTOR>::Filled(const BlockStructure& B, const SCALAR t){
    BlockedVector<VECTOR> v(B.v.size(),0);
    for(int i=0; i<B.size(); i++) {v.blocks[i]=new VECTOR(VECTOR::Filled(B[i],t)); v.n+=B[i];}
    return v;
  }

  template<class VECTOR>
  BlockedVector<VECTOR> BlockedVector<VECTOR>::Gaussian(const BlockStructure& B){
    BlockedVector<VECTOR> v(B.v.size(),0);
    for(int i=0; i<B.size(); i++) {v.blocks[i]=new VECTOR(VECTOR::Gaussian(B[i])); v.n+=B[i];}
    return v;
  }


  // ---- In-place Arithmetic --------------------------------------------------------------------------------------


  template<class VECTOR>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator+=(const SCALAR x){
    for(int i=0; i<nb; i++) (*blocks[i])+=x; return *this;}

  template<class VECTOR>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator-=(const SCALAR x){
    for(int i=0; i<nb; i++) (*blocks[i])-=x; return *this;}

  template<class VECTOR>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator*=(const SCALAR x){
    for(int i=0; i<nb; i++) (*blocks[i])*=x; return *this;}

  template<class VECTOR>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator/=(const SCALAR x){
    for(int i=0; i<nb; i++) (*blocks[i])/=x; return *this;}

  template<class VECTOR>
  template<class VECTOR2>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator+=(const BlockedVector<VECTOR2>& x){
    assert(x.nb==nb); for(int i=0; i<nb; i++) (*blocks[i])+=(*x.blocks[i]); return *this;}
  
  template<class VECTOR>
  template<class VECTOR2>
  BlockedVector<VECTOR>& BlockedVector<VECTOR>::operator-=(const BlockedVector<VECTOR2>& x){
    assert(x.nb==nb);
    for(int i=0; i<nb; i++) (*blocks[i])-=(*x.blocks[i]); return *this;}


  // ---- Arithmetic -----------------------------------------------------------------------------------------------


  template<class VECTOR>
  template<class VECTOR2>
  auto BlockedVector<VECTOR>::operator+(const BlockedVector<VECTOR2>& x) const ->decltype((*blocks[0])+x){
    assert(x.nb==nb); assert(x.n==n);
    BlockedVector<decltype((*blocks[0])+x)> r(nb,n);
    for(int i=0; i<nb; i++) 
      r.blocks[i]= new auto((*blocks[i])+(*x.blocks[i]));
    return r;
  }

  template<class VECTOR>
  template<class VECTOR2>
  auto BlockedVector<VECTOR>::operator-(const BlockedVector<VECTOR2>& x) const ->decltype((*blocks[0])-x){
    assert(x.nb==nb); assert(x.n==n);
    BlockedVector<decltype((*blocks[0])-x)> r(nb,n);
    for(int i=0; i<nb; i++) 
      r.blocks[i]= new auto((*blocks[i])-(*x.blocks[i]));
    return r;
  }

  template<class VECTOR>
  template<class VECTOR2>
  SCALAR BlockedVector<VECTOR>::dot(const BlockedVector<VECTOR2>& x) const{
    assert(x.nb==nb);
    SCALAR r=0; 
    for(int i=0; i<nb; i++) r+=blocks[i]->dot(*blocks[i]);
    return r;
  }


  // ---- I/O ------------------------------------------------------------------------------------------------------


  template<class VECTOR>
  string BlockedVector<VECTOR>::str() const{
    ostringstream oss; oss<<"{";
    for(int i=0; i<nb-1; i++) oss<<blocks[i]->str()<<",";
    if(nb>0) oss<<blocks[nb-1]->str(); oss<<"}";
    return oss.str();
  }
  
  
  // ---- Global ---------------------------------------------------------------------------------------------------


  template<class VECTOR>
  BlockedVector<VECTOR> operator*(const SCALAR c, const BlockedVector<VECTOR>& v){return v*c;}

  

} // namespace Mondrian

#endif
