#ifndef _AccumulatedVector
#define _AccumulatedVector

#include "ActiveVector.hpp"
#include "AccumulationTree.hpp"

namespace Mondrian{

  template<class VECTOR>
  class AccumulatedVector: public ActiveVector<VECTOR>{
  public:
    
    AccumulationTree<ivpair> tree;

    using ActiveVector<VECTOR>::n;
    using ActiveVector<VECTOR>::read;
    using ActiveVector<VECTOR>::set;


  public:
    
    using ActiveVector<VECTOR>::ActiveVector;

    AccumulatedVector(const initializer_list<SCALAR> list): 
      ActiveVector<VECTOR>(list){
      makeTree();
    }

    AccumulatedVector(const int _n, const initializer_list<ivpair> list): 
      ActiveVector<VECTOR>(_n){
      makeTree();
    }


  public: // change function ----------------------------------------------------------------------------------------
    

    void makeTree(){
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){
	  tree.addLeaf(i,ivpair(i,v));});
    }

    virtual void changed(const INDEX i, const SCALAR& v){
      tree.set(i,ivpair(i,v));
    }


  public: // copying ------------------------------------------------------------------------------------------------


    explicit AccumulatedVector(const AccumulatedVector<VECTOR>& x): 
      ActiveVector<VECTOR>(x,_NoWarn()),
      tree(x.tree){
      COPY_WARNING("AccumulatedVector<VECTOR>");
    } 

    explicit AccumulatedVector(const AccumulatedVector<VECTOR>& x, const _NoWarn dummy): 
      ActiveVector<VECTOR>(x,_NoWarn()),
      tree(x.tree){
    } 

    explicit AccumulatedVector(AccumulatedVector<VECTOR>&& x): 
      ActiveVector<VECTOR>(std::move(x)),
      tree(std::move(x.tree)){
      MOVE_WARNING("AccumulatedVector<VECTOR>");
    }

    AccumulatedVector& operator=(const AccumulatedVector<VECTOR>& x){
      VECTOR::operator=(x);
      tree=x.tree;
      return *this;
    }
    
    AccumulatedVector& operator=(AccumulatedVector<VECTOR>&& x){
      VECTOR::operator=(std::move(x));
      tree=std::move(x.tree);
      return *this;
    }
    
    AccumulatedVector<VECTOR> copy() const {
      AccumulatedVector<VECTOR> r(VECTOR::copy());
      r.tree=tree;
      return r;
    }
    
    AccumulatedVector<VECTOR> shallow() const { // check!
      AccumulatedVector<VECTOR> r(VECTOR::shallow());
      r.tree=tree;
    }
    
    void assign(const AccumulatedVector<VECTOR>& x) const { // check! 
      VECTOR::assign(x);
      tree=x.tree;
    }
    
    void detach(){VECTOR::detach();}


  public: // downcasting --------------------------------------------------------------------------------------------


    AccumulatedVector(const VECTOR& x):
      ActiveVector<VECTOR>(x){
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){tree.set(i,v);});
    }

    AccumulatedVector(VECTOR&& x):
      ActiveVector<VECTOR>(std::move(x)){
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){tree.set(i,v);});
    }

    AccumulatedVector& operator=(const VECTOR& x){
      VECTOR::operator=(std::move(x));
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){tree.set(i,v);});
      return *this;
    }
    
    AccumulatedVector& operator=(VECTOR&& x){
      VECTOR::operator=(std::move(x));
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){tree.set(i,v);});
      return *this;
    }
    

  public: // access -------------------------------------------------------------------------------------------------


    SCALAR sum() const{  // const added
      return tree.root->decor.second;
    }

 
  public: // I/O ----------------------------------------------------------------------------------------------------


    string str() const{
      return ActiveVector<VECTOR>::str();
 
	//ostringstream oss;
	//return oss.str();
    }


  };

  template<class VECTOR>
  inline ostream& operator<<(ostream& stream, const AccumulatedVector<VECTOR>& v){stream<<v.str(); return stream;}


}


#endif
