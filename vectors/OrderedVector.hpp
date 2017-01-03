#ifndef _OrderedVector
#define _OrderedVector

#include "ActiveVector.hpp"
#include "OrderedSet.hpp"

namespace Mondrian{

  template<class VECTOR>
  class OrderedVector: public ActiveVector<VECTOR>{
  public:
    
    OrderedSet<> ordering;

    using ActiveVector<VECTOR>::n;
    using ActiveVector<VECTOR>::read;
    using ActiveVector<VECTOR>::set;


  public:
    
    using ActiveVector<VECTOR>::ActiveVector;

    OrderedVector(const initializer_list<SCALAR> list): 
      ActiveVector<VECTOR>(list.size()){
      int i=0; for(SCALAR v:list) set(i++,v);}

    OrderedVector(const int _n, const initializer_list<ivpair> list): 
      ActiveVector<VECTOR>(_n){
      for(const ivpair& v:list) {assert(v.first<n); set(v.first)=v.second;}
    }


  public: // change function ----------------------------------------------------------------------------------------
    

    virtual void changed(const INDEX i, const SCALAR& v){
      ordering.set(i,v);
      //cout<<"v("<<i<<")<-"<<v<<endl;
    }


  public: // copying ------------------------------------------------------------------------------------------------


    explicit OrderedVector(const OrderedVector<VECTOR>& x): 
      ActiveVector<VECTOR>(x,_NoWarn()),
      ordering(x.ordering){
      COPY_WARNING("OrderedVector<VECTOR>");
    } 

    explicit OrderedVector(const OrderedVector<VECTOR>& x, const _NoWarn dummy): 
      ActiveVector<VECTOR>(x,_NoWarn()),
      ordering(x.ordering){
    } 

    explicit OrderedVector(OrderedVector<VECTOR>&& x): 
      ActiveVector<VECTOR>(std::move(x)),
      ordering(std::move(x.ordering)){
      MOVE_WARNING("OrderedVector<VECTOR>");
    }

    OrderedVector& operator=(const OrderedVector<VECTOR>& x){
      VECTOR::operator=(x);
      ordering=x.ordering;
      return *this;
    }
    
    OrderedVector& operator=(OrderedVector<VECTOR>&& x){
      VECTOR::operator=(std::move(x));
      ordering=std::move(x.ordering);
      return *this;
    }
    
    OrderedVector<VECTOR> copy() const {
      OrderedVector<VECTOR> r(VECTOR::copy());
      r.ordering=ordering;
      return r;
    }
    
    OrderedVector<VECTOR> shallow() const { // check!
      OrderedVector<VECTOR> r(VECTOR::shallow());
      r.ordering=ordering;
    }
    
    void assign(const OrderedVector<VECTOR>& x) const { // check! 
      VECTOR::assign(x);
      ordering=x.ordering;
    }
    
    void detach(){VECTOR::detach();}


  public: // downcasting --------------------------------------------------------------------------------------------


    OrderedVector(const VECTOR& x):
      ActiveVector<VECTOR>(x){
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){ordering.set(i,v);});
    }

    OrderedVector(VECTOR&& x):
      ActiveVector<VECTOR>(std::move(x)){
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){ordering.set(i,v);});
    }

    OrderedVector& operator=(const VECTOR& x){
      VECTOR::operator=(std::move(x));
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){ordering.set(i,v);});
      return *this;
    }
    
    OrderedVector& operator=(VECTOR&& x){
      VECTOR::operator=(std::move(x));
      static_cast<const VECTOR&>(*this).for_each_filled([this](const INDEX i, const SCALAR v){ordering.set(i,v);});
      return *this;
    }
    

  public: // access -------------------------------------------------------------------------------------------------


    INDEX best(const int c=0){
      return ordering.best(c);
    }


    INDEX worst(const int c=0){
      return ordering.worst(c);
    }


  public: // I/O ----------------------------------------------------------------------------------------------------


    string str() const{
      ostringstream oss;
      for(auto& p:ordering.ordered){
	int i=p.first;
	oss<<"("<<i<<","<<read(i)<<")"<<endl;
      }
      return oss.str();
    }


  };

  template<class VECTOR>
  inline ostream& operator<<(ostream& stream, const OrderedVector<VECTOR>& v){stream<<v.str(); return stream;}


}


#endif
