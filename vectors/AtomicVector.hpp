#ifndef _AtomicVector
#define _AtomicVector

namespace Mondrian{


  template<class VECTOR>
  class AtomicVector<VECTOR>: public VECTOR{
  public:
    
    using VECTOR::VECTOR;

    mutex mx;
    

  public: // element access

    void for_each(std::function<void(const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx);
      VECTOR::for_each(lambda);
    }


  public: // in-place arithmetic

    VECTOR& operator+=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator+=(x);
      return *this;
    }

    VECTOR& operator-=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator-=(x);
      return *this;
    }

    VECTOR& operator*=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator+=(x);
      return *this;
    }

    VECTOR& operator/=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator-=(x);
      return *this;
    }

    VECTOR& operator+=(const VECTOR& x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator+=(x);
      return *this;
    }

    VECTOR& operator-=(const VECTOR& x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator-=(x);
      return *this;
    }


  public: // in-place operators

    template<class OPERATOR>
    void apply(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      VECTOR::apply(op);
    }

    template<class OPERATOR>
    void applyInv(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      VECTOR::apply(op);
    }


  };

}

#endif
