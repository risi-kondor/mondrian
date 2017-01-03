#ifndef _AtomicMatrix
#define _AtomicMatrix

namespace Mondrian{


  template<class MATRIX>
  class AtomicMatrix<MATRIX>: public MATRIX{
  public:
    
    using MATRIX::MATRIX;

    mutex mx;

    
  public: // element access

    void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx); VECTOR::for_each(lambda);}
    void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx); VECTOR::for_eachin_row(i,lambda);}
    void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx); VECTOR::for_each_in_column(j,lambda);}


  public: // in-place arithmetic 

    MATRIX& operator+=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator+=(x);
      return *this;
    }

    MATRIX& operator-=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator-=(x);
      return *this;
    }

    MATRIX& operator*=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator+=(x);
      return *this;
    }

    MATRIX& operator/=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator-=(x);
      return *this;
    }

    MATRIX& operator+=(const MATRIX& x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator+=(x);
      return *this;
    }

    MATRIX& operator-=(const MATRIX& x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator-=(x);
      return *this;
    }


  public: // in-place operations

    void multiplyRowsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::multiplyRowsBy(v);}
    void divideRowsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::divideRowsBy(v);}
    void multiplyColsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::multiplyColsBy(v);}
    void divideColsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::divideColsBy(v);}
      
    template<class OPERATOR>
    void apply(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      MATRIX::apply(op);
    }

    template<class OPERATOR>
    void applyInv(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      MATRIX::apply(op);
    }


  };

}

#endif
