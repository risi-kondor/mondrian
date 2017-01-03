#ifndef _SparseVector
#define _SparseVector

#include "Vector.hpp"

namespace Mondrian{

class SparseVector: public Vector{
public:

  SparseVector(const int _n): Vector(_n){}


public: // attributes

  virtual bool isSparse() const {return true;}


public: // sparse vector methods

  virtual SCALAR* findptr(const int i)=0;
  virtual void insert(const int i, const SCALAR value)=0;
  virtual void append(const int i, const SCALAR value)=0;
  //virtual void zero(const int i)=0;
  virtual void sort() const {}
  virtual void tidy()=0;



};

} // Mondrian

#endif
