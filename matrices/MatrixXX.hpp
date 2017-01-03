#ifndef _MatrixXX
#define _MatrixXX

#include "MatrixX.hpp"

extern default_random_engine randomNumberGenerator;

namespace Mondrian{


template<class VECTOR>
class MatrixXX: public MatrixX<VECTOR>{
public:

  MatrixX<VECTOR> tr;

  using MatrixX<VECTOR>::cols;
  using MatrixX<VECTOR>::nrows;
  using MatrixX<VECTOR>::ncols;


public: // constructors

  MatrixXX(): MatrixXX(0,0){};

  MatrixXX(const int _nrows, const int _ncols): 
    MatrixX<VECTOR>(_nrows,_ncols), 
    tr(MatrixX<VECTOR>(_ncols,_nrows)){}
  

  MatrixXX(const int _nrows, const int _ncols, const _Zero dummy): 
    MatrixX<VECTOR>(_nrows,_ncols,dummy),
    tr(MatrixX<VECTOR>(_ncols,_nrows,dummy)){}

  ~MatrixXX(){}
  

public: // copying

  MatrixXX(const MatrixXX<VECTOR>& x): 
    MatrixX<VECTOR>(x,_NoWarn()){
    COPY_WARNING("MatrixXX<VECTOR>");
    tr=MatrixX<VECTOR>(x.tr,_NoWarn());
  }

  MatrixXX(const MatrixXX<VECTOR>& x, const _NoWarn dummy): 
    MatrixX<VECTOR>(x,_NoWarn()){
    tr=MatrixX<VECTOR>(x.tr,_NoWarn());
  }

  MatrixXX(MatrixXX<VECTOR>&& x): 
    MatrixX<VECTOR>(std::move(x)){
    MOVE_WARNING("MatrixXX<VECTOR>");
    //tr=std::move(x.tr);
  }

  MatrixXX(MatrixXX<VECTOR>&& x, const _NoWarn dummy): 
    MatrixX<VECTOR>(std::move(x)){
    tr=std::move(x.tr);
  }

  MatrixXX<VECTOR>& operator=(const MatrixXX<VECTOR>& x){
    ASSIGN_WARNING("MatrixXX<VECTOR>");
    MatrixXX<VECTOR>::operator=(x);
    tr=x.tr;
  }

  MatrixXX<VECTOR>& operator=(MatrixXX<VECTOR>&& x){
    MOVEASSIGN_WARNING("MatrixXX<VECTOR>");
    nrows=x.nrows; ncols=x.ncols;
    for(auto p:cols) delete p; 
    cols=x.cols; 
    x.nrows=0; x.ncols=0;
    tr=std::move(x);
    return *this;
  }

  MatrixXX<VECTOR> copy() {return MatrixXX(*this,_NoWarn());}

  void detach(){cols.clear(); tr.detach();}

  void assign(const MatrixXX& x){
    assert(x.nrows==nrows);
    assert(x.ncols==ncols);
    MatrixXX::assign(x);
    tr.assign(x.tr);
  }

  MatrixXX shallow() const{
    MatrixXX M(nrows,ncols); 
    M.cols=cols;
    M.tr=tr.shallow();
    return M;
  }


public: // named constructors

  static MatrixXX<VECTOR> Zero(const int _nrows, const int _ncols){
    return MatrixX<VECTOR>::Zero(_nrows,_ncols);}

  static MatrixXX<VECTOR> Filled(const int _nrows, const int _ncols, const SCALAR v){
    return MatrixX<VECTOR>::Filled(_nrows,_ncols,v);}

  static MatrixXX<VECTOR> Identity(const int n){
    return MatrixX<VECTOR>::Identity(n);}

  static MatrixXX<VECTOR> Uniform(const int _nrows, const int _ncols){
    return MatrixX<VECTOR>::Uniform(_nrows,_ncols);}

  static MatrixXX<VECTOR> Gaussian(const int _nrows, const int _ncols){
    //return MatrixXX<VECTOR>(0,0);}
    return MatrixX<VECTOR>::Gaussian(_nrows,_ncols);}

  static MatrixXX<VECTOR> Bernoulli(const int _nrows, const int _ncols, const double p=0.5){
    return MatrixX<VECTOR>::Bernoulli(_nrows,_ncols,p);}


public: // conversions 

  MatrixXX(const MatrixX<VECTOR>& x):
    MatrixX<VECTOR>(x),
    tr(MatrixTranspose<MatrixX<VECTOR> >(x)){}

  MatrixXX(MatrixX<VECTOR>&& x):
    MatrixX<VECTOR>(std::move(x)){
    tr=MatrixTranspose<MatrixX<VECTOR> >(x);
  }

  /*
  template<class VECTOR2>
  MatrixXX(const MatrixXX<VECTOR2>& x): 
    MatrixX<VECTOR>(x){
    CONVERT_WARNING("MatrixXX<VECTOR2>","MatrixXX<VECTOR>");
    tr=x.tr;
  }
  */

  MatrixXX(const Cmatrix& x):
    MatrixX<VECTOR>(x){
    CONVERT_WARNING("Cmatrix","MatrixXX<VECTOR>");
    mirror();
  }

  MatrixXX(const Cmatrix& x, const _NoWarn dummy):
    MatrixX<VECTOR>(x){
    mirror();
  }



public: // element access
  
  using MatrixX<VECTOR>::operator();
  using MatrixX<VECTOR>::read;
  using MatrixX<VECTOR>::isFilled;
  using MatrixX<VECTOR>::nFilled;
  using MatrixX<VECTOR>::column;

  SCALAR& operator()(const int i, const int j){
    cout<<"Warning: MatrixXX::operator()->SCALAR& should not be used"<<endl;
    assert(i<nrows); assert(j<ncols); return (*cols[j])(i);}

  void set(const int i, const int j, const SCALAR v){ 
    assert(i<nrows); assert(j<ncols); (*cols[j])(i)=v; (*tr.cols[i])(j)=v;} 

  template<class VECTOR2>
  VECTOR2 row(const int i) const {assert(i<nrows);
    return VECTOR2(tr->cols[i]);}


public: // iterators

  void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
    cout<<"Warning: MatrixXX::for_each_filled(...)->SCALAR& should not be used"<<endl;
    for(int j=0; j<ncols; j++) cols[j]->for_each_filled([j,&lambda](const int i, SCALAR& v){lambda(i,j,v);});}

  void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR)> lambda) const{
    for(int j=0; j<ncols; j++) cols[j]->for_each_filled([j,&lambda](const int i, const SCALAR v){lambda(i,j,v);});}

  void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
    cout<<"Warning: MatrixXX::for_each_filled_in_row(...)->SCALAR& should not be used"<<endl;
    tr.cols[i]->for_each_filled(lambda);}

  void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR)> lambda) const{
    static_cast<const VECTOR*>(tr.cols[i])->for_each_filled(lambda);}

  void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
    cout<<"Warning: MatrixXX::for_each_filled_in_column(...)->SCALAR& should not be used"<<endl;
    cols[j]->for_each_filled(lambda);}

  void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR)> lambda) const{
    static_cast<const VECTOR*>(cols[j])->for_each_filled(lambda);}


public: // sparse matrix methods

  void sort() {MatrixX<VECTOR>::sort(); tr.sort();}

  VectorSpattern nonemptyRows() const{
    return VectorSpattern(cols);
  }

  VectorSpattern nonemptyColumns() const{
    return VectorSpattern(tr.cols);
  }

  void mirror(){
    tr=MatrixTranspose<MatrixX<VECTOR> >(*this);
  }


public: // in-place operations

  void symmetrize(){
    MatrixX<VECTOR>::symmetrize();
    tr=MatrixX<VECTOR>::copy();
  }


public: // scalar valued operations

  using MatrixX<VECTOR>::nnz;
  using MatrixX<VECTOR>::norm2;
  using MatrixX<VECTOR>::diff2;
  using MatrixX<VECTOR>::spectralNorm;
  using MatrixX<VECTOR>::inp_of_columns;

  SCALAR inp_of_rows(const int j1, const int j2) const{
    return tr.cols[j1]->dot(*tr.cols[j2]);}


public: // vector valued operations 

  using MatrixX<VECTOR>::dot;

  template<class VECTOR2>
  VECTOR2 operator*(const VECTOR2& x) const{
    VECTOR2 r = Cvector::Zero(nrows);
    for(int i=0; i<nrows; i++)
      r(i)=tr.vecs[i]->dot(x);
    return r;
  }


public: // in place operations

  MatrixXX<VECTOR>& operator+=(const MatrixXX<VECTOR>& x){ 
    assert(nrows==x.nrows); assert(ncols==x.ncols);
    MatrixX<VECTOR>::operator+=(x);
    tr+=x.tr;
    return *this;
  }

  MatrixXX<VECTOR>& operator-=(const MatrixXX<VECTOR>& x){ 
    assert(nrows==x.nrows); assert(ncols==x.ncols);
    MatrixX<VECTOR>::operator-=(x);
    tr-=x.tr;
    return *this;
  }

  template<class VECTOR2>
  MatrixX<VECTOR>& multiplyRowsBy(const VECTOR2& v){
    MatrixX<VECTOR>::multiplyRowsBy(v);
    tr.multiplyColumnsBy(v);
  }

  template<class VECTOR2>
  MatrixX<VECTOR>& multiplyColumnsBy(const VECTOR2& v){
    MatrixX<VECTOR>::multiplyColumnsBy(v);
    tr.multiplyRowsBy(v);
  }

  template<class VECTOR2>
  MatrixX<VECTOR>& divideRowsBy(const VECTOR2& v){
    MatrixX<VECTOR>::divideRowsBy(v);
    tr.divideColumnsBy(v);
  }

  template<class VECTOR2>
  MatrixX<VECTOR>& divideColumnsBy(const VECTOR2& v){
    MatrixX<VECTOR>::divideColumnsBy(v);
    tr.divideRowsBy(v);
  }

  MatrixX<VECTOR>& transpose(){
    MatrixX<VECTOR> t=std::move(tr);
    tr=std::move(*this);
    MatrixX<VECTOR>::operator=(std::move(t));
    return *this;
  }

  MatrixX<VECTOR>& tidy(){
    MatrixX<VECTOR>::tidy();
    tr.tidy();
  }






};

} // namespace Mondrian 

#endif
