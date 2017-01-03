#ifndef _MatrixX
#define _MatrixX

#include "SparseMatrix.hpp"
#include "Cmatrix.hpp"
#include "GivensRotation.hpp"
#include "VectorSpattern.hpp"

extern default_random_engine randomNumberGenerator;

namespace Mondrian{


  template<class VECTOR>
  class MatrixX: public Mondrian::SparseMatrix{
  public:

    vector<VECTOR*> cols; 
    

  public: // constructors -------------------------------------------------------------------------------------------

    
    MatrixX(): MatrixX(0,0){};
    
    MatrixX(const int _nrows, const int _ncols): 
      SparseMatrix(_nrows,_ncols), cols(ncols){
      for(int i=0; i<ncols; i++) cols[i]=new VECTOR(nrows);
    }
    
    MatrixX(const int _nrows, const int _ncols, const _Zero dummy): 
      SparseMatrix(_nrows,_ncols), cols(ncols){
      for(int i=0; i<ncols; i++) cols[i]=new VECTOR(VECTOR::Zero(nrows));
    }
    
    MatrixX(const int _nrows, const int _ncols, const initializer_list<iivtriple>& list): 
      MatrixX<VECTOR>(_nrows,_ncols){
      for(const iivtriple& p:list){
	assert(p.first<nrows); assert(p.second<ncols); 
	set(p.first,p.second,p.third);
      }
    }
    
    ~MatrixX(){for(int i=0; i<cols.size(); i++) delete cols[i];}
    
    
  public: // copying ------------------------------------------------------------------------------------------------

    
    MatrixX(const MatrixX<VECTOR>& x): 
      SparseMatrix(x.nrows,x.ncols), 
      cols(x.ncols){
      COPY_WARNING("MatrixX<VECTOR>");
      for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j],_NoWarn());
    }
    
    MatrixX(const MatrixX<VECTOR>& x, const _NoWarn dummy): 
      SparseMatrix(x.nrows,x.ncols), 
      cols(x.ncols){
      for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j],_NoWarn());
    }

    MatrixX(MatrixX<VECTOR>&& x): 
      SparseMatrix(x.nrows,x.ncols), 
      cols(x.ncols){
      MOVE_WARNING("MatrixX<VECTOR>");
      cols=x.cols; 
      x.cols.clear();
      x.nrows=0; x.ncols=0;
    }

    MatrixX(MatrixX<VECTOR>&& x, const _NoWarn dummy): 
      SparseMatrix(x.nrows,x.ncols), 
      cols(x.ncols){
      cols=x.cols; 
      x.cols.clear();
      x.nrows=0; x.ncols=0;
    }
    
    MatrixX<VECTOR>& operator=(const MatrixX<VECTOR>& x){
      ASSIGN_WARNING("MatrixX<VECTOR>");
      nrows=x.nrows; ncols=x.ncols;
      for(auto p:cols) delete p; cols.resize(ncols); 
      for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j]);
      return *this;
    }
    
    MatrixX<VECTOR>& operator=(MatrixX<VECTOR>&& x){
      MOVEASSIGN_WARNING("MatrixX<VECTOR>");
      nrows=x.nrows; ncols=x.ncols;
      for(auto p:cols) delete p; 
      cols=x.cols; 
      x.nrows=0; x.ncols=0;
      return *this;
    }

    MatrixX<VECTOR> copy(){
      MatrixX M(nrows,ncols);
      for(int j=0; j<ncols; j++) M.cols[j]=new VECTOR(cols[j]->copy());
      return *this;
    }

    void detach(){cols.clear();}

    void assign(const MatrixX& x){
      assert(x.nrows==nrows);
      assert(x.ncols==ncols);
      for(int j=0; j<ncols; j++) *cols[j]=*x.cols[j]; 
    }

    MatrixX shallow() const{
      MatrixX M(nrows,ncols); 
      M.cols=cols;
      return M;
    }


  public: // conversions --------------------------------------------------------------------------------------------


    template<class VECTOR2>
    MatrixX(const MatrixX<VECTOR2>& x): 
      SparseMatrix(x.nrows,x.ncols), cols(x.ncols){
      CONVERT_WARNING("MatrixX<VECTOR2>","MatrixX<VECTOR>");
      for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j],_NoWarn());
    }

    MatrixX(const Cmatrix& x):
      MatrixX<VECTOR>(x.nrows,x.ncols){
      CONVERT_WARNING("Cmatrix","MatrixX<VECTOR>");
      for(int j=0; j<ncols; j++)
	for(int i=0; i<nrows; i++)
	  if(x.array[j*nrows+i]!=0) cols[j]->insert(i,x.array[j*nrows+i]);
    }

    MatrixX(const Cmatrix& x, const _NoWarn dummy):
      MatrixX<VECTOR>(x.nrows,x.ncols){
      for(int j=0; j<ncols; j++)
	for(int i=0; i<nrows; i++)
	  if(x.array[j*nrows+i]!=0) cols[j]->insert(i,x.array[j*nrows+i]);
    }

    template<class MATRIX>
    MatrixX(const MatrixTranspose<MATRIX>& x):
      MatrixX(x.ncols,x.nrows,_Zero()){
      TRANSPOSE_WARNING("MatrixX<VECTOR>");
      if(x.isSparseFormat()){
	x.obj.for_each_filled([this](const int i, const int j, const SCALAR v){cols[i]->append(j,v);}); return;}
      for(int j=0; j<ncols; j++)
	for(int i=0; i<nrows; i++){
	  SCALAR t=x.obj(j,i); if(t!=0) cols[j]->append(i,t);}
    }

    operator Cmatrix(){
      CONVERT_WARNING("MatrixX<VECTOR>","Cmatrix");
      Cmatrix M=Cmatrix::Zero(nrows,ncols);
      for_each_filled([&M](const int i, const int j, const SCALAR v){M(i,j)=v;});
      return M;
    }
  

  public: // named constructors ------------------------------------------------------------------------------------


    static MatrixX Zero(const int _nrows, const int _ncols){
      return MatrixX<VECTOR>(_nrows,_ncols,_Zero());}

    static MatrixX Filled(const int _nrows, const int _ncols, const SCALAR v){
      MatrixX<VECTOR> M(_nrows,_ncols);
      for(int j=0; j<_ncols; j++) (*M.cols[j])=VECTOR::Filled(_nrows,v); 
      return M;
    }

    static MatrixX Identity(const int n){
      MatrixX<VECTOR> M(n,n,_Zero());
      for(int i=0; i<n; i++) (*M.cols[i])(i)=1;
      return M;
    }

    static MatrixX Uniform(const int _nrows, const int _ncols){
      MatrixX<VECTOR> M(_nrows,_ncols);
      for(int j=0; j<_ncols; j++) (*M.cols[j])=VECTOR::Uniform(_nrows); 
      return M;
    }

    static MatrixX Gaussian(const int _nrows, const int _ncols){
      MatrixX<VECTOR> M(_nrows,_ncols);
      for(int j=0; j<_ncols; j++) (*M.cols[j])=VECTOR::Gaussian(_nrows); 
      return M;
    }

    static MatrixX Bernoulli(const int _nrows, const int _ncols, const double p=0.5){
      MatrixX<VECTOR> M(_nrows,_ncols);
      for(int j=0; j<_ncols; j++) (*M.cols[j])=VECTOR::Bernoulli(_nrows,p); 
      return M;
    }


  public: // comparisons --------------------------------------------------------------------------------------------


    template<class VECTOR2>
    bool operator==(const MatrixX<VECTOR2>& X) const{ 
      if(X.nrows!=nrows) return false; if(X.ncols!=ncols) return false;
      for(int j=0; j<ncols; j++) if(!(*cols[j]==*X.cols[j])) return false;
      return true;
    }


  public: // element access -----------------------------------------------------------------------------------------
  

    SCALAR operator()(const int i, const int j) const{
      assert(i<nrows); assert(j<ncols); return (*cols[j])(i);} 

    SCALAR read(const int i, const int j) const{
      assert(i<nrows); assert(j<ncols); return cols[j]->read(i);} 

    SCALAR& operator()(const int i, const int j){
      assert(i<nrows); assert(j<ncols); return (*cols[j])(i);}

    void set(const int i, const int j, const SCALAR v){ 
      assert(i<nrows); assert(j<ncols); (*cols[j])(i)=v;} 

    bool isFilled(const int i, const int j) const{
      return cols[j]->isFilled(i);}

    int nFilled() const{
      int t=0; for(int j=0; j<ncols; j++) t+=cols[j]->nFilled(); return t;}

    template<class VECTOR2>
    VECTOR2 column(const int j) const {assert(j<ncols);
      return VECTOR2(*cols[j]);}


  public: // function mappings --------------------------------------------------------------------------------------


    void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
      for(int j=0; j<ncols; j++) cols[j]->for_each_filled([j,&lambda](const int i, SCALAR& v){lambda(i,j,v);});}

    void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR)> lambda) const{
      for(int j=0; j<ncols; j++) cols[j]->for_each_filled([j,&lambda](const int i, const SCALAR v){lambda(i,j,v);});}

    void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      for(int j=0; j<ncols; j++) if(cols[j]->isFilled(i)) lambda(j,(*cols[j])(i));} // improve!

    void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR)> lambda) const{
      for(int j=0; j<ncols; j++) if(cols[j]->isFilled(i)) lambda(j,(*cols[j])(i));} // improve!

    void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      cols[j]->for_each_filled(lambda);}

    void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR)> lambda) const{
      static_cast<const VECTOR*>(cols[j])->for_each_filled(lambda);}

    template<class TYPE>
    TYPE accumulate_over_columns(std::function<TYPE(const VECTOR&)> lambda, 
      std::function<TYPE(const TYPE&, const TYPE&)> accumulator=std::plus<TYPE>(),
      const TYPE& t0=TYPE()) const{
      TYPE t(t0);
      for(int j=0; j<ncols; j++) 
	t=accumulator(t,lambda(*cols[j]));
      return t;
    }

    iipair find_best_over_columns(std::function<ivpair(const VECTOR&)> lambda, 
      std::function<bool(const SCALAR&, const SCALAR&)> selector) const{
      assert(ncols>0);
      ivpair p=lambda(*cols[0]);
      INDEX bestj=0;
      INDEX besti=p.first;
      SCALAR bestv=p.second;
      for(int j=1; j<ncols; j++){
	ivpair p=lambda(*cols[j]);	      
	if(selector(bestv,p.second)){bestj=j; besti=p.first; bestv=p.second;}
      }
      return iipair(besti,bestj);
    }


  public: // scalar valued operations -------------------------------------------------------------------------------


    SCALAR max() const{
      return accumulate_over_columns<SCALAR>([](const VECTOR& v){return v.max();},
	[](const SCALAR& a, const SCALAR& b){return std::max(a,fabs(b));});}
    SCALAR max_abs() const{
      return accumulate_over_columns<SCALAR>([](const VECTOR& v){return v.max_abs();},
	[](const SCALAR& a, const SCALAR& b){return std::max(a,fabs(b));});}
    iipair argmax() const{
      return find_best_over_columns([](const VECTOR& v){INDEX i=v.argmax(); return ivpair(i,v.read(i));}, 
	std::less<SCALAR>());}
    iipair argmax_abs() const{
      return find_best_over_columns([](const VECTOR& v){INDEX i=v.argmax_abs(); return ivpair(i,fabs(v.read(i)));}, 
	std::less<SCALAR>());}

   SCALAR min() const{
      return accumulate_over_columns<SCALAR>([](const VECTOR& v){return v.min();},
	[](const SCALAR& a, const SCALAR& b){return std::min(a,fabs(b));});}
    SCALAR min_abs() const{
      return accumulate_over_columns<SCALAR>([](const VECTOR& v){return v.min_abs();},
	[](const SCALAR& a, const SCALAR& b){return std::min(a,fabs(b));});}
    iipair argmin() const{
      return find_best_over_columns([](const VECTOR& v){INDEX i=v.argmin(); return ivpair(i,v.read(i));}, 
	std::greater<SCALAR>());}
    iipair argmin_abs() const{
      return find_best_over_columns([](const VECTOR& v){INDEX i=v.argmin_abs(); return ivpair(i,fabs(v.read(i)));}, 
	std::greater<SCALAR>());}

    SCALAR sum() const{
      SCALAR t=0; for(int j=0; j<ncols; j++) t+=cols[j]->sum(); return t;}
    SCALAR norm1() const{
      SCALAR t=0; for(int j=0; j<ncols; j++) t+=cols[j]->norm1(); return t;}
    SCALAR norm2() const{
      SCALAR t=0; for(int j=0; j<ncols; j++) t+=cols[j]->norm2(); return t;}
    SCALAR diff2(const MatrixX<VECTOR>& X) const {assert(X.ncols==ncols);
      SCALAR t=0; for(int j=0; j<ncols; j++) t+=cols[j]->diff2(*X.cols[j]); return t;}
    int nnz() const{
      int result=0; for(int j=0;j<ncols;j++) result+=cols[j]->nnz(); return result;}
    SCALAR spectralNorm() const;

    SCALAR inp_of_columns(const int j1, const int j2) const{
      return cols[j1]->dot(*cols[j2]);}

    SCALAR inp_of_rows(const int i1, const int i2) const{
      SCALAR t=0; for(int j=0; j<ncols; j++) t+=(cols[j]->read(i1))*(cols[j]->read(i2)); return t;}


  public: // sparse matrix methods ----------------------------------------------------------------------------------


    void sort() {for(auto p:cols) p->sort();}

    VectorSpattern nonemptyRows() const{
      return VectorSpattern(cols);
    }


  public: // in place operations ------------------------------------------------------------------------------------


    MatrixX<VECTOR>& operator+=(const MatrixX<VECTOR>& x){ assert(ncols==x.ncols);
      for(int j=0; j<ncols; j++) (*cols[j])+=(*x.cols[j]); return *this;}

    MatrixX<VECTOR>& operator-=(const MatrixX<VECTOR>& x){ assert(ncols==x.ncols);
      for(int j=0; j<ncols; j++) (*cols[j])-=(*x.cols[j]); return *this;}

    template<class VECTOR2>
    MatrixX<VECTOR>& multiplyRowsBy(const VECTOR2& v){
      for(int j=0; j<ncols; j++) (*cols[j])*=v; return *this;}

    template<class VECTOR2>
    MatrixX<VECTOR>& multiplyColumnsBy(const VECTOR2& v){
      assert(v.n==ncols); for(int j=0; j<ncols; j++) (*cols[j])*=v(j); return *this;}

    template<class VECTOR2>
    MatrixX<VECTOR>& divideRowsBy(const VECTOR2& v){
      for(int j=0; j<ncols; j++) (*cols[j])/=v; return *this;}

    template<class VECTOR2>
    MatrixX<VECTOR>& divideColumnsBy(const VECTOR2& v){
      assert(v.n==ncols); for(int j=0; j<ncols; j++) (*cols[j])*=(1.0/v(j)); return *this;}

    //void transpose(); 

    MatrixX<VECTOR>& tidy(){
      for(int j=0; j<ncols; j++) cols[j]->tidy(); return *this;}

    void symmetrize();


  public: // construct operations -----------------------------------------------------------------------------------

    
    MatrixX(const xcgram<Cmatrix>& e): 
      MatrixX(e.A.ncols,e.A.ncols,_Zero()){
      const Cmatrix& M=e.A;
      const int n=M.ncols;
      const int m=M.nrows;
      for(int i=0; i<n; i++)
	for(int j=0; j<=i; j++){
	  SCALAR t=0;
	  for(int k=0; k<m; k++)
	    t+=M.array[i*m+k]*M.array[j*m+k];
	  if(t==0) continue;
	  (*this)(i,j)=t;
	  (*this)(j,i)=t;
	}
    }

    template<class VECTOR2>
    MatrixX(const xcgram<MatrixX<VECTOR2> >& e): 
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()){
      const MatrixX<VECTOR2>& M=e.A;
      const int n=M.ncols;
      for(int i=0; i<n; i++)
	for(int j=0; j<=i; j++){
	  SCALAR t=M.cols[i]->dot(*M.cols[j]);
	  if(t==0) continue;
	  (*this)(i,j)=t;
	  (*this)(j,i)=t;
	}
    }
    
   template<class VECTOR2>
    MatrixX(const xcgram<SymmMatrixX<VECTOR2> >& e): 
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()){
      const MatrixX<VECTOR2>& M=e.A;
      const int n=M.ncols;
      for(int i=0; i<n; i++){
	VectorSpattern s;
	M.cols[i]->for_each_filled([&M,&s,i](const int u, const SCALAR x){s.add(*M.cols[u],i);});
	s.for_each_filled([this,&M,i](const int j){
			    SCALAR t=M.cols[i]->dot(*M.cols[j]);
			    if(t!=0){
			      (*this)(i,j)=t;
			      (*this)(j,i)=t;
			    }
			  });
      }
    }

    
   template<class VECTOR2>
   MatrixX(const xcgram<SymmMatrixXm<VECTOR2> >& e): 
     MatrixX<VECTOR>(gram(as_symmetric(MatrixX<VECTOR>(e.A)))){}

    
  public: // vector valued operations ------------------------------------------------------------------------------


    template<class VECTOR2>
    VECTOR2 operator*(const VECTOR2& x) const{
      VECTOR2 r = VECTOR2::Zero(nrows);
      assert(x.n==ncols);
      for(int j=0; j<ncols; j++){
	SCALAR t=x(j); //.array[j];
	cols[j]->for_each_filled([&r,t](const INDEX i, const SCALAR v){r(i)+=v*t;});
      }
      return r;
    }

    Cvector dot(const VECTOR& v) const{
      assert(v.n==nrows);
      Cvector r(ncols);
      for(int j=0; j<ncols; j++) r.array[j]=cols[j]->dot(v);
      return r;
    }


  public: // matrix valued operations --------------------------------------------------------------------------------


    //Cmatrix operator*(const MatrixX<VECTOR>& x) const{
    //  Cmatrix M(nrows,x.ncols);
    //  assert(ncols==nrows); 
    //  return M;
    //}

    Cmatrix dot(const MatrixX<VECTOR>& x) const{
      Cmatrix M(ncols,x.ncols);
      assert(x.nrows==nrows); 
      if(this==&x) for(int i=0; i<ncols; i++) for(int j=0; j<=i; j++) {M(i,j)=cols[i]->dot(*x.cols[j]); M(j,i)=M(i,j);}
      else for(int i=0; i<ncols; i++) for(int j=0; j<x.ncols; j++) M(i,j)=cols[i]->dot(*x.cols[j]);
      return M;
    }

    template<class MATRIX>
    MATRIX dot(const MatrixX<VECTOR>& x) const{}


  public: // Givens rotations ---------------------------------------------------------------------------------------


    MatrixX<VECTOR>& applyFromLeft(const GivensRotation& Q){
      for(int i=0; i<ncols; i++) cols[i]->apply(Q);
      return *this;
    }

    MatrixX<VECTOR>& applyFromLeftT(const GivensRotation& Q){
      for(int i=0; i<ncols; i++) cols[i]->applyT(Q);
      return *this;
    }

    MatrixX<VECTOR>& applyFromRightT(const GivensRotation& Q){
      VECTOR* newcol1=new VECTOR(nrows);
      VECTOR* newcol2=new VECTOR(nrows);
      assert(Q.i1<ncols); assert(Q.i2<ncols);
      newcol1->add(*cols[Q.i1],Q.cos); newcol1->add(*cols[Q.i2],-Q.sin);
      newcol2->add(*cols[Q.i2],Q.cos); newcol2->add(*cols[Q.i1],Q.sin);
      delete cols[Q.i1]; cols[Q.i1]=newcol1;
      delete cols[Q.i2]; cols[Q.i2]=newcol2;
      return *this;
    }

    MatrixX<VECTOR>& applyFromRight(const GivensRotation& Q){
      VECTOR* newcol1=new VECTOR(nrows);
      VECTOR* newcol2=new VECTOR(nrows);
      assert(Q.i1<ncols); assert(Q.i2<ncols);
      newcol1->add(*cols[Q.i1],Q.cos); newcol1->add(*cols[Q.i2],Q.sin);
      newcol2->add(*cols[Q.i2],Q.cos); newcol2->add(*cols[Q.i1],-Q.sin);
      delete cols[Q.i1]; cols[Q.i1]=newcol1;
      delete cols[Q.i2]; cols[Q.i2]=newcol2;
      return *this;
    }

    MatrixX<VECTOR>& conjugate(const GivensRotation& Q){
      applyFromLeft(Q);
      applyFromRightT(Q);
      return *this;
    }

    MatrixX<VECTOR>& conjugateT(const GivensRotation& Q){
      applyFromLeftT(Q);
      applyFromRight(Q);
      return *this;
    }

  
  public: // KpointOp<k> -------------------------------------------------------------------------------------------


    template<int k>
    MatrixX<VECTOR>& applyFromLeft(const KpointOp<k>& Q){
      for(int i=0; i<ncols; i++) cols[i]->apply(Q);
      return *this;
    }

    template<int k>
    MatrixX<VECTOR>& applyFromLeftT(const KpointOp<k>& Q){
      for(int i=0; i<ncols; i++) cols[i]->applyT(Q);
      return *this;
    }

    template<int k>
    MatrixX<VECTOR>& applyFromRightT(const KpointOp<k>& Q){
      VectorSpattern spattern;
      for(int i=0; i<k; i++){
	assert(Q.map(i)<ncols);
	spattern.add(*cols[Q.map(i)]);
      }
      for(auto r:spattern.filled){
	SCALAR* vptr[k];
	for(int i=0; i<k; i++)
	  vptr[i]=cols[Q.map(i)]->ptr(r);
	Q.applyTo(vptr);
      }      
      return *this;
    }

    template<int k>
    MatrixX<VECTOR>& applyFromRight(const KpointOp<k>& Q){
      VectorSpattern spattern;
      for(int i=0; i<k; i++){
	assert(Q.map(i)<ncols);
	spattern.add(*cols[Q.map(i)]);
      }
      for(auto r:spattern.filled){
	SCALAR* vptr[k];
	for(int i=0; i<k; i++)
	  vptr[i]=cols[Q.map(i)]->ptr(r);
	Q.applyToT(vptr);
      }      
      return *this;
    }

    template<int k>
    MatrixX<VECTOR>& conjugate(const KpointOp<k>& Q){
      applyFromLeft(Q);
      applyFromRightT(Q);
      return *this;
    }

    template<int k>
    MatrixX<VECTOR>& conjugateT(const KpointOp<k>& Q){
      applyFromLeftT(Q);
      applyFromRight(Q);
      return *this;
    }

  
  public: // I/O ----------------------------------------------------------------------------------------------------

    
    string str() const{
      ostringstream oss;
      oss.precision(3); 
      oss.setf(ios_base::fixed, ios_base::floatfield);
      for(int i=0; i<nrows; i++){oss<<"[ ";
	for(int j=0; j<ncols; j++) 
	  if(cols[j]->isFilled(i)) {oss.width(6); oss<<cols[j]->read(i)<<" ";}
	  else oss<<" 0.0   ";
	oss<<" ]\n";}
      return oss.str();
    }


    string str(const _Sparse dummy) const{
      ostringstream stream;
      stream.precision(5);
      for(int j=0; j<ncols; j++) 
	cols[j]->for_each_filled([&stream,j](const INDEX i, const SCALAR v){
				   stream<<"("<<i<<","<<j<<") : "<<v<<endl;
				 });
      return stream.str();
    }


  public: // Python interface --------------------------------------------------------------------------------------


    MatrixX(double* numpyDblArray, int numpyDim1, int numpyDim2): 
      MatrixX<VECTOR>(numpyDim1, numpyDim2){
      for(int j=0; j<ncols; j++)
	for(int i=0; i<nrows; i++)
	  if(numpyDblArray[i*ncols+j]!=0) cols[j]->insert(i,numpyDblArray[i*ncols+j]);
    }
    
    void np(double** numpyDblArray, int* numpyDim1, int* numpyDim2){
      *numpyDim1=nrows;
      *numpyDim2=ncols;
      *numpyDblArray=new double[nrows*ncols];
      std::fill(*numpyDblArray,*numpyDblArray+nrows*ncols,0);
      for_each_filled([this,numpyDblArray](const int i, const int j, const SCALAR v){
	  (*numpyDblArray)[i*ncols+j]=v;});
    }
    
    const char* __str__() { // why doesn't this work?
      ostringstream ostream; 
      cout<<str(_Sparse())<<endl;
      return ostream.str().c_str();
    }
   

  };






  // ---- Operations ----------------------------------------------------------------------------------------------


 
  /*
  template<class VECTOR>
  void MatrixX<VECTOR>::transpose(){ 
    vector<VECTOR*> newcols(nrows); 
    for(int i=0; i<nrows; i++) newcols[i]=new VECTOR(ncols);
    for(int i=0; i<ncols; i++){
      for(auto& p:*cols[i])
	if(p.second!=0) newcols[p.first]->append(i,p.second);
      delete cols[i];
    }
    cols=newcols; swapp(nrows,ncols);
  }
  */

  template<class VECTOR>
  void MatrixX<VECTOR>::symmetrize() {
    cout<<"Symmetrization unimplemented"<<endl;
    /*
      for(int j=0;j<ncols;j++){
      for(auto& p: *cols[j]) {
      if (p.first<=j){
      SCALAR& val=(*cols[p.first])(j);
      val += p.second;
      p.second = val;
      } else {
      SCALAR& val=(*cols[p.first])(j);
      if(val==0)
      val += p.second;
      }
      }
      }
    */
  }


  // ---- Scalar valued operations ---------------------------------------------------------------------------------



  template<class VECTOR>
  SCALAR MatrixX<VECTOR>::spectralNorm() const{  
    SCALAR norm=0;
    Cvector v = Cvector::Gaussian(ncols);
    SCALAR oldnorm=0;
  
    for(int i=0; i<200; i++){
      cout<<"test"<<endl;
      v=(*this)*v; //-M*v;
      norm=v.norm2();
      if(fabs((norm-oldnorm)/norm)<0.001) break;
      oldnorm=norm;
      v*=1.0/sqrt(norm);
    }
    return norm;
  }

  // ---- KpointOp<k> ----------------------------------------------------------------------------------------------

  /*
  template<>
  template<int k>
  MatrixX<Vectorv>& MatrixX<Vectorv>::applyFromRightT(const KpointOp<k>& Q){
    VectorSpattern spattern;
    for(int i=0; i<k; i++){
      assert(Q.map(i)<ncols);
      spattern.add(*cols[Q.map(i)]);
    }
    for(auto r:spattern.filled){
      SCALAR* vptr[k];
      for(int i=0; i<k; i++)
	vptr[i]=cols[Q.map(i)]->ptr(r);
      Q.applyTo(vptr);
    }      
    return *this;
  }
  */


  /*
  {
    VECTOR* col[Q.k]; 
    for(int j=0; j<Q.k; j++){
      assert(Q.map(j)<ncols);
      col[j]=new VECTOR(nrows);
      for(int i=0; i<Q.k; i++) col[j]->add(*cols[Q.map(i)],Q.q[i*Q.k+j]);}
    for(int j=0; j<Q.k; j++) {delete cols[Q.map(j)]; cols[Q.map(j)]=col[j];}
    return *this;
  }
  */


} // namespace Mondrian

#endif 



// =================================================================================================================


// ---- I/O --------------------------------------------------------------------------------------------------------


/*
template<class VECTOR>
MatrixX<VECTOR>::MatrixX(MatrixIF& file): MatrixX<VECTOR>(file.nrows,file.ncols){
  file.rewind();
  IndexValueTriple t;
  file>>t;
  while(t.i>=0){
    cols[t.j]->insert(t.i,t.value);
    file>>t;
  }
}


template<class VECTOR>
void MatrixX<VECTOR>::saveTo(MatrixOF& file) const{
  assert(file.nrows==nrows);
  assert(file.ncols==ncols);
  if(file.sparse){
    (foreach)([&file](const int i, const int j, const SCALAR v){file<<IndexValueTriple(i,j,v);});
  }else{
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	file<<read(i,j);
  }
}


template<class VECTOR>
string MatrixX<VECTOR>::classname(){return "MatrixX<"+VECTOR::classname()+">";}


template<class VECTOR>
MatrixX<VECTOR>::MatrixX(Bifstream& ifs):SparseMatrix(0,0){
  ifs.check(classname().c_str(),0);
  ifs.read(nrows);
  ifs.read(ncols);
  for(int j=0; j<ncols; j++)
    cols.push_back(new VECTOR(ifs));
};


template<class VECTOR>
void MatrixX<VECTOR>::serialize(Bofstream& ofs) const{
  ofs.tag(classname().c_str(),0);
  ofs.write(nrows);
  ofs.write(ncols);
  for(int j=0; j<ncols; j++) 
    cols[j]->serialize(ofs);
};


template<class VECTOR>
void MatrixX<VECTOR>::serialize(Rstream& rstream) const{
  rstream<<"MatrixX{"<<Rstream::endl;
  rstream.var("nrows",nrows);
  rstream.var("ncols",ncols);
  for(int j=0; j<ncols; j++){
    rstream<<"  cols["<<j<<"]=";
    cols[j]->serialize(rstream);
  }
  rstream<<"}"<<Rstream::endl;

};
*/

/*
template<class VECTOR>
string MatrixX<VECTOR>::str(const Sparse dummy) const{
  ostringstream stream;
  stream.precision(5);
  for(int j=0; j<ncols; j++) 
    for(auto& it: *cols[j]) stream<<"("<<it.first<<","<<j<<") : "<<it.second<<endl;
      return stream.str();
  }
*/

  /*
template<class VECTOR>
CSCmatrix MatrixX<VECTOR>::cscformat() {
 CSCmatrix result;
  result.nnz = nnz();
  result.ir = new INDEX[result.nnz];
  result.jc = new INDEX[ncols+1];
  result.val = new SCALAR[result.nnz];
  result.nrows = nrows; result.ncols = ncols; 

  result.jc[ncols] = result.nnz;

  int counter=0;
  for(int j=0;j<ncols;j++) {
    result.jc[j]=counter;
    for (auto& it: *cols[j]) {
      result.ir[counter] = it.first; 
      result.val[counter] = it.second; 
      counter++;
    }
  }
  return result;
  }

template<class VECTOR>
void MatrixX<VECTOR>::dump(SCALAR* result)  {
    int counter=0;
    for (int j=0;j<ncols;j++) {
      int i=0;
      for(auto&it: *cols[j]) {
        while(i!=it.first) {
          result[counter]=0; counter++;
          i++;
        }
        result[counter]=it.second; counter++;
      }
    }
  }
  */

// #ifdef HAVE_EIGEN
// template<class VECTOR>
// EigenSparseMatrix MatrixX<VECTOR>::convertToEigen() const{
//   EigenSparseMatrix M(nrows,ncols);
//   vector<Eigen::Triplet<double>> triplets;
//   for(int j=0;j<ncols;j++) {
//     for (auto& it: *cols[j]) {
//       triplets.push_back(Eigen::Triplet<double>(it.first,j,it.second));
//     }
//   }
// 	M.setFromTriplets(triplets.begin(), triplets.end());
// 	return M;
// }
// #endif


// #ifdef HAVE_EIGEN
// template<class VECTOR>
// MatrixX<VECTOR>::operator EigenSparseMatrix() const{
//   EigenSparseMatrix M(nrows,ncols);
//   vector<Eigen::Triplet<double>> triplets;
//   for(int j=0;j<ncols;j++) {
//     for (auto& it: *cols[j]) {
//       triplets.push_back(Eigen::Triplet<double>(it.first,j,it.second));
//     }
//   }
//   M.setFromTriplets(triplets.begin(), triplets.end());
//   return M;
// }
// #endif




//template<class VECTOR>
//ostream& operator<<(ostream& stream, const MatrixX<VECTOR>& x){stream<<x.str(Sparse()); return stream;}




  /*
template<class VECTOR>
template<class VECTOR2> // Warning: won't work for mapping unsorted v or h to l 
MatrixX<VECTOR>::MatrixX(const MatrixX<VECTOR2>& x, 
					 const Remap& remap1, const Remap& remap2, 
					 const int _nrows, const int _ncols):
  MatrixX(_nrows,_ncols){
  for(int j=0; j<x.ncols; j++){
    const VECTOR& xcol=*x.cols[j];
    VECTOR& col=*cols[remap2.forward[j]];
    for(auto& p:xcol)
      col.append(remap1.forward[p.first],p.second);
  }
}
  */
  //void applyFromLeft(const KpointRotation& Q){
  //  for(int i=0; i<ncols; i++) cols[i]->apply(Q);}

  //void applyFromLeftT(const KpointRotation& Q){
  //  for(int i=0; i<ncols; i++) cols[i]->applyInverse(Q);}

  /*
  void applyFromRightT(const KpointRotation& Q){
    VECTOR* col[Q.k]; 
    for(int j=0; j<Q.k; j++){
      assert(Q.ix[j]<ncols);
      col[j]=new VECTOR(nrows);
      for(int i=0; i<Q.k; i++) col[j]->add(*cols[Q.ix[i]],Q.q[i*Q.k+j]);}
    for(int j=0; j<Q.k; j++) {delete cols[Q.ix[j]]; cols[Q.ix[j]]=col[j];}
  }
  */

  /*
  void applyFromRight(const KpointRotation& Q){
    VECTOR* col[Q.k]; 
    for(int j=0; j<Q.k; j++){
      assert(Q.ix[j]<ncols);
      col[j]=new VECTOR(nrows);
      for(int i=0; i<Q.k; i++) col[j]->add(*cols[Q.ix[i]],Q.q[j*Q.k+i]);}
    for(int j=0; j<Q.k; j++) {delete cols[Q.ix[j]]; cols[Q.ix[j]]=col[j];}
  }
  */

  /*
  void symmetrize(){ // improve!
    for(int j=0; j<ncols; j++)
      for(int i=0; i<j; j++){
	SCALAR t=(cols(j)->read(i)+cols(i)->read(j))/2.0;
	(*this)(i,j)=t;
	(*this)(j,i)=t;
      }
  }
  */


  //MatrixX(DenseMatrixFile& file);
  //MatrixX(SparseMatrixFile& file);

  //MatrixX(MatrixIF& file);
  //void saveTo(MatrixOF& file) const;
  //void dump(SCALAR* result);
  //CSCmatrix cscformat();

  //static string classname();
  //MatrixX<VECTOR>(Bifstream& ifs);
  //void serialize(Bofstream& ofs) const;
  //void serialize(Rstream& rstream) const;

  // ---- Conversions -----------------------------------------------------------------------------------------------


/*
template<class VECTOR> 
template<class VECTOR2>
MatrixX<VECTOR>::MatrixX(const MatrixX<VECTOR2>& x): 
  SparseMatrix(x.nrows,x.ncols), cols(x.ncols){
  for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j]);
  cout<<"Warning: MatrixX convert-copied."<<endl;
}

template<class VECTOR> 
template<class VECTOR2>
MatrixX<VECTOR>::MatrixX(const MatrixX<VECTOR2>& x, const int _nrows, const int _ncols): 
  SparseMatrix(_nrows,_ncols), cols(_ncols){
  for(int j=0; j<ncols; j++){
    cols[j]=new VECTOR(nrows);
    if (j<x.ncols) {
    for(auto& p:*x.cols[j])
      if(p.first<nrows) cols[j]->insert(p.first,p.second);
  }
  }
}

template<class VECTOR> 
template<class VECTOR2>
MatrixX<VECTOR>::MatrixX(const MatrixX<VECTOR2>& x, const Transpose& dummy):
  SparseMatrix(x.ncols,x.nrows), cols(x.nrows){
  for(int j=0; j<ncols; j++) cols[j]=new VECTOR(ncols);
  for(int i=0; i<nrows; i++){
    x.cols[i]->sort();
    for(auto& p:*x.cols[i])
      if(p.second!=0) cols[p.first]->append(i,p.second);
  }
}

template<class VECTOR>
template<class VECTOR2> // Warning: won't work for mapping unsorted v or h to l 
MatrixX<VECTOR>::MatrixX(const MatrixX<VECTOR2>& x, 
					 const Remap& remap1, const Remap& remap2, 
					 const int _nrows, const int _ncols):
  MatrixX(_nrows,_ncols){
  for(int j=0; j<x.ncols; j++){
    const VECTOR& xcol=*x.cols[j];
    VECTOR& col=*cols[remap2.forward[j]];
    for(auto& p:xcol)
      col.append(remap1.forward[p.first],p.second);
  }
}

template<class VECTOR>
MatrixX<VECTOR>::MatrixX(const class Cmatrix& x): 
  SparseMatrix(x.nrows,x.ncols), cols(x.ncols){
  for(int j=0; j<ncols; j++)
    cols[j]=new VECTOR(x.vcols(j));
  cout<<"WARNING: dense->sparse conversion"<<endl;
}
*/

/*
  private: // remapping ---------------------------------------------------------------------------------------------


    //MatrixX<VECTOR> remapRows(const IndexMap& map) const{
    //  MatrixX<VECTOR> M(map.ndest,ncols);
    //  for(int j=0; j<ncols; j++) *M.cols[j]=cols[j]->remap(map);
    //  return M;
    //}

    //MatrixX<VECTOR> remapRows(const Inverse<IndexMap>& map) const{
    //  MatrixX<VECTOR> M(map.obj.nsource,ncols);
    //  for(int j=0; j<ncols; j++) *M.cols[j]=cols[j]->remap(map.obj);
    //  return M;
    //}

    //MatrixX<VECTOR> remapCols(const IndexMap& map) const{
    //  assert(map.nsource==ncols);
    //  MatrixX<VECTOR> M(nrows,map.ndest); 
    //  for(int j=0; j<ncols; j++) M.cols[map(j)]=new VECTOR(*cols[j]);
    //  return M;
    //}

    //MatrixX<VECTOR> remapCols(const Inverse<IndexMap>& map) const{
    //  assert(map.obj.ndest==ncols);
    //  MatrixX<VECTOR> M(nrows,map.obj.nsource); 
    //  for(int j=0; j<ncols; j++) M.cols[j]=new VECTOR(*cols[map.obj(j)]);
    //  return M;
    //}

    //MatrixX<VECTOR>& applyFromLeft(const IndexMap& map) const{
    //  for(auto p:cols) p->apply(map);
    //  return *this;
    //}

    //MatrixX<VECTOR>& applyFromLeft(const Inverse<IndexMap>& map) const{
    //  for(auto p:cols) p->apply(map);
    //  return *this;
    //}

    //MatrixX<VECTOR>& applyFromRight(const IndexMap& map) const{
    //  assert(map.nsource==ncols); assert(map.ndest==ncols);
    //  map.applyTo(cols);
    //  symmetricp=false;
    //  return *this;
    //}

    //MatrixX<VECTOR>& applyFromRight(const Inverse<IndexMap>& map) const{
    //  assert(map.obj.nsource==ncols); assert(map.obj.ndest==ncols);
    //  map.obj.applyToInv(cols);
    //  symmetricp=false;
    //  return *this;
    //}
    */

    /*
  public: // KpointOperator

    MatrixX<VECTOR>& applyFromLeft(const KpointOperator& Q){
      for(int i=0; i<ncols; i++) cols[i]->apply(Q);
      return *this;
    }

    MatrixX<VECTOR>& applyFromLeftT(const KpointOperator& Q){
      for(int i=0; i<ncols; i++) cols[i]->applyT(Q);
      return *this;
    }

    MatrixX<VECTOR>& applyFromRightT(const KpointOperator& Q){
      VECTOR* col[Q.k]; 
      for(int j=0; j<Q.k; j++){
	assert(Q.map(j)<ncols);
	col[j]=new VECTOR(nrows);
	for(int i=0; i<Q.k; i++) col[j]->add(*cols[Q.map(i)],Q.q[i*Q.k+j]);}
      for(int j=0; j<Q.k; j++) {delete cols[Q.map(j)]; cols[Q.map(j)]=col[j];}
      return *this;
    }

    MatrixX<VECTOR>& applyFromRight(const KpointOperator& Q){
      VECTOR* col[Q.k]; 
      for(int j=0; j<Q.k; j++){
	assert(Q.map(j)<ncols);
	col[j]=new VECTOR(nrows);
	for(int i=0; i<Q.k; i++) col[j]->add(*cols[Q.map(i)],Q.q[j*Q.k+i]);}
      for(int j=0; j<Q.k; j++) {delete cols[Q.map(j)]; cols[Q.map(j)]=col[j];}
      return *this;
    }

    MatrixX<VECTOR>& conjugate(const KpointOperator& Q){
      applyFromLeft(Q);
      applyFromRightT(Q);
      return *this;
    }

    MatrixX<VECTOR>& conjugateT(const KpointOperator& Q){
      applyFromLeftT(Q);
      applyFromRight(Q);
      return *this;
    }
    */

     /*
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()){
      const MatrixX<VECTOR2>& M=e.A;
      const int n=M.ncols;
      for(int i=0; i<n; i++){
	VectorSpattern s;
	M.cols[i]->for_each_filled([&M,&s,i](const int u, const SCALAR x){s.add(*M.cols[u],i);});
	s.for_each_filled([this,&M,i](const int j){
			    SCALAR t=M.cols[i]->dot(*M.cols[j]);
			    if(t!=0){
			      (*this)(i,j)=t;
			      (*this)(j,i)=t;
			    }
			  });
      }
    }
    */

    /*
    template<class MATRIX>
    MatrixX(const MATRIX& x): 
      MatrixX<VECTOR>(x.nrows,x.ncols,_Zero()){
      CONVERT_WARNING("MATRIX","MatrixX<VECTOR>");
      if(x.isSparseFormat()) x.for_each_filled([this](const int i, const int j, const SCALAR v){cols[i]->append(j,v);});
      else
	for(int j=0; j<ncols; j++)
	  for(int i=0; i<nrows; i++){
	    SCALAR t=x(i,j); if(t!=0) cols[j]->append(i,t);}
    }
    */

