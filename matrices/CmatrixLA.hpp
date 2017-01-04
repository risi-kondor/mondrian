/* ---------------------------------------------------------------------------
 
  Mondrian: open source multithreaded matrix library 

  Copyright (C) 2017 Imre Risi Kondor
  Copyright (C) 2015 Imre Risi Kondor, Nedelina Teneva, Pramod K Mudrakarta

  Parts of the following code are derived from the pMMF library 
  (https://github.com/risi-kondor/pMMF) which is licensed under the 
  GNU Public License, version 3. This code therefore is also licensed 
  under the terms of the GNU Public License, version 3. 
 
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 3
  of the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

 --------------------------------------------------------------------------- */


#ifndef _CmatrixLA
#define _CmatrixLA

#include "Cmatrix.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

  class AsCmatrixLA;  


class CmatrixLA: public Cmatrix{
public:

  using Cmatrix::Cmatrix;

public: // downcasting 

  CmatrixLA(const Cmatrix& x): 
    Cmatrix(x,_NoWarn()){
    //DOWNCASTCOPY_WARNING("Cmatrix","CmatrixLA");
  }

  CmatrixLA(Cmatrix&& x): 
    Cmatrix(std::move(x),_NoWarn()){
    //DOWNCAST_WARNING("Cmatrix","CmatrixLA");
  }

  CmatrixLA(const Cmatrix& x, const _Downcast dummy): 
    Cmatrix(x,_NoWarn()){
    //DOWNCASTCOPY_WARNING("Cmatrix","CmatrixLA");
  }

  CmatrixLA(Cmatrix&& x, const _Downcast dummy): 
    Cmatrix(std::move(x),_NoWarn()){
    //DOWNCAST_WARNING("Cmatrix","CmatrixLA");
  }


public: // conversions

  CmatrixLA(const AsCmatrixLA& x);

  CmatrixLA(const EigenMatrixXdAdaptor& X);
  operator EigenMatrixXdAdaptor() const;


public: // Linear algebra operations

  /*
  CmatrixLA& GSorthogonalize(){
    for(int j=0; j<ncols; j++){
      Detached<Cvector> vj=view_of_column(j);
      for(int i=0; i<j; i++){
	Detached<Cvector> vi=view_of_column(i);
	vj-=vi*(vj.dot(vi));
      }
      vj/=sqrt(vj.norm2());
    }
    return *this;
  }
  */

  package<CmatrixLA,Cvector> symmetricEVD() const;

};


  // wrapper allowing a Cmatrix to pose as a CmatrixLA  
  class AsCmatrixLA: public CmatrixLA{
  public:
    AsCmatrixLA()=delete;
    AsCmatrixLA& operator=(const AsCmatrixLA& x)=delete;
    AsCmatrixLA(Cmatrix& x): CmatrixLA(x.shallow(),_Downcast()){}
    AsCmatrixLA(const Cmatrix& x): CmatrixLA(x.shallow(),_Downcast()){}
    ~AsCmatrixLA(){array=nullptr;}
  };

  inline CmatrixLA::CmatrixLA(const AsCmatrixLA& x): 
    Cmatrix(x,_NoWarn()) {
    COPY_WARNING("CmatrixLA");
  }

  inline AsCmatrixLA linalg(Cmatrix& x){
    return AsCmatrixLA(x);}
  inline const AsCmatrixLA linalg(const Cmatrix& x){
    return AsCmatrixLA(x);}
  inline CmatrixLA linalg(Cmatrix&& x){
    MOVE_WARNING("Cmatrix");
    return CmatrixLA(std::move(x),_Downcast());}




} // namespace Mondrian

#endif
