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


#ifndef _GivensRotation
#define _GivensRotation

#include "Mondrian_base.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{


  class GivensRotation{
  public:
    
    int i1;
    int i2;
    double cos;
    double sin;
    

  public:
    
    GivensRotation(){}
    
    GivensRotation(const GivensRotation& x){i1=x.i1; i2=x.i2; cos=x.cos; sin=x.sin;}
    
    GivensRotation operator=(GivensRotation&& x){
      i1=x.i1; i2=x.i2; cos=x.cos; sin=x.sin; 
      return *this;
    }
    
    GivensRotation(int _i1, int _i2, double _cos, double _sin){
      if(_i1<_i2){i1=_i1;i2=_i2;cos=_cos;sin=_sin;}
      else{i1=_i2;i2=_i1;cos=_cos;sin=-_sin;}}
    
    GivensRotation(int _i1, int _i2, double _alpha){
      if(_i1<_i2){i1=_i1; i2=_i2; cos=::cos(_alpha); sin=::sin(_alpha);}
      else{i2=_i1; i1=_i2; cos=::cos(_alpha); sin=-::sin(_alpha);}
    }

  public:

    static GivensRotation Random(const int n){
      uniform_int_distribution<int> distri(0,n-1);
      uniform_real_distribution<SCALAR> distrr(0,2*M_PI);
      int i1=distri(randomNumberGenerator); 
      int i2=i1; while(i2==i1){i2=distri(randomNumberGenerator);}
      if(i2>i1) swapp(i1,i2);
      SCALAR cos=std::cos(distrr(randomNumberGenerator)); 
      SCALAR sin=sqrt(1.0-cos*cos);
      return GivensRotation(i1,i2,cos,sin);
    }

    
  public:
    
    string str() const{
      ostringstream stream;
      stream.precision(3); 
      stream.setf(ios_base::fixed, ios_base::floatfield);
      stream<<"Rotation on ("<<i1<<","<<i2<<"):"<<endl;
      stream<<"[ "; stream.width(6); stream<<cos<<" "; stream.width(6); stream<<-sin<<"  ]"<<endl;
      stream<<"[ "; stream.width(6); stream<<sin<<" "; stream.width(6); stream<<cos<<"  ]"<<endl;
      return stream.str();
    }

    //string str() const;


  };


  inline ostream& operator<<(ostream& stream, const GivensRotation& x){stream<<x.str(); return stream;}

}

#endif
