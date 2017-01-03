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
