/* -----------------------------------------------------------------------------
 
  Mondrian: open source multithreaded matrix library 

  Copyright (C) 2017 Imre Risi Kondor

 
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
 
----------------------------------------------------------------------------- */


#ifndef _BindexMap
#define _BindexMap

#include "Mondrian_base.hpp"


namespace Mondrian{

class BindexMap: public Serializable{
public:

  BindexMap(const int _k): k(_k), forward(_k) {}

  BindexMap(const initializer_list<iipair> list): BindexMap(list.size()){
    int i=0; for(auto& j:list) forward[i++]=j;}

  int k;
  vector<iipair> forward;
  
  
public: // copying 
  
  BindexMap(const BindexMap& x): k(x.k), forward(x.forward){}
  
  BindexMap(BindexMap&& x): k(x.k), forward(std::move(x.forward)) {}
					    
  BindexMap& operator=(const BindexMap& x){
    k=x.k; forward=x.forward; 
    return *this;}

  BindexMap& operator=(BindexMap&& x) {
    k=snatch(x.k); forward=std::move(x.forward); 
    return *this;
  }

  ~BindexMap(){}


public:

  iipair& operator()(const int i) {assert(i<k); return forward[i];}
  iipair operator()(const int i) const {assert(i<k); return forward[i];}


public:

  string str() const{
    ostringstream result; 
    //result<<"("; for(int i=0; i<k-1; i++) result<<forward[i]<<","; result<<forward[k-1]<<")"; 
    return result.str();
  }


public:


};


}

#endif

