/* ---------------------------------------------------------------------------
 
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
 
--------------------------------------------------------------------------- */


#ifndef _OrderedSet
#define _OrderedSet

#include <set>
#include <unordered_map>

#include "Mondrian_base.hpp"

// todo: copy and move!!

namespace Mondrian{

  class OrderedSetLess{
  public:
    bool operator()(const ivpair& a, const ivpair& b){return a.second<b.second;} 
  };


  template<class COMPARATOR=OrderedSetLess>
  class OrderedSet{
  public:
    

    typedef typename multiset<ivpair,COMPARATOR>::iterator ordered_iterator;

    unordered_map< INDEX,  ordered_iterator> lookup;
    multiset<ivpair,COMPARATOR> ordered;


  public: 


    void set(const INDEX& i, const SCALAR& v){
      auto found=lookup.find(i);
      if(found!=lookup.end()){
	ordered.erase(found->second);
	lookup.erase(found);
      }
      auto it=ordered.insert(ivpair(i,v));
      lookup[i]=it;
    }


    INDEX best(int c=0){
      c=std::min(c,static_cast<int>(ordered.size())-1);
      if(c<0) return 0;
      auto it=ordered.begin();
      for(int i=0; i<c ; i++) it++;
      return it->first;
    }

    INDEX worst(int c=0){
      c=std::min(c,static_cast<int>(ordered.size())-1);
      if(c<0) return 0;
      auto it=ordered.rbegin();
      for(int i=0; i<c ; i++) it++;
      return it->first;
    }

    
  public: // I/O ----------------------------------------------------------------------------------------------------

    string str() const{
      ostringstream oss;
      for(auto& p:ordered)
	oss<<p.str()<<endl;
      return oss.str();
    }

  };

  template<class COMPARATOR>
  inline ostream& operator<<(ostream& stream, const OrderedSet<COMPARATOR>& x){stream<<x.str(); return stream;}


}

#endif
