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


#ifndef _BinaryTreeNode
#define _BinaryTreeNode

#include "Mondrian_base.hpp"

namespace Mondrian{

  template<class DECORATION>
  class BinaryTreeNode{
  public:

    typedef BinaryTreeNode<DECORATION> Node;
    
    BinaryTreeNode<DECORATION>* parent=nullptr; 
    BinaryTreeNode<DECORATION>* child1=nullptr; 
    BinaryTreeNode<DECORATION>* child2=nullptr;

    int insertion_depth=0;

    DECORATION decor;

  public:

    BinaryTreeNode(const DECORATION& _decor, Node* _parent=nullptr): decor(_decor){
      if(_parent!=nullptr) _parent->addChild(this);
    }

    ~BinaryTreeNode(){
      if(child1!=nullptr) delete child1;
      if(child2!=nullptr) delete child2;
      if(parent!=nullptr) parent->unlink(this);
    }
    
  public:
      
    bool isFull() const{
      return (child1!=nullptr&&child2!=nullptr);}
    
    int nChildren() const{
      return static_cast<int>(child1!=nullptr)+(child2!=nullptr);}

    Node* shallowestChild() const{
      if(child2->insertion_depth<child1->insertion_depth) return child2; 
      else return child1;
    }
    
    void addChild(Node* newchild){
      if(child1==nullptr) child1=const_cast<Node*>(newchild); else{
	if(child2==nullptr) child2=const_cast<Node*>(newchild); else{
	  cout<<"Error: node is full"<<endl; return;}}
      newchild->parent=this;
      computeInsertionDepth();	
    }
    
    void unlink(const Node* kid){
      if(child1==kid) child1=nullptr;
      if(child2==kid) child2=nullptr;
    }
    
    void computeInsertionDepth(){
      int t=0;
      if(isFull()) t=std::min(child1->insertion_depth,child2->insertion_depth)+1;
      if(insertion_depth!=t){
	insertion_depth=t;
	if(parent!=nullptr) parent->computeInsertionDepth();
      }
    }

  public:

    string str(const int indent=0) const{
      ostringstream oss;
      for(int i=0; i<indent; i++) oss<<"  ";
      oss<<"Node: "<<decor.str()<<endl;
      if(child1!=nullptr) oss<<child1->str(indent+1);
      if(child2!=nullptr) oss<<child2->str(indent+1);
      return oss.str();
    }

  };


}


#endif
