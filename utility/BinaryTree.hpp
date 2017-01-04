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


#ifndef _BinaryTree
#define _BinaryTree

#include "BinaryTreeNode.hpp"

namespace Mondrian{

  template<class DECORATION>
  class BinaryTree{
  public:

    typedef BinaryTreeNode<DECORATION> Node;

    Node* root=nullptr;


  public:

    ~BinaryTree(){delete root;}

  public:

    Node* addRoot(const DECORATION& decor){
      root=new Node(decor);
      return root;
    }

    Node* addNodeBalanced(const DECORATION& decor){
      if(root==nullptr) return addRoot(decor);
      Node* node=root;
      while(node->isFull())
	node=node->shallowestChild();
      return new Node(decor,node);
    }

    Node* addNodeBalancedPreservingLeaves(const DECORATION& decor){
      if(root==nullptr) return addRoot(decor);
      Node* node=root;
      while(node->isFull())
	node=node->shallowestChild();
 
      if(node->nChildren()==1){
	cout<<"Error"<<endl;
	return nullptr;
      }

     Node* grandparent=node->parent;
     if(grandparent!=nullptr) grandparent->unlink(node);
     Node* replacement=new Node(DECORATION(),grandparent);
     if(grandparent==nullptr) root=replacement;
     replacement->addChild(node);

     return new Node(decor,replacement);
    }

  public:

    string str(const int indent=0) const{
      if(root==nullptr) return "";
      return root->str(0);
    }


};


}
#endif
