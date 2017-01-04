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


#ifndef _AccumulationTree
#define _AccumulationTree

#include <unordered_map>

#include "BinaryTree.hpp"


namespace Mondrian{

  template<class DECORATION, class ACCUMULATOR=std::plus<DECORATION> >
  class AccumulationTree: public BinaryTree<DECORATION>{
  public:

    typedef BinaryTreeNode<DECORATION> Node;

    unordered_map<INDEX,Node*> lookup;

    // const std::function<DECORATION(const DECORATION&, const DECORATION&)> accumulator;

  public:


    //AccumulationTree(const  std::function<DECORATION(const DECORATION&, const DECORATION&)>& _accumulator):
    //  accumulator(_accumulator){};


  public:

    Node* addLeaf(const INDEX ix, const DECORATION& decor){
      Node* node=BinaryTree<DECORATION>::addNodeBalancedPreservingLeaves(decor);
      recomputeFrom(node);
      lookup[ix]=node;
      return node;
    }


    void set(const INDEX ix, const DECORATION& decor){
      auto found=lookup.find(ix);
      if(found==lookup.end())
	addLeaf(ix,decor);
      found->second->decor=decor;
      recomputeFrom(found->second);
    }


    void recomputeFrom(const Node* leaf){
      Node* node=leaf->parent;
      while(node!=nullptr){
	if(node->nChildren()==1){
	  if(node->child1!=nullptr) node->decor=node->child1->decor;
	  else node->decor=node->child2->decor;
	}else{
	  ACCUMULATOR accumulator;
	  node->decor=accumulator(node->child1->decor,node->child2->decor);
	}
	node=node->parent;
      }
    }

  };

}

#endif
