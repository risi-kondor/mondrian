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
