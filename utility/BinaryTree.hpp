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
