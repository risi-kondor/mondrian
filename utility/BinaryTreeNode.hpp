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
