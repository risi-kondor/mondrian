#ifndef _Interruptable
#define _Interruptable

class Interruptable{
public:

  void halt(){interrupt=true;}

public:

  bool interrupt=false;

};


#endif
