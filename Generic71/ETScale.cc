// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ETScale class.
//

#include "ETScale.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

ETScale::ETScale() {}

ETScale::~ETScale() {}

IBPtr ETScale::clone() const {
  return new_ptr(*this);
}

IBPtr ETScale::fullclone() const {
  return new_ptr(*this);
}

Energy2 ETScale::renormalizationScale() const { 

  // mePartonData() returns a vector which contains the incoming, followed by
  // the outgoing particles presented to this scale object, meMomenta()
  // returns a vector with the momenta in the same order

  size_t k = 2;
  int top = -1;
  int antitop = -1;

  while ( (top == -1 || antitop == -1) && k < mePartonData().size() ){
    if ( mePartonData()[k]->id() == 6 ) {
      if ( top < 0 )
        top = k;
      else
        assert(false);
    } else if ( mePartonData()[k]->id() == -6 ) {
      if ( antitop < 0 )
        antitop = k;
      else
        assert(false);
    }
    k++;
  }

  if ( top < 2 || antitop < 2 ){
    throw Exception() << "MatchboxTopMTScale: Could not find a top-antitop-pair in the final state!\n"
                      << Exception::runerror;
  }
  // cerr << " sqrt(TopMTScale)  = "
  //   //      << sqrt(meMomenta()[top].mt2()+meMomenta()[antitop].mt2())/GeV
  //     //      << "\n" << flush;
  return (meMomenta()[top].mt()*meMomenta()[antitop].mt());
}

Energy2 ETScale::factorizationScale() const { 
  return renormalizationScale();
}

Energy2 ETScale::showerScale() const {
  return factorizationScale();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeNoPIOClass<ETScale,Herwig::MatchboxScaleChoice>
  describeHerwigETScale("Herwig::ETScale", "ETScale.so");

void ETScale::Init() {

  static ClassDocumentation<ETScale> documentation
    ("There is no documentation for the ETScale class");

}

