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

  pair<Lorentz5Momentum,Lorentz5Momentum> topMomenta;
  pair<bool,bool> topFound(false,false);

  auto dataIt = mePartonData().begin();
  auto momentumIt = meMomenta().begin();

  while ( !topFound.first && !topFound.second &&
	  dataIt != mePartonData().end() ) {
    if ( ((**dataIt).id() == 6 || (**dataIt).id() == -6) && !topFound.first ) {
      topMomenta.first = *momentumIt;
      topFound.first = true;
    }
    if ( ((**dataIt).id() == 6 || (**dataIt).id() == -6) &&
	 topFound.first && !topFound.second ) {
      topMomenta.second = *momentumIt;
      topFound.second = true;
    }
    ++dataIt; ++momentumIt;
  }

  if ( !topFound.first || !topFound.second )
    throw Exception()
      << "ETScale is expecting to find a top pair -- please check your process setup."
      << Exception::runerror;

  double result = topMomenta.first.mt()*topMomenta.second.mt();
  if (result<=0) throw Exception() << "Negative transverse mass!" << Exception::runerror;
  return result;

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

