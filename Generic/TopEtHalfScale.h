// -*- C++ -*-
#ifndef Herwig_TopEtHalfScale_H
#define Herwig_TopEtHalfScale_H
//
// This is the declaration of the TopEtHalfScale class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the TopEtHalfScale class.
 *
 * @see \ref TopEtHalfScaleInterfaces "The interfaces"
 * defined for TopEtHalfScale.
 */
class TopEtHalfScale: public Herwig::MatchboxScaleChoice {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TopEtHalfScale();

  /**
   * The destructor.
   */
  virtual ~TopEtHalfScale();
  //@}

public:

  /**
   * Return the renormalization scale. This default version returns
   * shat.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Return the factorization scale. This default version returns
   * shat.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Return the shower hard scale. This default implementation returns the
   * factorization scale.
   */
  virtual Energy2 showerScale() const;

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TopEtHalfScale & operator=(const TopEtHalfScale &);

};

}

#endif /* Herwig_TopEtHalfScale_H */
