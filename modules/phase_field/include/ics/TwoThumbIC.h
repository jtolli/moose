/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef TWOTHUMBIC_H
#define TWOTHUMBIC_H

#include "Kernel.h"
#include "InitialCondition.h"

// System includes
#include <string>

// Forward Declarations
class TwoThumbIC;

template<>
InputParameters validParams<TwoThumbIC>();

/**
 * ThumbIC creates a rectangle with a half circle on top
 */
class TwoThumbIC : public InitialCondition
{
public:
  /**
   * Constructor
   *
   * @param name The name given to the initial condition in the input file.
   * @param parameters The parameters object holding data for the class to use.
   * @param var_name The variable this InitialCondtion is supposed to provide values for.
   */
  TwoThumbIC(const std::string & name,
                 InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);

protected:
  Real _xcoord;
  Real _ycoord;
  Real _width;
  Real _height;
  Real _invalue;
  Real _outvalue;
};

#endif //TWOTHUMBIC_H
