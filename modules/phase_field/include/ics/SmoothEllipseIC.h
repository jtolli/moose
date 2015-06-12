/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef SMOOTHELLIPSEIC_H
#define SMOOTHELLIPSEIC_H

#include "Kernel.h"
#include "InitialCondition.h"

// System includes
#include <string>

// Forward Declarations
class SmoothEllipseIC;

template<>
InputParameters validParams<SmoothEllipseIC>();

/**
 * SmoothEllipseIC creates an ellipse
 */
class SmoothEllipseIC : public InitialCondition
{
public:
  /**
   * Constructor
   *
   * @param name The name given to the initial condition in the input file.
   * @param parameters The parameters object holding data for the class to use.
   * @param var_name The variable this InitialCondtion is supposed to provide values for.
   */
  SmoothEllipseIC(const std::string & name,
                 InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */

protected:
  MooseMesh & _mesh;

  Real _invalue;
  Real _outvalue;
  Real _int_width;
  bool _3D_spheres;
  unsigned int _num_dim;
  Real _x1;
  Real _y1;
  Real _z1;
  Real _x2;
  Real _y2;
  Real _z2;
  Real _distance;
  Point _focalPoint1;
  Point _focalPoint2;
    
  

  

 
  virtual Real value(const Point & p);

  virtual RealGradient gradient(const Point & p);

  virtual Real computeEllipseValue(const Point & p, const Point & fp1, const Point & fp2, const Real & distance);

  virtual Point computeEllipseGradient(const Point & p, const Point & fp1, const Point &fp2, const Real & radius);

};

#endif //SMOOTHELLIPSEIC_H
