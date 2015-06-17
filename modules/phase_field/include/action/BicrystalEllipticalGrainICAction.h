/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef BICRYSTALELLIPTICALGRAINICACTION_H
#define BICRYSTALELLIPTICALGRAINICACTION_H

#include "InputParameters.h"
#include "Action.h"

/**
 * Automatically generates all variables to model a polycrystal with op_num orderparameters
 */
class BicrystalEllipticalGrainICAction: public Action
{
public:
  BicrystalEllipticalGrainICAction(const std::string & name, InputParameters params);

  virtual void act();

private:
  static const Real _abs_zero_tol;

  std::string _var_name_base;
  unsigned int _op_num;

  Real _distance;
  Real _x1, _y1, _z1;
  Real _x2, _y2, _z2;
  Real _int_width;

  bool _3D_sphere;
};

template<>
InputParameters validParams<BicrystalEllipticalGrainICAction>();

#endif //BICRYSTALELLIPTICALGRAINICACTION_H
