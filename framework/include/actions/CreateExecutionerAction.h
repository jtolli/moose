/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef CREATEEXECUTIONERACTION_H
#define CREATEEXECUTIONERACTION_H

#include "MooseObjectAction.h"

class CreateExecutionerAction;

template<>
InputParameters validParams<CreateExecutionerAction>();


class CreateExecutionerAction : public MooseObjectAction
{
public:
  CreateExecutionerAction(const std::string & name, InputParameters params);

  virtual void act();

  static MultiMooseEnum getCommonPetscOptions();

  static void populateCommonExecutionerParams(InputParameters & params);
  static void    storeCommonExecutionerParams(FEProblem & problem, InputParameters & params);
};

#endif // CREATEEXECUTIONERACTION_H
