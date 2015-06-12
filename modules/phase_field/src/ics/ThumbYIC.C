/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "ThumbYIC.h"

template<>
InputParameters validParams<ThumbYIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("ycoord", "The x coordinate of the circle center");
  params.addRequiredParam<Real>("width", "The y coordinate of the circle center");
  params.addRequiredParam<Real>("height", "The z coordinate of the circle center");
  params.addRequiredParam<Real>("invalue", "The variable value inside the circle");
  params.addRequiredParam<Real>("outvalue", "The variable value outside the circle");

  return params;
}

ThumbYIC::ThumbYIC(const std::string & name,
                 InputParameters parameters) :
    InitialCondition(name, parameters),
    _ycoord(parameters.get<Real>("ycoord")),
    _width(parameters.get<Real>("width")),
    _height(parameters.get<Real>("height")),
    _invalue(parameters.get<Real>("invalue")),
    _outvalue(parameters.get<Real>("outvalue"))
{}

Real
ThumbYIC::value(const Point & p)
{
  Real value = 0.0;

  if (p(0) > _height)
  {
    Real rad = 0.0;
    Point center(_height, _ycoord, 0.0);
    for (unsigned int i = 0; i < 2; ++i)
      rad += (p(i) - center(i)) * (p(i) - center(i));

    rad = sqrt(rad);

    if (rad <= _width/2.0)
      value = _invalue;
    else
      value = _outvalue;
  }
  else
  {
    if (p(1) > _ycoord - _width/2.0 && p(1) < _ycoord + _width/2.0)
      value = _invalue;
    else
      value = _outvalue;
  }

  return value;
}
