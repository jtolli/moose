/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "TwoThumbIC.h"

template<>
InputParameters validParams<TwoThumbIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("xcoord", "The x coordinate of the first thumb's center");
  params.addRequiredParam<Real>("ycoord", "The y coordinate of the second thumb's center");
  params.addRequiredParam<Real>("width", "the width of the thumbs");
  params.addRequiredParam<Real>("height", "The height of the thumbs");
  params.addRequiredParam<Real>("invalue", "The variable value inside the thumbs");
  params.addRequiredParam<Real>("outvalue", "The variable value outside the thumbs");

  return params;
}

TwoThumbIC::TwoThumbIC(const std::string & name,
                 InputParameters parameters) :
    InitialCondition(name, parameters),
    _xcoord(parameters.get<Real>("xcoord")),
    _ycoord(parameters.get<Real>("ycoord")),
    _width(parameters.get<Real>("width")),
    _height(parameters.get<Real>("height")),
    _invalue(parameters.get<Real>("invalue")),
    _outvalue(parameters.get<Real>("outvalue"))
{}

Real
TwoThumbIC::value(const Point & p)
{
  Real value = _outvalue;

  if (p(1) > _height)
  {
    Real rad = 0.0;
    Point center(_xcoord, _height, 0.0);
    for (unsigned int i = 0; i < 2; ++i)
      rad += (p(i) - center(i)) * (p(i) - center(i));

    rad = sqrt(rad);

    if (rad <= _width/2.0)
      value = _invalue;
  }
  else
  {
    if (p(0) > _xcoord - _width/2.0 && p(0) < _xcoord + _width/2.0)
      value = _invalue;
  }

  if (p(0) > _height)
  {
    Real rad = 0.0;
    Point center(_height, _ycoord, 0.0);
    for (unsigned int i = 0; i < 2; ++i)
      rad += (p(i) - center(i)) * (p(i) - center(i));

    rad = sqrt(rad);

    if (rad <= _width/2.0)
      value = _invalue;
  }

  else
  {
    if (p(1) > _ycoord - _width/2.0 && p(1) < _ycoord + _width/2.0)
      value = _invalue;
  }

  return value;
}
