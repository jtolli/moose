/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "SmoothEllipseIC.h"

template<>
InputParameters validParams<SmoothEllipseIC>()
{
  InputParameters params = validParams<InitialCondition>();

  params.addRequiredParam<Real>("invalue", "The variable value inside the ellipse");
  params.addRequiredParam<Real>("outvalue", "The variable value outside the ellipse");
  params.addParam<Real>("int_width", 0.0, "The interfacial width of the void surface.  Defaults to sharp interface");
  params.addParam<bool>("3D_spheres", true, "in 3D, whether the objects are spheres or columns");
    params.addRequiredParam<Real>("x1"," The x value of the 1st focal point");
    params.addRequiredParam<Real>("y1"," The y value of the 1st focal point");
    params.addParam<Real>("z1"," The z value of the 1st focal point");
    params.addRequiredParam<Real>("x2"," The x value of the 2nd focal point");
    params.addRequiredParam<Real>("y2"," The y value of the 2nd focal point");
    params.addParam<Real>("z2"," The z value of the 2nd focal point");
    params.addRequiredParam<Real>("distance","the distance between the two points and the ellipse");

  return params;
}

SmoothEllipseIC::SmoothEllipseIC(const std::string & name,
                               InputParameters parameters) :
    InitialCondition(name, parameters),
    _mesh(_fe_problem.mesh()),
    _invalue(parameters.get<Real>("invalue")),
    _outvalue(parameters.get<Real>("outvalue")),
    _int_width(parameters.get<Real>("int_width")),
    _3D_spheres(parameters.get<bool>("3D_spheres")),
    _num_dim(_3D_spheres ? 3 : 2),
    _x1(parameters.get<Real>("x1")),
    _y1(parameters.get<Real>("y1")),
    _z1(parameters.get<Real>("z1")),
    _x2(parameters.get<Real>("x2")),
    _y2(parameters.get<Real>("y2")),
    _z2(parameters.get<Real>("z2")),
    _distance(parameters.get<Real>("distance")),
    _focalPoint1(_x1, _y1, _z1),
    _focalPoint2(_x2, _y2, _z2)
{
}

Real
SmoothEllipseIC::value(const Point & p)
{
  Real value = _outvalue;
  Real val2 = 0.0;

    val2 = computeEllipseValue(p, _focalPoint1, _focalPoint2, _distance);
    if ( (val2 > value && _invalue > _outvalue) || (val2 < value && _outvalue > _invalue) )
      value = val2;


  return value;
}

RealGradient
SmoothEllipseIC::gradient(const Point & p)
{
  Point gradient = Gradient(0.0, 0.0, 0.0);
  Real value = _outvalue;
  Real val2 = 0.0;

    val2 = computeEllipseValue(p, _focalPoint1, _focalPoint2, _distance);
    if ( (val2 > value && _invalue > _outvalue) || (val2 < value && _outvalue > _invalue) )
    {
      value = val2;
      gradient = computeEllipseGradient(p, _focalPoint1, _focalPoint2, _distance);
    }

  return gradient;
}

Real
SmoothEllipseIC::computeEllipseValue(const Point & p, const Point & fp1, const Point & fp2, const Real & distance)
{
  Point l_fp1 = fp1;
  Point l_fp2 = fp2;
  Point l_p = p;
  if (!_3D_spheres) //Create 3D cylinders instead of spheres
  {
    l_p(2) = 0.0;
    l_fp1(2) = 0.0;
    l_fp2(2) = 0.0;
  }
  //Compute the distance between the current point and the center
  Real dist = _mesh.minPeriodicDistance(_var.number(), l_p, l_fp1) +
    _mesh.minPeriodicDistance(_var.number(), l_p, l_fp2);

  //Return value
  Real value = _outvalue;//Outside circle

    if (dist <= distance - _int_width/2.0) //Inside circle
      value = _invalue;
    else if (dist < distance + _int_width/2.0) //Smooth interface
    {
      Real int_pos = (dist - distance + _int_width/2.0)/_int_width;
      value = _outvalue + (_invalue - _outvalue) * (1.0 + std::cos(int_pos * libMesh::pi)) / 2.0;
    }

    return value;
}

Point
SmoothEllipseIC::computeEllipseGradient(const Point & p, const Point & fp1, const Point & fp2, const Real & distance)
{
  Point l_fp1 = fp1;
  Point l_fp2 = fp2;
  Point l_p = p;
    Point center((_x1+_x2)/2, (_y1+_y2)/2, (_z1+_z2)/2);
  if (!_3D_spheres) //Create 3D cylinders instead of spheres
  {
    l_p(2) = 0.0;
    l_fp1(2) = 0.0;
    l_fp2(2) = 0.0;
  }
  //Compute the distance between the current point and the center
    Real dist = _mesh.minPeriodicDistance(_var.number(), l_p, l_fp1) +
    _mesh.minPeriodicDistance(_var.number(), l_p, l_fp2);
    
  Real DvalueDr = 0.0;

  if (dist < distance + _int_width/2.0 && dist > distance - _int_width/2.0)
  {
    Real int_pos = (dist - distance + _int_width / 2.0) / _int_width;
    Real Dint_posDr = 1.0 / _int_width;
    DvalueDr = Dint_posDr * (_invalue - _outvalue) * (-std::sin(int_pos * libMesh::pi) * libMesh::pi) / 2.0;
  }

  //Set gradient over the smooth interface
  if (dist != 0.0)
    return _mesh.minPeriodicVector(_var.number(), center, p) * (DvalueDr / dist);
  else
    return Gradient(0.0, 0.0, 0.0);

}
