/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "BicrystalEllipticalGrainICAction.h"
#include "Factory.h"
#include "Parser.h"
#include "FEProblem.h"

#include <sstream>
#include <stdexcept>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"

const Real BicrystalEllipticalGrainICAction::_abs_zero_tol = 1e-12;

template<>
InputParameters validParams<BicrystalEllipticalGrainICAction>()
{
  InputParameters params = validParams<Action>();

  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addRequiredParam<unsigned int>("op_num", "Number of grains, should be 2");
  params.addRequiredParam<Real>("distance", "distance between the two points and the ellipse");
  params.addRequiredParam<Real>("x1", "The x coordinate of the 1st ellipse grain focal point");
  params.addRequiredParam<Real>("y1", "The y coordinate of the 1st ellipse grain  focal point");
  params.addParam<Real>("z1", 0.0, "The z coordinate of the 1st ellipse grain focal point");
  params.addRequiredParam<Real>("x2", "The x coordinate of the 2nd ellipse grain focal point");
  params.addRequiredParam<Real>("y2", "The y coordinate of the 2st ellipse grain  focal point");
  params.addParam<Real>("z2", 0.0, "The z coordinate of the 2nd ellipse grain focal point");
  params.addParam<Real>("int_width", 0.0, "The interfacial width of the void surface.  Defaults to sharp interface");
  params.addParam<bool>("3D_sphere", true, "in 3D, whether the smaller grain is a spheres or columnar grain");
  return params;
}

BicrystalEllipticalGrainICAction::BicrystalEllipticalGrainICAction(const std::string & name, InputParameters params) :
    Action(name, params),
    _var_name_base(getParam<std::string>("var_name_base")),
    _op_num(getParam<unsigned int>("op_num")),
    _distance(getParam<Real>("distance")),
    _x1(getParam<Real>("x1")),
    _y1(getParam<Real>("y1")),
    _z1(getParam<Real>("z1")),
    _x2(getParam<Real>("x2")),
    _y2(getParam<Real>("y2")),
    _z2(getParam<Real>("z2")),
    _int_width(getParam<Real>("int_width")),
    _3D_sphere(getParam<bool>("3D_sphere"))
{
  if (_op_num != 2)
    mooseError("op_num must equal 2 for bicrystal ICs");
}

void
BicrystalEllipticalGrainICAction::act()
{
#ifdef DEBUG
  Moose::err << "Inside the BicrystalEllipticalGrainICAction Object\n";
#endif

  // Loop through the number of order parameters
  for (unsigned int op = 0; op < _op_num; op++)
  {
    //Create variable names
    std::string var_name = _var_name_base;
    std::stringstream out;
    out << op;
    var_name.append(out.str());

    //Set parameters for SmoothEllipseIC
    InputParameters poly_params = _factory.getValidParams("SmoothEllipseIC");
    poly_params.set<VariableName>("variable") = var_name;
    poly_params.set<Real>("x1") = _x1;
    poly_params.set<Real>("y1") = _y1;
    poly_params.set<Real>("z1") = _z1;
    poly_params.set<Real>("x2") = _x2;
    poly_params.set<Real>("y2") = _y2;
    poly_params.set<Real>("z2") = _z2;
    poly_params.set<Real>("distance") = _distance;
    poly_params.set<Real>("int_width") = _int_width;
    poly_params.set<bool>("3D_spheres") = _3D_sphere;
    if (op == 0)
    {
      //Values for ellipse grain
      poly_params.set<Real>("invalue") = 1.0;
      poly_params.set<Real>("outvalue") = 0.0;
    }
    else
    {
      //Values for matrix grain
      poly_params.set<Real>("invalue") = 0.0;
      poly_params.set<Real>("outvalue") = 1.0;
    }

    //Add initial condition
    _problem->addInitialCondition("SmoothEllipseIC", "InitialCondition", poly_params);
  }
}
