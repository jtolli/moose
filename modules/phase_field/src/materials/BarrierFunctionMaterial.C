/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "BarrierFunctionMaterial.h"

template<>
InputParameters validParams<BarrierFunctionMaterial>()
{
  InputParameters params = validParams<OrderParameterFunctionMaterial>();
  params.addClassDescription("Helper material to provide g(eta) and its derivative in a polynomial.\nSIMPLE: eta^2*(1-eta)^2\nLOW: eta*(1-eta)");
  MooseEnum g_order("SIMPLE=0 LOW", "SIMPLE");
  params.addParam<MooseEnum>("g_order", g_order, "Polynomial order of the barrier function g(eta)");
  params.addParam<bool>("well_only", false, "Make the g zero in [0:1] so it only contributes to enforcing the eta range and not to the phase transformation berrier.");
  params.set<std::string>("function_name") = std::string("g");
  return params;
}

BarrierFunctionMaterial::BarrierFunctionMaterial(const std::string & name,
                                                 InputParameters parameters) :
    OrderParameterFunctionMaterial(name, parameters),
    _g_order(getParam<MooseEnum>("g_order")),
    _well_only(getParam<bool>("well_only"))
{
}

void
BarrierFunctionMaterial::computeQpProperties()
{
  const Real n = _eta[_qp];

  if (_well_only && n >= 0.0 && n <= 1.0) {
    _prop_f[_qp] = 0.0;
    _prop_df[_qp] = 0.0;
    _prop_d2f[_qp] = 0.0;
    return;
  }

  switch (_g_order)
  {
    case 0: // SIMPLE
      _prop_f[_qp]   =  n*n * (1.0 - n) * (1.0 - n);
      _prop_df[_qp]  =  2.0 * n * (n - 1.0) * (2.0 * n - 1.0);
      _prop_d2f[_qp] = 12.0 * (n * n - n) + 2.0;
      break;

    case 1: // LOW
      _prop_f[_qp]   = n * (1.0 - n);
      _prop_df[_qp]  = 1.0 - 2.0 * n;
      _prop_d2f[_qp] = - 2.0;
      break;

    default:
      mooseError("Internal error");
  }
}
