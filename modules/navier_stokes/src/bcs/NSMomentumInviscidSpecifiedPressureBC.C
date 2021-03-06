/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "NSMomentumInviscidSpecifiedPressureBC.h"

template<>
InputParameters validParams<NSMomentumInviscidSpecifiedPressureBC>()
{
  InputParameters params = validParams<NSMomentumInviscidBC>();

  // Required parameters.
  params.addRequiredParam<Real>("specified_pressure", "The specified pressure for this boundary");

  return params;
}




NSMomentumInviscidSpecifiedPressureBC::NSMomentumInviscidSpecifiedPressureBC(const std::string & name, InputParameters parameters)
    : NSMomentumInviscidBC(name, parameters),

      // Parameters to be specified in input file block...
      _specified_pressure(getParam<Real>("specified_pressure"))
{
}




Real NSMomentumInviscidSpecifiedPressureBC::computeQpResidual()
{
  // Velocity vector object
  RealVectorValue vel(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  // Velocity vector dotted with normal
  Real u_dot_n = vel * _normals[_qp];

  // The current value of the vector (rho*u)(u.n)
  RealVectorValue rhou_udotn = u_dot_n * _rho[_qp] * vel;

  return
    this->pressure_qp_residual(_specified_pressure) +
    this->convective_qp_residual( rhou_udotn(_component) );
}



Real NSMomentumInviscidSpecifiedPressureBC::computeQpJacobian()
{
  // There is no Jacobian for the pressure term when the pressure is specified,
  // so all we have left is the convective part.  The on-diagonal variable number
  // is _component+1
  return this->convective_qp_jacobian(_component+1);
}



Real NSMomentumInviscidSpecifiedPressureBC::computeQpOffDiagJacobian(unsigned jvar)
{
  return this->convective_qp_jacobian( this->map_var_number(jvar) );
}
