#include "GB100Test.h"
#include <cmath>

template<>
InputParameters validParams<GB100Test>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("T","Temperature in Kelvin");
  params.addParam<Real>("f0s", 0.125, "The GB energy constant ");
  params.addRequiredParam<Real>("wGB", "Diffuse GB width in nm ");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-9, "Time scale in s, where default is ns");
  params.addParam<Real>("GBMobility", -1, "GB mobility input in m^4/(J*s), that overrides the temperature dependent calculation");
  params.addParam<Real>("GBmob0", 2.5e-6, "Grain boundary mobility prefactor in m^4/(J*s)");
  params.addParam<Real>("Q", 0.23, "Grain boundary migration activation energy in eV");
  params.addParam<Real>("GBenergy", 0.708, "Grain boundary energy in J/m^2");
  params.addParam<Real>("molar_volume", 7.11e-6, "Molar volume in m^3/mol, needed for temperature gradient driving force");
  params.addRequiredCoupledVar("v", "Order parameter values");

  return params;
}


GB100Test::GB100Test(const std::string & name,
                         InputParameters parameters) :
    Material(name, parameters),
    _ncrys(coupledComponents("v")),
    _f0s(getParam<Real>("f0s")),
    _wGB(getParam<Real>("wGB")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _GBmob0(getParam<Real>("GBmob0")),
    _Q(getParam<Real>("Q")),
    _GBenergy(getParam<Real>("GBenergy")),
    _GBMobility(getParam<Real>("GBMobility")),
    _molar_vol(getParam<Real>("molar_volume")),
    _T(coupledValue("T")),
    _sigma(declareProperty<Real>("sigma")),
    _M_GB(declareProperty<Real>("M_GB")),
    _kappa(declareProperty<Real>("kappa_op")),
    _gamma(declareProperty<Real>("gamma_asymm")),
    _L(declareProperty<Real>("L")),
    _l_GB(declareProperty<Real>("l_GB")),
    _mu(declareProperty<Real>("mu")),
    _entropy_diff(declareProperty<Real>("entropy_diff")),
    _molar_volume(declareProperty<Real>("molar_volume")),
    _act_wGB(declareProperty<Real>("act_wGB")),
    _tgrad_corr_mult(declareProperty<Real>("tgrad_corr_mult")),
    _kb(8.617343e-5) //Boltzmann constant in eV/K
{
  if (_GBMobility == -1 && _GBmob0 == 0)
    mooseError("Either a value for GBMobility or for GBmob0 and Q must be provided");

  // Resize and initialize values to calculate grain boundary energy
  _v.resize(_ncrys);
  _grad_v.resize(_ncrys);

  for (unsigned int crys = 0; crys < _ncrys; ++crys)
  {
    _v[crys] = &coupledValue("v",crys);
    _grad_v[crys] = &coupledGradient("v",crys);
  }


}

void
GB100Test::computeQpProperties()
{
  Real M0 = _GBmob0;
  Real JtoeV = 6.24150974e18; // joule to eV conversion
  Real xgp, ygp; // quadrature point location to be used for the test
  std::vector<VariableValue *> vfirst; // values of the order parameter


  //Testing to see if I can modify material properties based on quadrature point location
  xgp = _q_point[_qp](0);
  ygp = _q_point[_qp](1);


//  Moose::out << (*_grad_v[1])[_qp](0) << std::endl;

// This is the quadranttest setup for _Q, _GBmob0, and _GBenergy
//  if (xgp > 500 && ygp > 500)
//    {
//    _Q = 0.23;
//    _GBmob0 = 2.6e-06;
//    _GBenergy = 0.508;
//    }
//  else
//    {
//    _Q = 0.23;
//    _GBmob0 = 2.6e-06;
//    _GBenergy = 0.908;
//    }

// This is the test to use the gradient of the order parameter to pick a gb energy

  _grad_size = std::sqrt((*_grad_v[1])[_qp]*(*_grad_v[1])[_qp]);
  _mag_grad_comp = std::sqrt((*_grad_v[1])[_qp](0)*(*_grad_v[1])[_qp](0));

  if (_grad_size <= 0.01) //so we don't have to deal with division by zero
    _grad_size = 0.01;
  //if (_mag_grad_comp <= 0.01) //so we don't have to deal with log of zero
  //  _mag_grad_comp = 0.01;

  _GBenergy = .70 + (.92-.70)*std::pow((2/M_PI * std::asin(_mag_grad_comp/_grad_size)),3.2);


  M0 *= _time_scale/(JtoeV*(_length_scale*_length_scale*_length_scale*_length_scale));//Convert to lengthscale^4/(eV*timescale);

  _sigma[_qp] = _GBenergy*JtoeV*(_length_scale*_length_scale);// eV/nm^2

  if (_GBMobility < 0)
    _M_GB[_qp] = M0*std::exp(-_Q/(_kb*_T[_qp]));
  else
    _M_GB[_qp] = _GBMobility*_time_scale/(JtoeV*(_length_scale*_length_scale*_length_scale*_length_scale)); //Convert to lengthscale^4/(eV*timescale)

  _l_GB[_qp] = _wGB; //in the length scale of the system
  _L[_qp] = 4.0/3.0*_M_GB[_qp]/_l_GB[_qp];
  _kappa[_qp] = 3.0/4.0*_sigma[_qp]*_l_GB[_qp];
  _gamma[_qp] = 1.5;
  _mu[_qp] = 3.0/4.0*1/_f0s*_sigma[_qp]/_l_GB[_qp];
  _entropy_diff[_qp] = 8.0e3*JtoeV; //J/(K mol) converted to eV(K mol)
  _molar_volume[_qp] = _molar_vol/(_length_scale*_length_scale*_length_scale); //m^3/mol converted to ls^3/mol
  _act_wGB[_qp] = 0.5e-9/_length_scale;
  _tgrad_corr_mult[_qp] = _mu[_qp]*9.0/8.0;

}
