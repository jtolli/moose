#include "ACGrGrAnis.h"

template<>
InputParameters validParams<ACGrGrAnis>()
{
  InputParameters params = validParams<ACBulk>();
  params.addClassDescription("Grain-Boundary model poly crystaline interface Allen-Cahn Kernel");
  params.addRequiredCoupledVar("v", "Array of coupled variable names");
  params.addCoupledVar("T", "temperature");
  params.addRequiredParam<unsigned int>("op","The order parameter number this is acting on");
  params.addRequiredParam<FileName>("Anisotropic_GB_file_name", "Name of the file containing: 1)GB mobility prefactor; 2) GB migration activation energy; 3)GB energy");
  return params;
}

ACGrGrAnis::ACGrGrAnis(const InputParameters & parameters) :
    ACBulk(parameters),
    _Anisotropic_GB_file_name(getParam<FileName>("Anisotropic_GB_file_name")),
    _mu(getMaterialProperty<Real>("mu")),
    _gamma(getMaterialProperty<Real>("gamma_asymm")),
    _tgrad_corr_mult(getMaterialProperty<Real>("tgrad_corr_mult")),
    _l_GB(getMaterialProperty<Real>("l_GB")),
    _has_T(isCoupled("T")),
    _grad_T(_has_T ? &coupledGradient("T") : NULL),
    _op(getParam<unsigned int>("op")),
    _ncrys(coupledComponents("v"))
{
  // Initialize values for crystals
  _vals.resize(_ncrys);
  _grad_vals.resize(_ncrys);

  // reshape vectors
  _sigma.resize(_ncrys);

  for (unsigned int crys = 0; crys < _ncrys; ++crys)
  {
    // Initialize variables
    _vals[crys] = &coupledValue("v", crys);
    _grad_vals[crys] = &coupledGradient("v", crys);

    _sigma[crys].resize(_ncrys);

  }

  // Read in data from "Anisotropic_GB_file_name"
  std::ifstream inFile(_Anisotropic_GB_file_name.c_str());

  if (!inFile)
    mooseError("Can't open GB anisotropy input file");

  for (unsigned int i = 0; i < 2; ++i)
    inFile.ignore(255, '\n'); // ignore line

  Real data;
  for (unsigned int i = 0; i < 3*_ncrys; ++i)
  {
    std::vector<Real> row; // create an empty row of double values
    for (unsigned int j = 0; j < _ncrys; ++j)
    {
      inFile >> data;
      row.push_back(data);
    }

    if (i < _ncrys)
      _sigma[i] = row; // unit: J/m^2
  }

  inFile.close();
}

Real
ACGrGrAnis::computeDFDOP(PFFunctionType type)
{
  Real SumEtaj = 0.0;
  Real SumEtaij = 0.0;
  Real SumEtaSigmaj = 0.0;
  Real SumEtaSigmaij = 0.0;
  Real Sum2Eta4Eta = (_u[_qp] * _u[_qp] * _u[_qp] * _u[_qp]) / 4 * (_u[_qp] * _u[_qp]) / 2;
  Real SumGradEta = _grad_u[_qp] * _grad_u[_qp];
  for (unsigned int i = 0; i < _ncrys; ++i)
  {
    SumEtaj += ((*_vals[i])[_qp] * (*_vals[i])[_qp]); //Sum all other order parameters
    SumEtaij += ((*_vals[i])[_qp] * (*_vals[i])[_qp] * _u[_qp] * _u[_qp]);
    SumEtaSigmaj += ((*_vals[i])[_qp] * (*_vals[i])[_qp]) * _sigma[i][_op];
    SumEtaSigmaij += ((*_vals[i])[_qp] * (*_vals[i])[_qp] * _u[_qp] * _u[_qp] *_sigma[i][_op]);
    Sum2Eta4Eta += ((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[i])[_qp]) / 4 * ((*_vals[i])[_qp] * (*_vals[i])[_qp]) / 2;
    SumGradEta += ((*_grad_vals[i])[_qp] * (*_grad_vals[i])[_qp]);

    for (unsigned int j = i; j < _ncrys; ++j)
    {
      SumEtaij += ((*_vals[i])[_qp]*(*_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp]);
      SumEtaSigmaij += ((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp]*(*_vals[j])[_qp]) * _sigma[i][j];
    }
  }

  Real Dsigma_Deta = 2 * _u[_qp] * (SumEtaj * SumEtaSigmaij - SumEtaSigmaj * SumEtaSigmaij) / (SumEtaij * SumEtaij);

  Real tgrad_correction = 0.0;

  //Calcualte either the residual or jacobian of the grain growth free energy
  switch (type)
  {
    case Residual:
      if (_has_T)
        tgrad_correction = _tgrad_corr_mult[_qp]*_grad_u[_qp]*(*_grad_T)[_qp];
      return Dsigma_Deta * (6.0 / _l_GB[_qp] * (Sum2Eta4Eta + _gamma[_qp] * SumEtaij + 0.25) + 0.75 * _l_GB[_qp] *SumGradEta);

    case Jacobian:
      if (_has_T)
        tgrad_correction = _tgrad_corr_mult[_qp]*_grad_phi[_j][_qp]*(*_grad_T)[_qp];
      return _phi[_j][_qp]*Dsigma_Deta * (6.0 / _l_GB[_qp] * (Sum2Eta4Eta + _gamma[_qp] * SumEtaij + 0.25) + 0.75 * _l_GB[_qp] *SumGradEta);
  }

  mooseError("Invalid type passed in");
}

Real
ACGrGrAnis::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}
