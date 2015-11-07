#include "ACBulk.h"

#ifndef ACGRGRANIS_H
#define ACGRGRANIS_H

//Forward Declarations
class ACGrGrAnis;

template<>
InputParameters validParams<ACGrGrAnis>();

/**
 * This kernel calculates the residual for grain growth.
 * It calculates the residual of the ith order parameter, and the values of
 * all other order parameters are coupled variables and are stored in vals.
 */
class ACGrGrAnis : public ACBulk
{
public:
  ACGrGrAnis(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
  std::vector<VariableValue *> _vals;
  std::vector<unsigned int> _vals_var;
  std::string _Anisotropic_GB_file_name;

  std::vector<std::vector<Real> > _sigma;


  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _gamma;
  const MaterialProperty<Real> & _tgrad_corr_mult;
  const MaterialProperty<Real> & _l_GB;

  bool _has_T;
  VariableGradient * _grad_T;

  unsigned int _op;
  unsigned int _ncrys;
  // Real _gamma;
  std::vector<VariableGradient *> _grad_vals;
};

#endif //ACGrGrAnis_H
