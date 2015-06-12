#ifndef GBQUADRANTTEST_H
#define GBQUADRANTTEST_H

#include "Material.h"

//Forward Declarations
class GBQuadrantTest;

template<>
InputParameters validParams<GBQuadrantTest>();

class GBQuadrantTest : public Material
{
public:
  GBQuadrantTest(const std::string & name,
              InputParameters parameters);

protected:
  virtual void computeQpProperties();
  unsigned int _ncrys;

private:
  Real _f0s;
  Real _wGB;
  Real _length_scale;
  Real _time_scale;
  Real _GBmob0;
  Real _Q;
  Real _GBenergy;
  Real _GBMobility;
  Real _molar_vol;
  Real _grad_size;
  Real _mag_grad_comp;


  VariableValue & _T;
  std::vector<VariableValue *> _v;
  std::vector<VariableGradient *> _grad_v;

  MaterialProperty<Real> & _sigma;
  MaterialProperty<Real> & _M_GB;
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _gamma;
  MaterialProperty<Real> & _L;
  MaterialProperty<Real> & _l_GB;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _entropy_diff;
  MaterialProperty<Real> & _molar_volume;
  MaterialProperty<Real> & _act_wGB;
  MaterialProperty<Real> & _tgrad_corr_mult;

  const Real _kb;
};

#endif //GBQUADRANTTEST_H
