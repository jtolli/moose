/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef GBDISCOUNTINOUSENERGY_H
#define GBDISCOUNTINOUSENERGY_H

#include "Material.h"

//Forward Declarations
class GBDiscountinousEnergy;

template<>
InputParameters validParams<GBDiscountinousEnergy>();

class GBDiscountinousEnergy : public Material
{
public:
  GBDiscountinousEnergy(const std::string & name,
              InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
  Real _f0s;
  Real _wGB;
  Real _length_scale;
  Real _time_scale;
  Real _GBmob0;
  Real _Q;
  Real _GBenergy;
  Real _GBenergyIn;
  Real _GBenergyOut;
  Real _GBMobility;
  Real _molar_vol;
  //Real _ncrys;
  Real _cut_off;

  VariableValue & _T;
  VariableValue & _v;

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

#endif //GBEVOLUTION_H
