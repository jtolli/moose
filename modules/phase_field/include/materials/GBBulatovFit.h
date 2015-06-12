#ifndef GBBULATOVFIT_H
#define GBBULATOVFIT_H

#include "Material.h"
#include "RotationMatrix.h"
#include "RotationTensor.h"

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

//Forward Declarations
class GBBulatovFit;

template<>
InputParameters validParams<GBBulatovFit>();

class GBBulatovFit : public Material
{
public:
  GBBulatovFit(const std::string & name,
              InputParameters parameters);

protected:
  virtual void computeQpProperties();
  unsigned int _ncrys;

private:
  Real gB5DOF(RealTensorValue _P, RealTensorValue _S);
  void distancesToSet(RealTensorValue _P, RealTensorValue _S, unsigned int axis,
     std::vector<std::vector<Real> > & geom);
  void quat2Mat(std::vector<Real> q, RealTensorValue & M);
  void mat2Quat(RealTensorValue M, std::vector<Real> & q);
  void makeParVec();
  Real weightedMeanEnergy(std::vector<std::vector <Real> > geom100,
    std::vector<std::vector <Real> > geom110,
    std::vector<std::vector <Real> > geom111);
  void set100(std::vector<std::vector<Real> > geom100, std::vector<Real> & e100);
  void set110(std::vector<std::vector<Real> > geom110, std::vector<Real> & e110);
  void set111(std::vector<std::vector<Real> > geom111, std::vector<Real> & e111);
  void twist100(std::vector<Real> ksi, std::vector<Real> & entwist);
  void twist110(std::vector<Real> ksi, std::vector<Real> & entwist);
  void twist111(std::vector<Real> ksi, std::vector<Real> & entwist);
  void aTGB100(std::vector<Real> eta, std::vector<Real> ksi,
    std::vector<Real> & entilt);
  void aTGB110(std::vector<Real> eta, std::vector<Real> ksi,
    std::vector<Real> & entilt);
  void aTGB111(std::vector<Real> eta, std::vector<Real> ksi,
    std::vector<Real> & entilt);
  void sTGB100(std::vector<Real> ksi, std::vector<Real> & _en);
  void sTGB110(std::vector<Real> ksi, std::vector<Real> & _en);
  Real rSW(Real theta, Real thetaMin, Real thetaMax, Real a);


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
  std::string _material;
  std::vector<RealVectorValue> _euler_angles;
  std::vector<RealTensorValue> _rotation_tensors;
  std::vector<Real> _par_vec;

  VariableValue & _bnds;
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

#endif //GBBULATOVFIT_H
