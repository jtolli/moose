#include "GBBulatovFit.h"

template<>
InputParameters validParams<GBBulatovFit>()
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
  params.addParam<std::string>("material", "Cu", "Which material parameters to use");
  params.addRequiredCoupledVar("v", "Order parameter values");
  params.addRequiredCoupledVar("bnds", "Grain Boundary locations");
  //params.addRequiredCoupledVar("euler_angles","euler angles for the grains");
  return params;
}


GBBulatovFit::GBBulatovFit(const std::string & name,
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
  _material(getParam<std::string>("material")),
  _bnds(coupledValue("bnds")),
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

    _par_vec.resize(43);
    makeParVec();
    // Resize and initialize values to calculate grain boundary energy
    _v.resize(_ncrys);
    _grad_v.resize(_ncrys);
    _euler_angles.resize(2);
    _rotation_tensors.resize(_ncrys);
    _euler_angles[0](0) = 0;
    _euler_angles[0](1) = 0;
    _euler_angles[0](2) = 0;
    _euler_angles[1](0) = 0;
    _euler_angles[1](1) = 8;
    _euler_angles[1](2) = 0;
    for (unsigned int crys = 0; crys < _ncrys; ++crys)
    {
      _v[crys] = &coupledValue("v",crys);
      _grad_v[crys] = &coupledGradient("v",crys);
      _rotation_tensors[crys] = RotationTensor(_euler_angles[crys]);
    }

  }

  void
  GBBulatovFit::computeQpProperties()
  {
    Real M0 = _GBmob0;
    Real JtoeV = 6.24150974e18; // joule to eV conversion
    std::vector<VariableValue *> vfirst; // values of the order parameter
    RealTensorValue _rotation_matrix; // rotates _gbdir to <100>
    RealVectorValue _gbdir; // the normal to the grain boundary

    std::vector<Real> _op_in_bnd; //which order parameters are in this boundary
    _op_in_bnd.resize(0);
    for (unsigned int crys = 0; crys < _ncrys; crys++)
    {
      if ((*_v[crys])[_qp] > 0.001)
      _op_in_bnd.push_back(crys);
    }

    if (_op_in_bnd.size() == 2)
    {
      _grad_size = std::sqrt((*_grad_v[_op_in_bnd[0]])[_qp]*(*_grad_v[_op_in_bnd[0]])[_qp]);
      _gbdir = (*_grad_v[_op_in_bnd[0]])[_qp]/_grad_size; //RealGradient and RealVectorValue are really the same type


      RealVectorValue axis100(1,0,0);
      _rotation_matrix = RotationMatrix::rotVec1ToVec2(_gbdir, axis100); //currently there is an issue with this not being a number

      _GBenergy = gB5DOF(_rotation_matrix*_rotation_tensors[_op_in_bnd[0]],
        _rotation_matrix*_rotation_tensors[_op_in_bnd[1]]);
      }
      else
      {
        _GBenergy = getParam<Real>("GBenergy");
      }

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

    void
    GBBulatovFit::makeParVec()
    {
      /* makeParVec()

      Creates a 43-parameter vector as used by weightedMeanEnergy

      _material is a string "Al", "Ni", "Au", or "Cu", which                                                 then sets all other parameters automatically.

      AlCuparameter:  Position on the Al-Cu axis, where 0 is Al and 1 is Cu
      (this parameter is capital Phi in the journal article).
      This is related to eSF/eRGB, where eSF is the stacking-fault energy.

      eRGB:  Energy of a "random" grain boundary in J/m^2


      par42Al:  The 42 dimensionless parameters for Al

      par42Cu:  The 42 dimensionless parameters for Cu

      Note a majority of the entries for par42Al and par42Cu are normally identical.

      par42Al and par42Cu are the values found by numerical fitting to the 388*4 boundaries.
      */

      Real par42Al[42] = {0.405204179289160, 0.738862004021890, 0.351631012630026,
        2.40065811939667, 1.34694439281655, 0.352260396651516, 0.602137375062785,
        1.58082498976078, 0.596442399566661, 1.30981422643602, 3.21443408257354,
        0.893016409093743, 0.835332505166333, 0.933176738717594, 0.896076948651935,
        0.775053293192055, 0.391719619979054, 0.782601780600192, 0.678572601273508,
        1.14716256515278, 0.529386201144101, 0.909044736601838, 0.664018011430602,
        0.597206897283586, 0.200371750006251, 0.826325891814124, 0.111228512469435,
        0.664039563157148, 0.241537262980083, 0.736315075146365, 0.514591177241156,
        1.73804335876546, 3.04687038671309, 1.48989831680317, 0.664965104218438,
        0.495035051289975, 0.495402996460658, 0.468878130180681, 0.836548944799803,
        0.619285521065571, 0.844685390948170, 1.02295427618256};
        Real par42Cu[42] = {0.405204179289160, 0.738862004021890, 0.351631012630026,
          2.40065811939667, 1.34694439281655, 3.37892632736175, 0.602137375062785,
          1.58082498976078, 0.710489498577995, 0.737834049784765, 3.21443408257354,
          0.893016409093743, 0.835332505166333, 0.933176738717594, 0.896076948651935,
          0.775053293192055, 0.509781056492307, 0.782601780600192, 0.762160812499734,
          1.10473084066580, 0.529386201144101, 0.909044736601838, 0.664018011430602,
          0.597206897283586, 0.200371750006251, 0.826325891814124, 0.0226010533470218,
          0.664039563157148, 0.297920289861751, 0.666383447163744, 0.514591177241156,
          1.73804335876546, 2.69805148576400, 1.95956771207484, 0.948894352912787,
          0.495035051289975, 0.301975031994664, 0.574050577702240, 0.836548944799803,
          0.619285521065571, 0.844685390948170, 0.0491040633104212};

          Real AlCuParameter;
          Real eRGB;

          if (!_material.compare("Ni"))
          {
            eRGB = 1.44532834613925;
            AlCuParameter = 0.767911805073948;
          }
          else if (!_material.compare("Al"))
          {
            eRGB = 0.547128733614891;
            AlCuParameter = 0;
          }
          else if (!_material.compare("Au"))
          {
            eRGB = 0.529912885175204;
            AlCuParameter = 0.784289766313152;
          }
          else if (!_material.compare("Cu"))
          {
            eRGB = 1.03669431227427;
            AlCuParameter = 1;
          }

          _par_vec[0]=eRGB;

          for (unsigned int i = 1; i < 43; i++)
          {
            _par_vec[i] = par42Al[i-1]+AlCuParameter*(par42Cu[i-1]-par42Al[i-1]);
          }


        }

        Real
        GBBulatovFit::gB5DOF(RealTensorValue _P, RealTensorValue _S)
        {
          /* This function calculates the energy */
          std::vector<std::vector<Real> > geom100;
          std::vector<std::vector<Real> > geom110;
          std::vector<std::vector<Real> > geom111;
          // _P and _S are the rotation Matrices for the top and bottom grains
          distancesToSet(_P, _S, 100, geom100);
          distancesToSet(_P, _S, 110, geom110);
          distancesToSet(_P, _S, 111, geom111);
          return weightedMeanEnergy(geom100, geom110, geom111);
        }

        void
        GBBulatovFit::distancesToSet(RealTensorValue _P, RealTensorValue _S, unsigned int axis, std::vector<std::vector<Real> > & geom)
        {
          /* distancesToSet(P,S,axis,geom)

          Calculates the geometry parameters (geom) for a given grain boundary relative to
          a given set of axes.

          P and S are rotation matrices giving the orientations of the two grains.
          The grain boundary normal is fixed at [1,0,0].

          the S rotation matrix is Q in Bulatov's Code

          axis is one of 100, 110, or 111

          Result geom is a 4xn vector array of the 4 vectors distance, ksi, eta, and
          phi. It keeps all of the hits up to a distance of 1.

          distance is 2*sin(delta/2) where delta is the angle of closest approach
          between a misorientation axis and one of the axes.  Note there are 24
          representations of the rotations and 3, 6, or 4 equivalent high-symmetry
          axes, so it calculates as many as 144 distances.  But only ones below
          the cutoff of 1 are kept.

          Once it's picked the closest approximation to the boundary for a given
          axis and coset element, it finds the parameters ksi, eta, phi defining
          that idealized boundary (since the axis is defined, it's a 3-space).

          These are:
          phi, the angle between the rotation axis and the boundary plane normal
          (taken as the mean of the normals represented in the two actual grain
          orientations, which works when dismax is less than 1)

          ksi, the misorientation angle

          eta, a parameter giving the second axis of the boundary plane normal in
          terms of specified directions ('dirs') perpendicular to each
          high-symmetry axis.
          */

          Real dismax = 0.999999;  //force the distance to be strictly less than one, allowing for roundoff
          Real s2 = 1/std::sqrt(2);
          Real s3 = 1/std::sqrt(3);

          //determine the number of axes for a given rotation axis
          unsigned int naxes;
          if (axis == 110)
          {
            naxes = 6;
          }
          else if (axis == 111)
          {
            naxes = 4;
          }
          else
          {
            naxes = 3;
          }
          std::vector<RealVectorValue> axes(naxes);
          std::vector<RealVectorValue> dirs(naxes);
          if (axis == 110)
          {
            //define the 110 axes, normalize
            axes[0](0) = s2;
            axes[1](0) = s2;
            axes[2](0) = s2;
            axes[3](0) = s2;
            axes[4](0) = 0;
            axes[5](0) = 0;
            axes[0](1) = s2;
            axes[1](1) = -s2;
            axes[2](1) = 0;
            axes[3](1) = 0;
            axes[4](1) = s2;
            axes[5](1) = s2;
            axes[0](2) = 0;
            axes[1](2) = 0;
            axes[2](2) = s2;
            axes[3](2) = -s2;
            axes[4](2) = s2;
            axes[5](2) = -s2;
            //define the crystal direction perpendicular to each rotation axis.
            //The formalism demands that this be an axis of at least two-fold symmetry.
            dirs[0](0) = 0;
            dirs[1](0) = 0;
            dirs[2](0) = 0;
            dirs[3](0) = 0;
            dirs[4](0) = 1;
            dirs[5](0) = 1;
            dirs[0](0) = 0;
            dirs[1](1) = 0;
            dirs[2](1) = 1;
            dirs[3](1) = 1;
            dirs[4](1) = 0;
            dirs[5](1) = 0;
            dirs[0](2) = 1;
            dirs[1](2) = 1;
            dirs[2](2) = 0;
            dirs[3](2) = 0;
            dirs[4](2) = 0;
            dirs[5](2) = 0;
          }
          else if (axis == 111)
          {
            //define the 111 axes, normalize.
            axes[0](0) = s3;
            axes[1](0) = s3;
            axes[2](0) = -s3;
            axes[3](0) = -s3;
            axes[0](1) = s3;
            axes[1](1) = -s3;
            axes[2](1) = s3;
            axes[3](1) = -s3;
            axes[0](2) = s3;
            axes[1](2) = -s3;
            axes[2](2) = -s3;
            axes[3](2) = s3;
            //define crystal direction
            dirs[0](0) = s2;
            dirs[1](0) = s2;
            dirs[2](0) = s2;
            dirs[3](0) = s2;
            dirs[0](1) = -s2;
            dirs[1](1) = s2;
            dirs[2](1) = s2;
            dirs[3](1) = -s2;
            dirs[0](2) = 0;
            dirs[1](2) = 0;
            dirs[2](2) = 0;
            dirs[3](2) = 0;
          }
          else
          {
            //define 100 axes, normalize
            axes[0](0) = 1;
            axes[1](0) = 0;
            axes[2](0) = 0;
            axes[0](1) = 0;
            axes[1](1) = 1;
            axes[2](1) = 0;
            axes[0](2) = 0;
            axes[1](2) = 0;
            axes[2](2) = 1;
            //define crystal direction
            dirs[0](0) = 0;
            dirs[1](0) = 0;
            dirs[2](0) = 1;
            dirs[0](1) = 1;
            dirs[1](1) = 0;
            dirs[2](1) = 0;
            dirs[0](2) = 0;
            dirs[1](2) = 1;
            dirs[2](2) = 0;
          }
          Real period = libMesh::pi*naxes/6;

          //define symmetry operators
          RealTensorValue rotX90(1, 0, 0, 0, 0, 1, 0, -1, 0);     //Rotation by +90 degrees around X axis.  Check to see if it works correctly since Bulatov vectors were treated as rows.
          RealTensorValue rotY90(0, 0, -1, 0, 1, 0, 1, 0, 0);     //Rotation by +90 degrees around Y axis
          RealTensorValue rotZ90(0, 1, 0, -1, 0, 0, 0, 0, 1);     //Rotation by +90 degrees around Z axis
          RealTensorValue rotZ90m(0, -1, 0, 1, 0, 0, 0, 0, 1);    //Rotation by -90 degrees around Z axis

          /* Create 24 symmetry equivalent variants of Q
          This is the coset appropriate for the rotation convention where Q*P
          is the misorientation represented in the grain frame.  If you're
          getting odd results, e.g. misorientations that you know are CSL are
          coming out entirely wrong, you may be using the opposite convention;
          try replacing _P and _S with P' and S'.*/

          std::vector< RealTensorValue > V(24);   //initialize 24x3x3 array to hold variants of Q
          V[0] = _S;  //input Q into slot 0 of array

          RealTensorValue Temp;
          RealTensorValue Temp2;
          Temp = _S;

          for (unsigned int i = 1; i < 4; i++)   //Rotate three times around X by +90 degrees
          {
            Temp2 = Temp * rotX90;
            V[i] = Temp2;
            Temp = Temp2;
          }
          for (unsigned int i = 4; i < 16; i++)   //Rotate three times around Y by +90 degrees
          {
            Temp = V[i-4];
            Temp2 = Temp*rotY90;
            V[i] = Temp2;
          }
          for (unsigned int i = 16; i < 20; i++)   //Rotate around Z by +90 degrees
          {
            Temp= V[i-16];
            Temp2 = Temp * rotZ90;
            V[i] = Temp2;
          }
          for (unsigned int i = 20; i < 24; i++)   //Rotate around Z by -90 degrees
          {
            Temp = V[i-20];
            Temp2 = Temp * rotZ90m;
            V[i] = Temp2;
          }


          //Preallocate all parameter lists at their maximum possible sizes
          //Redundant representations will be removed at the end
          std::vector<Real> distances(24*naxes);
          std::vector<Real> phis(24*naxes);
          std::vector<Real> ksis(24*naxes);
          std::vector<Real> etas(24*naxes);

          unsigned int thisindex = 0;   //Number of hits found so far

          //Step through all combinations of symmetrically-equivalent axes and coset elements V
          for (unsigned int i = 0; i < naxes; i++)
          {
            RealVectorValue ax = axes[i];   //high-symmetry axis

            //pivot vector used to partition the rotation around axis i
            RealVectorValue dir = dirs[i];

            //Completing the orthonormal coordinate set.
            //theta1 and theta2 are defined in the plane spanned by (dir,dir2)
            RealVectorValue dir2 = ax.cross(dir);

            RealTensorValue _R;
            for (unsigned int j = 0; j < 24; j++)   //for each symmetry-related variant of the second grain
            {
              _S = V[j];   //get matrix from V array
              _R = _S.transpose() * _P;
              //This rotates any vector in cube P into a vector in cube S

              std::vector<Real> q(4,0);

              mat2Quat(_R,q);   //Calculation from here on out is much easier with quaternions.

              Real lq = std::sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
              RealVectorValue axi;
              axi(0) = q[1]/lq;
              axi(1) = q[2]/lq;
              axi(2) = q[3]/lq;   //Normalized roation axis
              Real psi = 2*std::acos(q[0]);   //Rotation angle

              Real dot_p = axi*ax;

              //Compute rotational distance from boundary P/S to the rotation set i.
              //This formula produces 2*sin(delta/2), where delta is the angle of closest approach.
              Real dis = 2*std::sqrt(std::abs(1-dot_p*dot_p))*std::sin(psi/2);

              if (dis < dismax)
              {
                //angle of rotation about ax that most closely approximates R
                Real theta = 2*std::atan(dot_p*std::tan(psi/2));
                //Compute the normal of the best-fitting GB in grain 1
                RealVectorValue n1;
                n1(0) = _P(0,0);
                n1(1) = _P(0,1);
                n1(2) = _P(0,2);
                RealVectorValue n2;
                n2(0) = _S(0,0);
                n2(1) = _S(0,1);
                n2(2) = _S(0,2);

                RealTensorValue RA;
                std::vector<Real> q2(4);
                q2[0] = std::cos(theta/2);
                q2[1] = std::sin(theta/2)*ax(0);
                q2[2] = std::sin(theta/2)*ax(1);
                q2[3] = std::sin(theta/2)*ax(2);
                quat2Mat(q2,RA);   //Rotation matrix that most closely approximates R

                //from this point on we're dealing with the idealized rotation RA, not the original R
                RealTensorValue RAp;
                RAp = RA.transpose();
                RealVectorValue m1 = n1 + RAp*n2;

                if (m1.size() >= 0.000001)   //discard values for very large distances
                {
                  //halfway between the two normal vectors from the two grains
                  Real l = m1.size();
                  m1 /= l;

                  //same representation in the other grain
                  RealVectorValue m2 = RA*m1;

                  //compute the inclination angle for the common rotation axis
                  //      Moose::out << m1 << std::endl;
                  Real phi = std::acos(std::abs(m1*ax));
                  Real theta1;
                  Real theta2;

                  //partition the total rotation angle "theta"
                  if (std::abs(ax*m1) > 0.9999)   //check if best-fitting GB is pure twist
                  {
                    theta1 = -theta/2;   //eta is meaningless for a twist boundary
                    theta2 = theta/2;
                  }
                  else
                  {
                    //Project m1 and m2 into the plane normal to ax and determine
                    //the rotation angles of them relative to dir
                    theta1 = std::atan2(dir2*m1,dir*m1);
                    theta2 = std::atan2(dir2*m2,dir*m2);
                  }

                  //Reduce both angles to interval (-period/2,period/2)
                  //semi-open with a small numerical error
                  theta2 = theta2-std::round(theta2/period)*period;
                  theta1 = theta1-std::round(theta1/period)*period;

                  //implement the semi-open interval in order to avoid an annoying
                  //numerical problem where certain representations are Real-counted
                  if (std::abs(theta2+period/2) < 0.000001)
                  {
                    theta2 = theta2+period;
                  }
                  if (std::abs(theta1+period/2) < 0.000001)
                  {
                    theta1 = theta1+period;
                  }

                  /* Since this is only being run on fcc elements, which are
                  centrosymmetric, and all dir vectors are 2-fold axes, then
                  the operations of swapping theta1 and theta2, and of
                  multilying both by -1, are symmetries for the energy
                  function. This lets us fold everything into a small right
                  triangle in (ksi,eta) space: */

                  Real ksi = std::abs(theta2-theta1);
                  Real eta = std::abs(theta2+theta1);


                  //round everything to 1e-6, so that negligible numerical differences are dropped
                  distances[thisindex] = 0.000001*std::round(dis*1000000);
                  ksis[thisindex] = 0.000001*std::round(ksi*1000000);
                  etas[thisindex] = 0.000001*std::round(eta*1000000);
                  phis[thisindex] = 0.000001*std::round(phi*1000000);

                  thisindex = thisindex+1;
                }
                else
                {
                  //discard large distances
                  distances.erase(distances.begin()+thisindex);
                  ksis.erase(ksis.begin()+thisindex);
                  etas.erase(etas.begin()+thisindex);
                  phis.erase(phis.begin()+thisindex);
                }
              }
              else
              {
                //dump excess preallocated slots
                distances.erase(distances.begin()+thisindex);
                ksis.erase(ksis.begin()+thisindex);
                etas.erase(etas.begin()+thisindex);
                phis.erase(phis.begin()+thisindex);
              }
            }
          }
          //sort values by distance, ksi, eta, and phi.
          unsigned int change = 1;
          while (change > 0)
          {
            change = 0;
            for (unsigned int i = 0; i < (distances.size()-1); i++)
            {

              if (distances[i+1] < distances[i])
              {
                Real temp = distances[i];
                distances[i] = distances[i+1];
                distances[i+1] = temp;
                temp = ksis[i];
                ksis[i] = ksis[i+1];
                ksis[i+1] = temp;
                temp = etas[i];
                etas[i] = etas[i+1];
                etas[i+1] = temp;
                temp = phis[i];
                phis[i] = phis[i+1];
                phis[i+1] = temp;
                change++;
              }
              else if (distances[i+1] == distances[i])
              {
                if (ksis[i+1] < ksis[i])
                {
                  Real temp = distances[i];
                  distances[i] = distances[i+1];
                  distances[i+1] = temp;
                  temp = ksis[i];
                  ksis[i] = ksis[i+1];
                  ksis[i+1] = temp;
                  temp = etas[i];
                  etas[i] = etas[i+1];
                  etas[i+1] = temp;
                  temp = phis[i];
                  phis[i] = phis[i+1];
                  phis[i+1] = temp;
                  change++;
                }
                else if (ksis[i+1] == ksis[i])
                {
                  if (etas[i+1] < etas[i])
                  {
                    Real temp = distances[i];
                    distances[i] = distances[i+1];
                    distances[i+1] = temp;
                    temp = ksis[i];
                    ksis[i] = ksis[i+1];
                    ksis[i+1] = temp;
                    temp = etas[i];
                    etas[i] = etas[i+1];
                    etas[i+1] = temp;
                    temp = phis[i];
                    phis[i] = phis[i+1];
                    phis[i+1] = temp;
                    change++;
                  }
                  else if (etas[i+1] == etas[i])
                  {
                    if (phis[i+1] < phis[i])
                    {
                      Real temp = distances[i];
                      distances[i] = distances[i+1];
                      distances[i+1] = temp;
                      temp = ksis[i];
                      ksis[i] = ksis[i+1];
                      ksis[i+1] = temp;
                      temp = etas[i];
                      etas[i] = etas[i+1];
                      etas[i+1] = temp;
                      temp = phis[i];
                      phis[i] = phis[i+1];
                      phis[i+1] = temp;
                      change++;
                    }
                  }
                }
              }
            }
          }


          //remove duplicate values. Real-counting in the same representation of one boundary messes up the
          //weighing functions in weightedMeanEnergy()
          for (unsigned int i = 0; i < (distances.size()-1); i++)
          {
            for (unsigned int j = i+1; j < distances.size(); j++)
            {
              if (distances[i] == distances[j] && ksis[i] == ksis[j] && etas[i] == etas[j] && phis[i] == phis[j])
              {
                distances.erase(distances.begin()+j);
                ksis.erase(ksis.begin()+j);
                etas.erase(etas.begin()+j);
                phis.erase(phis.begin()+j);
                j--;
              }
            }
          }

          //remove excess zeroes
          for (unsigned int i = 0; i < distances.size(); i++)
          {
            if (distances[i]== 0 && ksis[i] == 0 && etas[i] == 0 && phis[i] == 0)
            {
              distances.erase(distances.begin()+i);
              ksis.erase(ksis.begin()+i);
              etas.erase(etas.begin()+i);
              phis.erase(phis.begin()+i);
            }
          }


          //store in 4xn geom vector
          geom.resize(4);
          for (unsigned int i = 0; i < distances.size(); i++)
          {
            geom[0].push_back(distances[i]);
            geom[1].push_back(ksis[i]);
            geom[2].push_back(etas[i]);
            geom[3].push_back(phis[i]);
          }

        }

        void
        GBBulatovFit::quat2Mat(std::vector<Real> q, RealTensorValue & _M)
        {
          /* converts a quaternion into a rotation matrix */
          Real d = q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]; //making sure it is normalized

          _M(0,0) = (q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3])/d;
          _M(0,1) = (2*q[1]*q[2]-2*q[0]*q[3])/d;
          _M(0,2) = (2*q[1]*q[3]+2*q[0]*q[2])/d;
          _M(1,0) = (2*q[1]*q[2]+2*q[0]*q[3])/d;
          _M(1,1) = (q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3])/d;
          _M(1,2) = (2*q[2]*q[3]-2*q[0]*q[1])/d;
          _M(2,0) = (2*q[1]*q[3]-2*q[0]*q[2])/d;
          _M(2,1) = (2*q[2]*q[3]+2*q[0]*q[1])/d;
          _M(2,2) = (q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3])/d;
        }

        void
        GBBulatovFit::mat2Quat(RealTensorValue m, std::vector<Real> & q)
        {
          /* converts a rotation matrix, assumed orthonormal, into a unit quaternion */
          /* mat2Quat(m,q)
          Auxiliary function converts a rotation matrix, assumed orthonormal, into
          a unit quaternion. */

          Real t = m(0,0)+m(1,1)+m(2,2);
          Real q0 = std::sqrt(1+t)/2;
          if (t > -0.999999999)
          {
            q[1] = (m(1,2)-m(2,1))/(4*q0);
            q[2] = (m(2,0)-m(0,2))/(4*q0);
            q[3] = (m(0,1)-m(1,0))/(4*q0);
          }
          else
          {
            q0 = 0;
            Real q3 = std::sqrt(-(m(0,0)+m(1,1))/2);
            if (std::abs(q3) > 0.00000002)   //Check for singularity, allowing numerical error
            {
              q[1] = m(0,2)/(2*q3);
              q[2] = m(1,2)/(2*q3);
              q[3] = q3;
            }
            else
            {
              Real q1 = std::sqrt((m(0,0)+1)/2);
              if (q1 != 0)
              {
                q[1] = q1;
                q[2] = m(1,0)/(2*q1);
                q[3] = 0;
              }
              else
              {
                q[1] = 0;
                q[2] = 1;
                q[3] = 0;
              }
            }
          }
          q[0] = q0;
          q[1] = -q[1];
          q[2] = -q[2];
          q[3] = -q[3];
        }

        Real
        GBBulatovFit::weightedMeanEnergy(std::vector<std::vector<Real> > geom100,
          std::vector<std::vector<Real> > geom110, std::vector<std::vector<Real> > geom111)
          {
            /* _en = weightedMeanEnergy(geom100,geom110,geom111)
            Calculate the energy for a single grain boundary.

            Input variables geom100, geom110, and geom111 are each matrices with 4
            columns, giving the non-redundant representations of the boundary about each
            set of axes as generated by distancesToSet().  See the comments in that
            function for further information.  The columns are distance;ksi;eta;phi.

            Input variable _par_vec is a 43-element vector as created by makeParVec().
            This specifies all of the parameters needed for the 5DOF energy function
            on a specific fcc metal.

            Return variable _en is the energy in J/m^2.

            The physical relevance of the parameters is commented wherever they
            appear, in this function and in the set functions.
            */

            //Pull out the parameters relevant to the weighting of the 100, 110, and 111 sets
            Real eRGB = _par_vec[0];   //The only dimensioned parameter. The energy scale,
            //set by the energy of a random boundary.
            Real d0100 = _par_vec[1];  //Maximum distance for the 100 set. Also the distance
            //scale for the RSW weighting function.
            Real d0110 = _par_vec[2];  //Same for the 110 set
            Real d0111 = _par_vec[3];  //Same for the 111 set
            Real weight100 = _par_vec[4];  //Weight for the 100 set, relative to the random boundary
            Real weight110 = _par_vec[5];  //Same for 110
            Real weight111 = _par_vec[6];  //Same for 111

            Real offset = 0.00001;  //Cutoff of weighting function at small d for numerical purposes

            //Calculate energy lists in units of eRGB
            std::vector<Real> e100;
            std::vector<Real> e110;
            std::vector<Real> e111;
            set100(geom100,e100);
            set110(geom110,e110);
            set111(geom111,e111);

            std::vector<Real> s100;
            std::vector<Real> s110;
            std::vector<Real> s111;

            //Calculate the weights, in a manner designed to give an RSW-like function of distance
            //Note it calculates a weight for every representation of the boundary in each set
            for (unsigned int i = 0; i < geom100[0].size(); i++)
            {
              if (geom100[0][i] > d0100)   //weight saturates at zero above d0
              {
                s100.push_back(1);
              }
              else if (geom100[0][i] < (offset*d0100))   //Avoid nan's, replace with something small but finite
              {
                s100.push_back(offset*libMesh::pi/2);
              }
              else
              {
                s100.push_back(std::sin((libMesh::pi/2)*geom100[0][i]/d0100));
              }
            }
            //same for 110
            for (unsigned int i = 0; i < geom110[0].size(); i++)
            {
              if (geom110[0][i] > d0110)
              {
                s110.push_back(1);
              }
              else if (geom110[0][i] < (offset*d0110))
              {
                s110.push_back(offset*libMesh::pi/2);
              }
              else
              {
                s110.push_back(std::sin((libMesh::pi/2)*geom110[0][i]/d0110));
              }
            }
            //same for 111
            for (unsigned int i = 0; i < geom111[0].size(); i++)
            {
              if (geom111[0][i] > d0111)
              {
                s111.push_back(1);
              }
              else if (geom111[0][i] < (offset*d0111))
              {
                s111.push_back(offset*libMesh::pi/2);
              }
              else
              {
                s111.push_back(std::sin((libMesh::pi/2)*geom111[0][i]/d0111));
              }
            }

            Real w100;
            Real w110;
            Real w111;
            Real sum100 = 0;
            Real sum110 = 0;
            Real sum111 = 0;
            Real sumw100 = 0;
            Real sumw110 = 0;
            Real sumw111 = 0;
            //calculate weights and sum up results
            for (unsigned int i = 0; i < geom100[0].size(); i++)
            {
              w100 = (1/(s100[i]*(1-0.5*std::log(s100[i])))-1)*weight100;
              sum100 += (e100[i]*w100);
              sumw100 += w100;
            }
            for (unsigned int i = 0; i < geom110[0].size(); i++)
            {
              w110 = (1/(s110[i]*(1-0.5*std::log(s110[i])))-1)*weight110;
              sum110 += (e110[i]*w110);
              sumw110 += w110;
            }
            for (unsigned int i = 0; i < geom111[0].size(); i++)
            {
              w111 = (1/(s111[i]*(1-0.5*std::log(s111[i])))-1)*weight111;
              sum111 += (e111[i]*w111);
              sumw111 += w111;
            }

            Real _en = eRGB*(sum100+sum110+sum111)/(sumw100+sumw110+sumw111);  //calculate energy

            return _en;
          }

          void
          GBBulatovFit::set100(std::vector<std::vector<Real> > geom100, std::vector<Real> & e100)
          {

            /* set100(geom100,e100)
            Calculate the dimensionless contribution to the boundary based on the
            nearby <100> rotations.  Meant to be called by weightedMeanEnergy(), but
            also can be a stand-alone function for purposes of plotting cross
            sections through the function.
            Input variables geom100 and _par_vec are as generated by distancesToSet()
            and makeParVec().  See comments in those functions for more information.
            */

            Real pwr1 = _par_vec[7];  //100 tilt/twist mix power law: twist
            Real pwr2 = _par_vec[8];  //100 tilt/twist mix power law: tilt

            std::vector<Real> entwist;
            std::vector<Real> entilt;
            twist100(geom100[1], entwist);
            aTGB100(geom100[2],geom100[1],entilt);

            for (unsigned int i = 0; i < geom100[1].size(); i++)
            {
              e100.push_back(entwist[i]*std::pow((1-(2*geom100[3][i]/libMesh::pi)),pwr1)+entilt[i]*std::pow((2*geom100[3][i]/libMesh::pi),pwr2));
              //  Moose::out << geom100[1][i] << ' ' << geom100[2][i] << ' ' << geom100[3][i] <<std::endl;
            }
          }

          void
          GBBulatovFit::twist100(std::vector<Real> ksi, std::vector<Real> & entwist)
          {

            /*
            twist100(ksi, entwist)
            Dimensionless 100 twist contribution to the energy
            */

            Real a = _par_vec[9];   //100 twist maximum energy
            Real b = _par_vec[9]*_par_vec[10];  //100 twist RSW shape factor. The unusual split into
            //two parameters is a holdover from an older version
            Real perio = libMesh::pi/2;  //twist period

            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              ksi[i] = fmod(std::abs(ksi[i]),perio);  //rotation symmetry
              if (ksi[i] > perio/2)
              {
                ksi[i] = perio-ksi[i];
              }
              //implement an RSW function of ksi
              Real xlogx = std::sin(2*ksi[i])*std::log(std::sin(2*ksi[i]));
              if (std::isnan(xlogx))  //force the limit to zero as x->0
              {
                xlogx = 0;
              }
              entwist.push_back(a*std::sin(2*ksi[i])-b*xlogx);
            }

          }

          void
          GBBulatovFit::aTGB100(std::vector<Real> eta, std::vector<Real> ksi, std::vector<Real> & entilt)
          {

            /*
            aTGB100(eta,ksi,_par_vec,entilt)
            This function is a fit to the energies of all 100-tilt boundaries
            */
            Real pwr = _par_vec[11];  //100 aTGB interpolation power law
            Real period = libMesh::pi/2;

            std::vector<Real> en1;
            std::vector<Real> en2;
            std::vector<Real> temp;

            sTGB100(ksi,en1);   //value at eta = 0
            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              temp.push_back(period-ksi[i]);
            }
            sTGB100(temp,en2);   //value at eta = pi/2

            //eta dependence is a power law that goes from the higher to the lower,
            //whichever direction that is
            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              if (en1[i] >= en2[i])
              {
                entilt.push_back(en1[i]-(en1[i]-en2[i])*std::pow((eta[i]/period),pwr));
              }
              else
              {
                entilt.push_back(en2[i]-(en2[i]-en1[i])*std::pow((1-(eta[i]/period)),pwr));
              }
            }

          }

          void
          GBBulatovFit::sTGB100(std::vector<Real> ksi, std::vector<Real> & _en)
          {

            /*
            sTGB100(ksi,_par_vec,_en)
            dimensionless 100 symmetric tilt energy
            This is implemented as a piecewise-RSW function, specified by energy parameters
            _en, angle breaks th, and shape factors a.
            */

            Real en2 = _par_vec[12];  //peak before first Sigma5
            Real en3 = _par_vec[13];  //first Sigma5
            Real en4 = _par_vec[14];  //peak between Sigma5's
            Real en5 = _par_vec[15];  //second Sigma5
            Real en6 = _par_vec[16];  //Sigma17
            Real th2 = _par_vec[17];  //position of peak before first Sigma5
            Real th4 = _par_vec[18];  //position of peak between Sigma5's
            Real th6 = 2*std::acos(5.0/std::sqrt(34.0));  //Sigma17 rotation angle
            Real a12 = 0.5;  //RSW shape factor. In previous versions, these were allowed to
            Real a23 = a12;  //vary, however there were too few vicinal boundaries in the
            Real a34 = a12;  //ensemble to constrain them. We found that forcing the great
            Real a45 = a12;  //majority of them to be 0.5 helped to constrain the fit and
            Real a56 = a12;  //produced reasonable results. This holds true for most of the
            Real a67 = a12;  //RSW shape factors throughout this code.
            Real en1 = 0;  //Sigma1 at left end
            Real en7 = 0;  //Sigma1 at right end
            Real th1 = 0;  //Sigma1 at left end
            Real th3 = std::acos(0.8);  //first Sigma5
            Real th5 = std::acos(0.6);  //second Sigma5
            Real th7 = libMesh::pi/2;  //Sigma1 at right end

            //piecewise RSW function
            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              if (ksi[i] <= th2)
              {
                _en.push_back(en1+(en2-en1)*rSW(ksi[i],th1,th2,a12));
              }
              else if (ksi[i] > th2 && ksi[i] <= th3)
              {
                _en.push_back(en3+(en2-en3)*rSW(ksi[i],th3,th2,a23));
              }
              else if (ksi[i] > th3 && ksi[i] <= th4)
              {
                _en.push_back(en3+(en4-en3)*rSW(ksi[i],th3,th4,a34));
              }
              else if (ksi[i] > th4 && ksi[i] <= th5)
              {
                _en.push_back(en5+(en4-en5)*rSW(ksi[i],th5,th4,a45));
              }
              else if (ksi[i] > th5 && ksi[i] <= th6)
              {
                _en.push_back(en6+(en5-en6)*rSW(ksi[i],th6,th5,a56));
              }
              else if (ksi[i] > th6 && ksi[i] <= th7)
              {
                _en.push_back(en7+(en6-en7)*rSW(ksi[i],th7,th6,a67));
              }
            }
          }

          void
          GBBulatovFit::set110(std::vector<std::vector<Real> > geom110, std::vector<Real> & e110)
          {

            /* set110(geom110, e100)
            Dimensionless contribution to energy from <110> rotations
            Very similar to set100; see comments therein for general information.
            Comments in this file will be limited to 110-specific information.
            */

            Real pwr1 = _par_vec[19];  //110 tilt/twist mix power law: twist
            Real pwr2 = _par_vec[20];  //110 tilt/twist mix power law: tilt

            std::vector<Real> entwist;
            std::vector<Real> entilt;
            twist110(geom110[1],entwist);
            aTGB110(geom110[2],geom110[1],entilt);

            for (unsigned int i = 0;i < geom110[1].size(); i++)
            {
              e110.push_back(entwist[i]*std::pow((1-(2*geom110[3][i]/libMesh::pi)),pwr1)+entilt[i]*std::pow((2*geom110[3][i]/libMesh::pi),pwr2));
            }
          }

          void
          GBBulatovFit::twist110(std::vector<Real> ksi, std::vector<Real> & entwist)
          {

            /*
            twist110(ksi,_par_vec,entwist)
            See comments on set110
            */

            Real th1 = _par_vec[21];  //110 twist peak position
            Real en1 = _par_vec[22];  //110 twist energy peak value
            Real en2 = _par_vec[23];  //Sigma3 energy (110 twist, so not a coherent twin)
            Real en3 = _par_vec[24];  //energy at the symmetry point
            Real a01 = 0.5;
            Real a12 = 0.5;
            Real a23 = 0.5;
            Real th2 = std::acos(1.0/3.0);  //Sigma3
            Real th3 = libMesh::pi/2;  //110 90-degree boundary is semi-special, although not a CSL
            Real perio = libMesh::pi;  //the twist period


            for (unsigned int i = 0;i < ksi.size(); i++)
            {
              ksi[i] = std::fmod(std::abs(ksi[i]),perio);  //rotation symmetry

              if (ksi[i] > perio/2)
              {
                ksi[i] = perio-ksi[i];
              }
              if (ksi[i] <= th1)
              {
                entwist.push_back(en1*rSW(ksi[i],0,th1,a01));
              }
              else if (ksi[i] > th1 && ksi[i] <= th2)
              {
                entwist.push_back(en2+(en1-en2)*rSW(ksi[i],th2,th1,a12));
              }
              else if (ksi[i] > th2)
              {
                entwist.push_back(en3+(en2-en3)*rSW(ksi[i],th3,th2,a23));
              }
            }

          }

          void
          GBBulatovFit::aTGB110(std::vector<Real> eta, std::vector<Real> ksi, std::vector<Real> & entilt)
          {

            /*
            aTGBs110(eta,ksi, entilt)
            See comments on set110
            */

            Real a  =  _par_vec[25];  //110 aTGB interpolation RSW shape factor
            Real period = libMesh::pi;

            std::vector<Real> en1;
            std::vector<Real> en2;
            std::vector<Real> temp;

            sTGB110(ksi,en1);
            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              temp.push_back(period-ksi[i]);
            }
            sTGB110(temp,en2);

            //Power-law interpolation did not work well in this case. Did an RSW function instead
            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              if (en1[i] >= en2[i])
              {
                entilt.push_back(en2[i]+(en1[i]-en2[i])*rSW(eta[i],libMesh::pi,0,a));
              }
              else
              {
                entilt.push_back(en1[i]+(en2[i]-en1[i])*rSW(eta[i],0,libMesh::pi,a));
              }
            }

          }

          void
          GBBulatovFit::sTGB110(std::vector<Real> ksi, std::vector<Real> & _en)
          {

            /*
            sTGB110(ksi, _en)
            See comments on set110
            */

            Real en2 = _par_vec[26];  //peak between Sigma1 and Sigma3
            Real en3 = _par_vec[27];  //Coherent Sigma3 twin relative energy; one of the more
            //important element-dependent parameters
            Real en4 = _par_vec[28];  //energy peak between Sigma3 and Sigma11
            Real en5 = _par_vec[29];  //Sigma11 energy
            Real en6 = _par_vec[30];  //energy peak between Sigma11 and Sigma1
            Real th2 = _par_vec[31];  //peak between Sigma1 and Sigma3
            Real th4 = _par_vec[32];  //peak between Sigma3 and Sigma11
            Real th6 = _par_vec[33];  //peak between Sigma11 and higher Sigma1
            Real a12 = 0.5;
            Real a23 = a12;
            Real a34 = a12;
            Real a45 = a12;
            Real a56 = a12;
            Real a67 = a12;
            Real en1 = 0;
            Real en7 = 0;
            Real th1 = 0;
            Real th3 = std::acos(1.0/3.0);  //Sigma3
            Real th5 = std::acos(-7.0/11.0);  //Sigma11
            Real th7 = libMesh::pi;

            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              ksi[i] = libMesh::pi-ksi[i];  //This is a legacy of an earlier (ksi,eta) mapping
              if (ksi[i] <= th2)
              {
                _en.push_back(en1+(en2-en1)*rSW(ksi[i],th1,th2,a12));
              }
              else if (ksi[i] > th2 && ksi[i] <= th3)
              {
                _en.push_back(en3+(en2-en3)*rSW(ksi[i],th3,th2,a23));
              }
              else if (ksi[i] > th3 && ksi[i] <= th4)
              {
                _en.push_back(en3+(en4-en3)*rSW(ksi[i],th3,th4,a34));
              }
              else if (ksi[i] > th4 && ksi[i] <= th5)
              {
                _en.push_back(en5+(en4-en5)*rSW(ksi[i],th5,th4,a45));
              }
              else if (ksi[i] > th5 && ksi[i] <= th6)
              {
                _en.push_back(en5+(en6-en5)*rSW(ksi[i],th5,th6,a56));
              }
              else if (ksi[i] > th6 && ksi[i] <= th7)
              {
                _en.push_back(en7+(en6-en7)*rSW(ksi[i],th7,th6,a67));
              }
            }

          }

          void
          GBBulatovFit::set111(std::vector<std::vector<Real> > geom111, std::vector<Real> & e111)
          {

            /* set111(geom111,e111)
            Dimensionless contribution to energy from <111> rotations
            Very similar to set100; see comments therein for general information.
            Comments in this file will be limited to 111-specific information.
            */

            Real a = _par_vec[34];   //linear part of 111 tilt/twist interpolation
            Real b = a-1;   //ensures correct value at x = 1

            std::vector<Real> entwist;
            std::vector<Real> entilt;
            twist111(geom111[1],entwist);
            aTGB111(geom111[2],geom111[1],entilt);

            //This one fit well enough with a simple one-parameter parabola that the more
            //complicated power laws in the other sets weren't needed
            for (unsigned int i = 0;i < geom111[1].size(); i++)
            {
              e111.push_back(entwist[i]+(entilt[i]-entwist[i])*(a*2*geom111[3][i]/libMesh::pi-b*std::pow((2*geom111[3][i]/libMesh::pi),2)));
            }

          }

          void
          GBBulatovFit::twist111(std::vector<Real> ksi, std::vector<Real> & entwist)
          {

            /*
            twist111(ksi, entwist)
            See comments on set111
            */

            Real thd = _par_vec[36];  //111 twist peak position
            Real enm = _par_vec[37];  //111 twist energy at the peak
            Real en2 = _par_vec[27];  //Coherent Sigma3 twin shows up in two distinct places in the code
            Real a1 = _par_vec[35];  //111 twist RSW shape parameter
            Real a2 = a1;

            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              if (ksi[i] > libMesh::pi/3)
              {
                ksi[i] = 2 * libMesh::pi / 3 - ksi[i];
              }
              if (ksi[i] <= thd)
              {
                entwist.push_back(enm*rSW(ksi[i],0,thd,a1));
              }
              else
              {
                entwist.push_back(en2+(enm-en2)*rSW(ksi[i],libMesh::pi/3,thd,a2));
              }
            }


          }

          void
          GBBulatovFit::aTGB111(std::vector<Real> eta, std::vector<Real> ksi, std::vector<Real> & entilt)
          {

            /* aTGBs111(eta,ksi,entilt)
            This function is a fit to the energies of all 111-tilt boundaries
            */

            /* There's an additional symmetry in 111 aTGBs that doesn't exist in 100 or
            110 aTGBs.  This is because a rotation about [111] equal to half the period
            (i.e. 60 degrees) is equivalent to a mirror reflection in the (111)
            plane.  Both are Sigma3 operations.  The same is not true of the
            45-degree [100] or the 90-degree [110] rotation.
            The following loop accounts for this extra symmetry.
            */
            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              if (ksi[i] > libMesh::pi/3)
              {
                ksi[i] = 2*libMesh::pi/3-ksi[i];
              }
              if (eta[i] > libMesh::pi/3)
              {
                eta[i] = 2*libMesh::pi/3-eta[i];
              }
            }

            /* Below the following value of ksi, we ignore the eta dependence.  This is
            because there's little evidence that it actually varies.  Above this
            value, we interpolate on an RSW function that follows the Sigma3 line,
            which is also a line of symmetry for the function.
            */
            Real ksim = _par_vec[38];  //111 aTGB ksi break
            Real enmax = _par_vec[39];  //Energy at the peak (ksi == ksim)
            Real enmin = _par_vec[40];  //energy at the minimum (Sigma3, eta == 0)
            Real encnt = _par_vec[41];  //energy at the symmetry point (Sigma3, eta == pi/3)
            Real a1 = 0.5;
            Real a2 = 0.5;
            /* eta scaling parameter for 111 aTGB RSW function on Sigma3 line
            This RSW function is unusual in that the change in shape of the
            function is much better captured by changing the angular scale rather
            than changing the dimensionless shape factor.
            */
            Real etascale = _par_vec[42];
            Real chi;

            for (unsigned int i = 0; i < ksi.size(); i++)
            {
              if (ksi[i] <= ksim)
              {
                entilt.push_back(enmax*rSW(ksi[i],0,ksim,a1));
              }
              else
              {
                //chi is the shape of the function along the Sigma3 line.
                chi = enmin+(encnt-enmin)*rSW(eta[i],0,libMesh::pi/(2*etascale),0.5);
                entilt.push_back(chi+(enmax-chi)*rSW(ksi[i],libMesh::pi/3,ksim,a2));
              }
            }


          }


          Real
          GBBulatovFit::rSW(Real theta, Real thetaMin, Real thetaMax, Real a)
          {

            /* _en = rSW(theta,theta1,theta2,a)
            This function computes the value of Read-Shockley-Wolf function at theta.
            The RSW function is normalized to be 1.0 at theta2 and 0.0 at theta1.

            theta             angle at which to compute the function
            thetaMin          the starting angle of the interval
            thetaMax          the end angle of the interval
            a                 parameter defining the shape of the RSW function
            */

            Real dtheta = thetaMax-thetaMin;  //Interval of angles where defined
            theta = (theta-thetaMin)*libMesh::pi/(dtheta*2);  //Normalized angle

            //RSW function evaluation
            Real _en;
            if (std::sin(theta) >= 0.000001)  //Cut off a small sins to avoid 0*infinity problem.
            //The proper limit is 0.
            {
              _en = std::sin(theta)-a*(std::sin(theta)*std::log(std::sin(theta)));
            }
            else
            {
              _en = std::sin(theta);
            }

            return _en;
          }
