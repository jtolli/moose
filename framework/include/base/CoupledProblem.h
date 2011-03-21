#ifndef COUPLEDPROBLEM_H
#define COUPLEDPROBLEM_H

#include <vector>
#include <string>
#include <map>

#include "Moose.h"
#include "Problem.h"
#include "ParallelUniqueId.h"
#include "InputParameters.h"
#include "Function.h"
#include "GeometricSearchData.h"
// libMesh
#include "equation_systems.h"
#include "numeric_vector.h"
#include "sparse_matrix.h"
#include "transient_system.h"
#include "nonlinear_implicit_system.h"

class SubProblem;
class MooseVariable;
class MooseMesh;

class CoupledProblem : public Problem
{
public:
  CoupledProblem(MooseMesh * mesh);
  virtual ~CoupledProblem();

  virtual EquationSystems & es() { return _eq; }
  MooseMesh & mesh() { return *_mesh; }

  virtual Problem * parent() { return NULL; }

  void addSubProblem(const std::string & file_name, SubProblem *subproblem);
  SubProblem *subProblem(const std::string & name);

  void solveOrder(const std::vector<std::string> & solve_order);

  // API /////
  virtual bool hasVariable(const std::string & var_name);
  virtual MooseVariable & getVariable(THREAD_ID tid, const std::string & var_name);

  virtual AssemblyData & assembly(THREAD_ID tid);
  virtual void prepare(const Elem * elem, THREAD_ID tid);
  virtual void reinitElem(const Elem * elem, THREAD_ID tid);
  virtual void reinitElemFace(const Elem * elem, unsigned int side, unsigned int bnd_id, THREAD_ID tid);
  virtual void reinitNode(const Node * node, THREAD_ID tid);
  virtual void reinitNodeFace(const Node * node, unsigned int bnd_id, THREAD_ID tid);

  virtual void subdomainSetup(unsigned int subdomain, THREAD_ID tid);

  // Solve /////
  virtual void init();
  virtual void update();

  virtual void computeResidual(NonlinearImplicitSystem & sys, const NumericVector<Number> & soln, NumericVector<Number> & residual);
  virtual void computeJacobian(NonlinearImplicitSystem & sys, const NumericVector<Number> & soln, SparseMatrix<Number> &  jacobian);

  // Initial conditions /////
  virtual Number initialValue (const Point & p, const Parameters & parameters, const std::string & sys_name, const std::string & var_name);
  virtual Gradient initialGradient (const Point & p, const Parameters & parameters, const std::string & sys_name, const std::string & var_name);

  virtual void initialCondition(EquationSystems & es, const std::string & system_name);

  // Postprocessors /////
  virtual void computePostprocessors(int pps_type = Moose::PPS_TIMESTEP);
  virtual void outputPostprocessors();
  virtual Real & getPostprocessorValue(const std::string & name, THREAD_ID tid = 0);

  // Materials /////
  virtual void reinitMaterials(unsigned int blk_id, THREAD_ID tid);
  virtual void reinitMaterialsFace(unsigned int blk_id, unsigned int side, THREAD_ID tid);

  // Transient /////
  virtual void transient(bool trans);
  virtual bool transient() { return _transient; }

  virtual Real & time() { return _time; }
  virtual int & timeStep() { return _t_step; }
  virtual Real & dt() { return _dt; }
  virtual Real & dtOld() { return _dt_old; }

  virtual void copySolutionsBackwards();

  virtual void dump();

  virtual Output & out() { return _out; }
  virtual void output();

  virtual GeometricSearchData & geomSearchData() { return *_geometric_search_data; }

protected:
  std::map<std::string, SubProblem *> _subproblems;
  std::vector<std::string> _solve_order;

  /// Keep track of the correspondence between libMesh objects and Moose objects
  std::map<std::string, SubProblem *> _map;

  MooseMesh * _mesh;
  EquationSystems _eq;

  bool _transient;
  Real & _time;
  int & _t_step;
  Real & _dt;
  Real _dt_old;

  // Output system
  Output _out;

  // FIXME: this should not be here, right now it is here for satisfying the interface (which is not very clean - because of this). understand?
  GeometricSearchData * _geometric_search_data;
};

#endif /* COUPLEDPROBLEM_H_ */
