[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 12 # Number of elements in the x-direction
  ny = 12 # Number of elements in the y-direction
  nz = 0 # Number of elements in the z-direction
  xmin = 0 # minimum x-coordinate of the mesh
  xmax = 1000 # maximum x-coordinate of the mesh
  ymin = 0 # minimum y-coordinate of the mesh
  ymax = 1000 # maximum y-coordinate of the mesh
  zmin = 0
  zmax = 0
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 3 # Initial uniform refinement of the mesh
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 4 # Number of grains
  var_name_base = gr # Base name of grains
  grain_num = 4
  length_scale = 1.0e-9
  time_scale = 1.0e-9
  wGB = 14
  Anisotropic_GB_file_name = anisotropy_mobility.txt   # anisotropy_energy.txt
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
    # Custom action that created all of the grain variables
    order = FIRST # element type used by each grain variable
    family = LAGRANGE
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalHexGrainIC]

    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = FIRST
    family = LAGRANGE
  [../]
  [./var_indices]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
  [./anisotropy0]
    type = ACGrGrAnis
    variable = gr0
    v = 'gr1 gr2 gr3'
    op = 0
  [../]
  [./anisotropy1]
    type = ACGrGrAnis
    variable = gr1
    v = 'gr0 gr2 gr3'
    op = 1
  [../]
  [./anisotropy2]
    type = ACGrGrAnis
    variable = gr2
    v = 'gr0 gr1 gr3'
    op = 2
  [../]
  [./anisotropy3]
    type = ACGrGrAnis
    variable = gr3
    v = 'gr0 gr1 gr2'
    op = 3
  [../]

[]

[AuxKernels]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./Periodic]
    [./top_bottom]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropy
    block = 0
    T = 600 # K

    # molar_volume_value = 7.11e-6 #Units:m^3/mol

    inclination_anisotropy = false # true
    outputs = exodus
  [../]
[]

[Postprocessors]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]

  [./gr1_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  [../]
  [./gr2_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr2
  [../]
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'
  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_abs_tol = 1e-11 # Relative tolerance for nonlienar solves
  nl_rel_tol = 1e-8 # Absolute tolerance for nonlienar solves
  start_time = 0
  end_time = 17810
  dt = 5
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 5 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 #Time step will adapt to maintain this number of nonlinear iterations
  [../]
  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 5 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]


[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  csv = true
[]
