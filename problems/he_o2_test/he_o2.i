[Mesh]
  type = GeneratedMesh
  dim = 1
  xmin = 0
  xmax = 1
  nx = 1
[]

[Variables]
  [neTe]
    family = SCALAR
    order = FIRST
    initial_condition = 1e6
  []
  # ODE variables
  [./e]
    family = SCALAR
    order = FIRST
    initial_condition = 1e6
    # scaling = 1e-2
  [../]

  [./O2+]
    family = SCALAR
    order = FIRST
    initial_condition = 1e6
    # scaling = 1e-2
  [../]

  [./O+]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./O2]
    family = SCALAR
    order = FIRST
    initial_condition = 1.1e17 
    #scaling = 1e-15
  [../]

  [./O]
    family = SCALAR
    order = FIRST
    initial_condition = 1e6
    # scaling = 1e-2
  [../]

  [./O-]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./O21Dg]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./O*]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./O1SIG+]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./O3]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./O3-]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]
  
  [./O2-]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]
  
  [./O1D]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]

  [./He]
    family = SCALAR
    order = FIRST
    initial_condition = 2.2e19
  [../]
  [./He*]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]
  [./He+]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]
  [./He2+]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]
  [./He2*]
    family = SCALAR
    order = FIRST
    initial_condition = 1
  [../]
[]

[ScalarKernels]
  [neTe_power]
    type = InputPower
    variable = neTe 
    value = 330 
  []
  [neTe_dt]
    type = ODETimeDerivative
    variable = neTe
  []
  [./de_dt]
    type = ODETimeDerivative
    variable = e
  [../]

  [./dO2p_dt]
    type = ODETimeDerivative
    variable = O2+
  [../]

  [./dOp_dt]
    type = ODETimeDerivative
    variable = O+
  [../]

  [./dO2_dt]
    type = ODETimeDerivative
    variable = O2
  [../]

  [./dO_dt]
    type = ODETimeDerivative
    variable = O 
  [../]

  [./dOm_dt]
    type = ODETimeDerivative
    variable = O-
  [../]

  [./dO21Dg_dt]
    type = ODETimeDerivative
    variable = O21Dg
  [../]

  [./dO*_dt]
    type = ODETimeDerivative
    variable = O*
  [../]

  [./dO1SIGp_dt]
    type = ODETimeDerivative
    variable = O1SIG+
  [../]

  [./dO3_dt]
    type = ODETimeDerivative
    variable = O3
  [../]

  [./dO3m_dt]
    type = ODETimeDerivative
    variable = O3-
  [../]

  [./dO2m_dt]
    type = ODETimeDerivative
    variable = O2-
  [../]

  [./dO1D_dt]
    type = ODETimeDerivative
    variable = O1D
  [../]

  [./dHe_dt]
    type = ODETimeDerivative
    variable = He
  [../]

  [./dHe*_dt]
    type = ODETimeDerivative
    variable = He*
  [../]

  [./dHep_dt]
    type = ODETimeDerivative
    variable = He+
  [../]

  [./dHe2p_dt]
    type = ODETimeDerivative
    variable = He2+
  [../]

  [./dHe2*_dt]
    type = ODETimeDerivative
    variable = He2* 
  [../]

[]

[GlobalReactions]
  [argon]
    species = 'e O2+ O+ O2 O2- O O- O21Dg O* O1SIG+ O3 O3- O1D He He* He+ He2+ He2*'
    #file_location = 'data'

    # These are parameters required equation-based rate coefficients
    equation_constants = 'Tgas J pi'
    equation_values = '300 2.405 3.141'
    equation_variables = 'Te'
    #sampling_variable = 'reduced_field'
    electron_energy = neTe


    #reactions = 'O2+ + e -> O + O       : 4.8e-7
    #             e + O2 -> O- + O       : {8.8e-11*exp(-4.4/Te)} [-4.4]
    #             e + O2 -> O2+ + e + e  : {9e-10*Te*exp(-12.6/Te)} [-12]
    #             O- + O -> O2 + e       : 2e-10 
    #             O2 + O+ -> O2+ + O     : 2e-11
    #             e + O2 -> O21Dg + e    : {1.7e-9*exp(-3.1/Te)} [-3.1]
    #             O- + O21Dg -> O3 + e   : 3e-10
    #             O- + O21Dg -> O2- + O  : 1e-10
    #             e + O2 -> O + O1D + e  : {5e-8*exp(-8.4/Te)} [-8.4]
    #             e + O -> O1D + e       : {4.2e-9*exp(-2.25/Te)} [-2.25]
    #             O1SIG+ + O3 -> O + O2 + O2  : 1.5e-11
    #             O* + O2 -> O + O1SIG+ : 2e-11
    #             e + O2 -> O + O + e    : {4.2e-9*exp(-5.6/Te)} [-5.6]
    #             O21Dg + O3 -> O2 + O2 + O*   : 1e-11
    #             e + O -> O+ + e + e      : {9e-9*Te^0.7}
    #             O- + O3 -> O3- + O     : 5.3e-10
    #             O3- + O -> O2 + O2 + e : 3e-10
    #             O3- + O -> O2- + O2    : 3.2e-10
    #             He* + O2 -> He + O2+ + e   : 2.4e-10
    #             He + O2 + O -> He + O3     : 6.27e-34
    #             e + He -> He* + e          : {2.3e-10*Te^0.31}
    #             e + He -> He+ + e + e      : {2.5e-12*Te^0.68*exp(-24.6/Te)} [-24.6]
    #             e + He2+ -> He* + He       : {5.386e-7*Te^0.5}
    #             He* + He + He -> He2* + He : 1.3e-33
    #             He+ + He + He -> He2+ + He  : 1.3e-31
    #             He2* + He2* -> He2+ + He + He + e  : 1.5e-9
    #             O1D + O2 -> O + O2             : 3e-11
    #             O1D + He -> O + He             : 1e-13
    #             O1SIG+ + O -> O21Dg + O       : 3e-12
    #             e + O3 -> O- + O2              : {9.3e-10*Te^-0.62}
    #             e + O3 -> O + O2-              : 2e-10
    #             O3- + O2+ -> O + O + O3        : {1.01e-7*(300/Tgas)^0.5}
    #             e + O3 -> O + O2 + e           : {1e-8*(300/Tgas)^0.5}
    #             He+ + O2 -> O+ + O + He        : {1.07e-9*(Tgas/300)^0.5}
    #             O2- + O2 -> O2 + O2 + e        : {2.7e-10*(Tgas/300)^0.5 * exp(-5590/Tgas)}
    #             O2- + He -> O2 + He + e        : {2.7e-10*(Tgas/300)^0.5 * exp(-5590/Tgas)}'
    reactions = 'O2+ + e -> O + O       : 4.8e-7
                 e + O2 -> O- + O       : {8.8e-11*exp(-4.4/Te)} [-4.4]
                 e + O2 -> O2+ + e + e  : {9e-10*Te*exp(-12.6/Te)} 
                 O- + O -> O2 + e       : 2e-10 
                 O2 + O+ -> O2+ + O     : 2e-11
                 e + O2 -> O21Dg + e    : {1.7e-9*exp(-3.1/Te)} 
                 O- + O21Dg -> O3 + e   : 3e-10
                 O- + O21Dg -> O2- + O  : 1e-10
                 e + O2 -> O + O1D + e  : {5e-8*exp(-8.4/Te)} 
                 e + O -> O1D + e       : {4.2e-9*exp(-2.25/Te)} 
                 O1SIG+ + O3 -> O + O2 + O2  : 1.5e-11
                 O* + O2 -> O + O1SIG+ : 2e-11
                 e + O2 -> O + O + e    : {4.2e-9*exp(-5.6/Te)} 
                 O21Dg + O3 -> O2 + O2 + O*   : 1e-11
                 e + O -> O+ + e + e      : {9e-9*Te^0.7}
                 O- + O3 -> O3- + O     : 5.3e-10
                 O3- + O -> O2 + O2 + e : 3e-10
                 O3- + O -> O2- + O2    : 3.2e-10
                 He* + O2 -> He + O2+ + e   : 2.4e-10
                 He + O2 + O -> He + O3     : 6.27e-34
                 e + He -> He* + e          : {2.3e-10*Te^0.31}
                 e + He -> He+ + e + e      : {2.5e-12*Te^0.68*exp(-24.6/Te)} 
                 e + He2+ -> He* + He       : {5.386e-7*Te^0.5}
                 He* + He + He -> He2* + He : 1.3e-33
                 He+ + He + He -> He2+ + He  : 1.3e-31
                 He2* + He2* -> He2+ + He + He + e  : 1.5e-9
                 O1D + O2 -> O + O2             : 3e-11
                 O1D + He -> O + He             : 1e-13
                 O1SIG+ + O -> O21Dg + O       : 3e-12
                 e + O3 -> O- + O2              : {9.3e-10*Te^-0.62}
                 e + O3 -> O + O2-              : 2e-10
                 O3- + O2+ -> O + O + O3        : {1.01e-7*(300/Tgas)^0.5}
                 e + O3 -> O + O2 + e           : {1e-8*(300/Tgas)^0.5}
                 He+ + O2 -> O+ + O + He        : {1.07e-9*(Tgas/300)^0.5}
                 O2- + O2 -> O2 + O2 + e        : {2.7e-10*(Tgas/300)^0.5 * exp(-5590/Tgas)}
                 O2- + He -> O2 + He + e        : {2.7e-10*(Tgas/300)^0.5 * exp(-5590/Tgas)}'
  []
[]

[AuxVariables]
  [./Te]
    order = FIRST
    family = SCALAR
    initial_condition = 4
  [../]
[]

[AuxScalarKernels]
  #[temperature_calculation]
  #  type = ScalarLinearInterpolation
  #  variable = Te
  #  sampler = reduced_field
  #  property_file = 'data/electron_temperature.txt'
  #  execute_on = 'INITIAL TIMESTEP_BEGIN'
  #[]
  [temperature_calculation]
    type = ParsedAuxScalar
    variable = Te
    args = 'e'
    function = '2/(3*e)'
    execute_on = 'TIMESTEP_END'
  []
[]

#[Debug]
#  show_var_residual_norms = true
#[]

[Executioner]
  type = Transient
  end_time = 1e-3
  solve_type = linear
  dtmin = 1e-16
  dtmax = 1e-8
  line_search = none
  steady_state_detection = true
  [./TimeStepper]
    type = IterationAdaptiveDT
    cutback_factor = 0.9
    dt = 1e-12
    growth_factor = 1.01
  [../]
  [TimeIntegrator]
    type = LStableDirk2
  []
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    #ksp_norm = none
  [../]
[]

[Outputs]
  [out]
    type = CSV
  []
  [console]
    type = Console
    execute_scalars_on = 'none'
  []
[]
