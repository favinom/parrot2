[Problem]
solve = false
[]

[Mesh]
[file]
  type = FileMeshGenerator
  file = AdvectionOut_5120_4.e
[]
[]

[Variables]
  active = 'u'
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./pressure_c]
    order = FIRST
    family = LAGRANGE
  [../]
[]


[Functions]
  [./fine_function]
    type = SolutionFunction
    solution = fine_solution
  [../]
[]

[UserObjects]
  [./fine_solution]
    type = SolutionUserObject
    mesh = AdvectionOut_5120_4.e
    es = AdvectionOut_5120_4.e
    system_variables = pressure
  [../]
[]


[Postprocessors]

  [./h1_error]
    type = ElementH1Error
    variable = pressure_c
    function = fine_function
    execute_on = 'timestep_end'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

[]

#[Outputs]
#file_base = output_sub
#exodus = true
#[]

