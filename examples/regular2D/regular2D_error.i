[Problem]
# solve = false
[]

[Mesh]
[gmg]
dim = 2
type = GeneratedMeshGenerator
nx = ${nxMoose}
ny = ${nxMoose}
xmin= 0.0
xmax= 1.0
ymin= 0.0
ymax= 1.0
[../]

[inclusionsList]
input = gmg 
type = InclusionsMeshModifier
fn = 6
fx_string = '0.5 0.5 0.75 0.75 0.625 0.625'
fy_string = '0.5 0.5 0.75 0.75 0.625 0.625'
fd1_string = '1.0 1.0 0.5 0.5 0.25 0.25'
fd2_string = '1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4 1.0e-4'
fa1_string = '0.0 90.0 0.0 90.0 0.0 90.0'
[../]


[inclusionRefinement]
input = inclusionsList
type  = InclusionRefinement
inclusions_list = inclusionsList
doBoundaryRefinement = true
refinements = '0 0'  # amr uniforme
#outputFileName = test.e
[../]

[]


[Variables]
[pressure] []
[]

[Kernels]
[./diffusion] 
type = PermeabilityDiffusion 
variable = pressure 
[../]
[]

[Materials]
[./flowAndTransport] 
type = FlowAndTransport 
k = 1.0 kFrac = 1.0e4
inclusions_list = inclusionsList 
[../]
[]

# observe that with the second BCs the stiffness matrix is SDP and we can use choleski factorization
[BCs]
[./inflowBC]  type = DirichletBC variable = pressure value = 1.0  boundary = right [../]
[./ouflowBC]  type = NeumannBC   variable = pressure value = 1.0  boundary = left  [../]
[]

[Preconditioning]
[./prec] type = SMP full = true ksp_norm = default [../]
[]

[Executioner]
type = Transient
num_steps = 1
dt = 1
solve_type = LINEAR
petsc_options_iname = '-pc_type -pc_hypre_type '
petsc_options_value = ' hypre    boomeramg     '
[./Quadrature] order = NINTH type = GRID [../]
[]

#[Outputs]
#file_base = output
#exodus = true
#[]

[MultiApps]
[./sub]
type = TransientMultiApp
app_type = Parrot2App
positions = '0 0 0 0.0 0 0'
input_files = elemenet_h1_error.i
execute_on = timestep_end
[../]
[]

[Transfers]
[./tosub]
type = MultiAppMeshFunctionTransfer
direction = to_multiapp
multi_app = sub
source_variable = pressure
variable = pressure_c
execute_on = timestep_end
[../]
[]
