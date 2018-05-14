#> \file
#> \author Chris Bradley
#> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
#>
#> \section LICENSE
#>
#> Version: MPL 1.1/GPL 2.0/LGPL 2.1
#>
#> The contents of this file are subject to the Mozilla Public License
#> Version 1.1 (the "License"); you may not use this file except in
#> compliance with the License. You may obtain a copy of the License at
#> http://www.mozilla.org/MPL/
#>
#> Software distributed under the License is distributed on an "AS IS"
#> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#> License for the specific language governing rights and limitations
#> under the License.
#>
#> The Original Code is OpenCMISS
#>
#> The Initial Developer of the Original Code is University of Auckland,
#> Auckland, New Zealand and University of Oxford, Oxford, United
#> Kingdom. Portions created by the University of Auckland and University
#> of Oxford are Copyright (C) 2007 by the University of Auckland and
#> the University of Oxford. All Rights Reserved.
#>
#> Contributor(s): Elias Ghadam Soltani
#>
#> Alternatively, the contents of this file may be used under the terms of
#> either the GNU General Public License Version 2 or later (the "GPL"), or
#> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
#> in which case the provisions of the GPL or the LGPL are applicable instead
#> of those above. If you wish to allow use of your version of this file only
#> under the terms of either the GPL or the LGPL, and not to allow others to
#> use your version of this file under the terms of the MPL, indicate your
#> decision by deleting the provisions above and replace them with the notice
#> and other provisions required by the GPL or the LGPL. If you do not delete
#> the provisions above, a recipient may use your version of this file under
#> the terms of any one of the MPL, the GPL or the LGPL.
#>

#> \example ClassicalField/Diffusion/DiffusionConstantSource/Python/DiffusionConstantSourceExample.py
## Example program to solve a diffusion equation using OpenCMISS calls.
## \htmlinclude ClassicalField/Diffusion/DiffusionConstantSource/history.html
#<

#!> Main program
#!PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE
# Add Python bindings directory to PATH
import numpy,sys,os,time

#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

#!  USE OpenCMISS
#!  USE OpenCMISS_Iron
# Intialise OpenCMISS-Iron start
from opencmiss.iron import iron
t=time.time()
# Diagnostics
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.IN,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",[""])
#iron.ErrorHandlingModeSet(iron.ErrorHandlingModes.TRAP_ERROR)
#iron.OutputSetOn("Testing")

# Get the computational nodes info
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber    = iron.ComputationalNodeNumberGet()

#!#ifndef NOMPIMOD
#!  USE MPI
#!#endif

#!#ifdef WIN32
#!  USE IFQWIN
#!#endif

#!  IMPLICIT NONE

#!#ifdef NOMPIMOD
#!#include "mpif.h"
#!#endif


#!Elias Variables
#!  REAL(CMISSRP) :: node_value
#!  INTEGER(CMISSIntg) :: node_number, N_row, N_col, i,j, k, adjacentElement
#!  Real(CMISSRP) ::   K_Conductivity, Rho_Density, Cp_SpecificHeat, alpha
#!  INTEGER(CMISSIntg),allocatable :: ElementUserNodes(:)
    
#================================================================================================================================
#  Start Program
#================================================================================================================================
 
#!INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
#!INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
#!INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
#!INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
#!INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
#!INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
#!INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
#!INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
#!INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
#!INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
#!INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
#!INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
#!INTEGER(CMISSIntg), PARAMETER :: ControlLoopNode=0
#!INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=13
#!INTEGER(CMISSIntg), PARAMETER :: AnalyticFieldUserNumber=14

controlLoopNode             = 0
coordinateSystemUserNumber  = 1
regionUserNumber            = 2
basisUserNumber             = 3
generatedMeshUserNumber     = 4
meshUserNumber              = 5
decompositionUserNumber     = 6
geometricFieldUserNumber    = 7
equationsSetFieldUserNumber = 8
dependentFieldUserNumber    = 9
materialsFieldUserNumber    = 10
equationsSetUserNumber      = 11
problemUserNumber           = 12
sourceFieldUserNumber       = 13
analyticFieldUserNumber     = 14  

#================================================================================================================================
#  Initial Data & Default values
#================================================================================================================================


#!  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
#!  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
#!  REAL(CMISSRP), PARAMETER :: LENGTH=4.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: HEIGHT=1.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: LENGTH=4.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: HEIGHT=0.14_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: WIDTH=0.7_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: LENGTH=0.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: WIDTH=2.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: LENGTH=0.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: HEIGHT=2.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: WIDTH=4.0_CMISSRP
#!!  REAL(CMISSRP), PARAMETER :: LENGTH=0.0_CMISSRP

#HEIGHT = 2.0
#WIDTH  = 2.0
#LENGTH = 4.0

# Set the time parameters
timeIncrement   = 0.001
startTime       = 0.0
stopTime  = 0.8001

# Set the output parameters
DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY = 10

# Set the solver parameters
#relativeTolerance   = 1.0E-05  # default: 1.0E-05
#absoluteTolerance   = 1.0E-08  # default: 1.0E-10
#DIVERGENCE_TOLERANCE = 1.0e+10  # default: 1.0e+05
MAXIMUM_ITERATIONS   = 1000   # default: 100000
#RESTART_VALUE        = 3000     # default: 30

#================================================================================================================================
#  Problem Control Panel
#================================================================================================================================



#!  !Program types
  
#!  !Program variables

#!  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS,NUMBER_GLOBAL_Z_ELEMENTS
#!  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS
  
#!  INTEGER(CMISSIntg) :: MPI_IERROR
  
#!  !CMISS variables

#!  TYPE(cmfe_BasisType) :: Basis
#!  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
#!  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
#!  TYPE(cmfe_DecompositionType) :: Decomposition
#!  TYPE(cmfe_EquationsType) :: Equations
#!  TYPE(cmfe_EquationsSetType) :: EquationsSet
#!  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField,SourceField,AnalyticField
#!  TYPE(cmfe_FieldsType) :: Fields
#!  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
#!  TYPE(cmfe_MeshType) :: Mesh
#!  TYPE(cmfe_ProblemType) :: Problem
#!  TYPE(cmfe_ControlLoopType) :: ControlLoop
#!  TYPE(cmfe_RegionType) :: Region,WorldRegion
#!  TYPE(cmfe_SolverType) :: Solver, LinearSolver
#!  TYPE(cmfe_SolverEquationsType) :: SolverEquations

#!  LOGICAL :: EXPORT_FIELD,IMPORT_FIELD
 

#!#ifdef WIN32
#!  !Quickwin type
#!  LOGICAL :: QUICKWIN_STATUS=.FALSE.
#!  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#!#endif
  
#!  !Generic CMISS variables
  
#!  INTEGER(CMISSIntg) :: EquationsSetIndex
#!  INTEGER(CMISSIntg) :: FirstNodeNumber,LastNodeNumber
#! INTEGER(CMISSIntg) :: Err
  
#!#ifdef WIN32
#!  !Initialise QuickWin
#!  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
#!  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
#!  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
#!  !Set the window parameters
#!  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#!  !If attempt fails set with system estimated values
#!  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#!#endif




#!SELECT CASE(1)
#!  CASE(1)
#!    NUMBER_GLOBAL_X_ELEMENTS=4
#!    NUMBER_GLOBAL_Y_ELEMENTS=4
#!    NUMBER_GLOBAL_Z_ELEMENTS=8
#!  CASE(2)
#!    NUMBER_GLOBAL_X_ELEMENTS=2
#!    NUMBER_GLOBAL_Y_ELEMENTS=2
#!    NUMBER_GLOBAL_Z_ELEMENTS=4
#!  CASE(3)
#!    NUMBER_GLOBAL_X_ELEMENTS=4
#!    NUMBER_GLOBAL_Y_ELEMENTS=2
#!    NUMBER_GLOBAL_Z_ELEMENTS=0
#!  CASE(4)
#!    NUMBER_GLOBAL_X_ELEMENTS=8
#!    NUMBER_GLOBAL_Y_ELEMENTS=4
#!    NUMBER_GLOBAL_Z_ELEMENTS=0
#!  CASE(5)
#!    NUMBER_GLOBAL_X_ELEMENTS=8
#!    NUMBER_GLOBAL_Y_ELEMENTS=8
#!    NUMBER_GLOBAL_Z_ELEMENTS=0
#!  CASE(6)    
#!    NUMBER_GLOBAL_X_ELEMENTS=2
#!    NUMBER_GLOBAL_Y_ELEMENTS=2
#!    NUMBER_GLOBAL_Z_ELEMENTS=0
#!  CASE(7)
#!    NUMBER_GLOBAL_X_ELEMENTS=4
#!    NUMBER_GLOBAL_Y_ELEMENTS=8
#!    NUMBER_GLOBAL_Z_ELEMENTS=0
#!  CASE(8)
#!    NUMBER_GLOBAL_X_ELEMENTS=10
#!    NUMBER_GLOBAL_Y_ELEMENTS=10
#!    NUMBER_GLOBAL_Z_ELEMENTS=0
#!END SELECT

#numberGlobalXElements = 4
#numberGlobalYElements = 4
#numberGlobalZElements = 8

#totalNumberOfElements  = numberGlobalXElements*numberGlobalYElements*numberGlobalZElements
meshOrigin = [0.0,0.0,0.0]
#!  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)

#!  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
#!  NUMBER_OF_DOMAINS=1


#!  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
#!  CALL MPI_BCAST(NUMBER_GLOBAL_Y_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
#!  CALL MPI_BCAST(NUMBER_GLOBAL_Z_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
#!  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
ProgressDiagnostics = False   # Set to diagnostics

#================================================================================================================================
#  Coordinate System
#================================================================================================================================

#!    !Start the creation of a new RC coordinate system
#!    CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
#!    CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
#!    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
#!      !Set the coordinate system to be 2D
#!      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
#!    ELSE
#!      !Set the coordinate system to be 3D
#!      CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,3,Err)
#!    ENDIF
#!    !Finish the creation of the coordinate system
#!   CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)


if (ProgressDiagnostics):
    print " == >> COORDINATE SYSTEM << == "
    
# Start the creation of RC coordinate system    
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()

#================================================================================================================================
#  Region
#================================================================================================================================

#!    !Start the creation of the region
#!    CALL cmfe_Region_Initialise(Region,Err)
#!    CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
#!    CALL cmfe_Region_LabelSet(Region,"MuscleRegion",Err)
#!    !Set the regions coordinate system to the 2D RC coordinate system that we have created
#!    CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
#!    !Finish the creation of the region
#!    CALL cmfe_Region_CreateFinish(Region,Err)
    
if (ProgressDiagnostics):
    print " == >> COORDINATE SYSTEM << == "

# Start the creation of SPACE region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "DiffusionRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()    

#================================================================================================================================
#  Bases
#================================================================================================================================

#!    !Start the creation of a basis (default is trilinear lagrange)
#!    CALL cmfe_Basis_Initialise(Basis,Err)
#!    CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
#!    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
#!      !Set the basis to be a bilinear Lagrange basis
#!      CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
#!    ELSE
#!      !Set the basis to be a trilinear Lagrange basis
#!      CALL cmfe_Basis_NumberOfXiSet(Basis,3,Err)
#!    ENDIF
#!    !Finish the creation of the basis
#!    CALL cmfe_Basis_CreateFinish(BASIS,Err)

if (ProgressDiagnostics):
    print " == >> BASIS << == "

# Create a tri-linear lagrange basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
#basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.type = iron.BasisTypes.SIMPLEX  
#if (numberGlobalZElements == 0):
#    basis.numberOfXi = 2
#    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*2
#    basis.QuadratureOrderSet(2)
#else:
basis.numberOfXi = 3	
#    basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*3
basis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]*3
basis.QuadratureOrderSet(2)

basis.CreateFinish()

#================================================================================================================================
#  Mesh
#================================================================================================================================


#!    !Start the creation of a generated mesh in the region
#!    CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
#!    CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
#!    !Set up a regular x*y*z mesh
#!    CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
#!    !Set the default basis
#!    CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
#!    !Define the mesh on the region
#!    IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
#!      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT],Err)
#!      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS],Err)
#!    ELSE
#!      CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[WIDTH,HEIGHT,LENGTH],Err)
#!      CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS,NUMBER_GLOBAL_Y_ELEMENTS, &
#!        & NUMBER_GLOBAL_Z_ELEMENTS],Err)
#!    ENDIF    
#!    !Finish the creation of a generated mesh in the region
#!    CALL cmfe_Mesh_Initialise(Mesh,Err)
#!    CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)

if (ProgressDiagnostics):
    print " == >> MESH << == "
    
# Create a generated mesh
#generatedMesh = iron.GeneratedMesh()
#generatedMesh.CreateStart(generatedMeshUserNumber,region)
#generatedMesh.type = iron.GeneratedMeshTypes.REGULAR
#generatedMesh.basis = [basis]
#if (numberGlobalZElements == 0):
#	generatedMesh.extent = [WIDTH,HEIGHT]
#	generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements]
#else:
#	generatedMesh.extent = [WIDTH,HEIGHT,LENGTH]
#	generatedMesh.numberOfElements = [numberGlobalXElements,numberGlobalYElements,numberGlobalZElements]

#mesh = iron.Mesh()
#generatedMesh.CreateFinish(meshUserNumber,mesh)



FileName_ele = "MaxVol1000/WholeLeftArm.1.ele"
numberOfNodes    = 40813
numberOfElements = 237070
numberOfLocalNodes = 4
offset = 1 # offset is 1 if nodes and elements number begin from 0, offset for default = 0

printHead = True
printTail = True
head = range(5)
tail = range(numberOfElements-4,numberOfElements+1)

muscleRegionLabel      = 2
leftRadiusRegionLabel  = 3
leftHumerusRegionLabel = 4
leftUlnaRegionLabel    = 5

muscleElements      = []
leftRadiusElements  = []
leftHumerusElements = []
leftUlnaElements    = []

mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber,region,3)
mesh.origin = meshOrigin
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(numberOfElements)
# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region,numberOfNodes)
nodes.CreateFinish()

elements = iron.MeshElements()
meshComponentNumber = 1
elements.CreateStart(mesh, meshComponentNumber, basis)


localNodes=[0]*numberOfLocalNodes
#elements.NodesSet(elementNumber,localNodes)
elementNumber=0
totalNumberOfElements=1
elementNodes = False
#with open("upperlimb_refined2.exelem", "r") as f:
print "Elapsed time before reading ele file is: ", time.time()-t
# reading elements and its local Nodes and setting elements.Nodes
with open(FileName_ele, "r") as f:
    target=f.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
    for lineNum,line in enumerate(target):
        if (lineNum !=0 and lineNum<=numberOfElements):
# reading elements and its localNodes
            elementNumber = int(target[lineNum][0])+offset
            localNodes[0] = int(target[lineNum][1])+offset
            localNodes[1] = int(target[lineNum][2])+offset
            localNodes[2] = int(target[lineNum][3])+offset
            localNodes[3] = int(target[lineNum][4])+offset
            materialType  = int(target[lineNum][5])
# Giving elements material label
            if materialType == muscleRegionLabel:
                muscleElements.append(elementNumber)
            elif materialType == leftRadiusRegionLabel:
                leftRadiusElements.append(elementNumber)
            elif materialType == leftHumerusRegionLabel:
                leftHumerusElements.append(elementNumber)
            elif materialType == leftUlnaRegionLabel:
                leftUlnaElements.append(elementNumber)
# printing the head and the tail of the elements
            if elementNumber in head and printHead:
              print elementNumber,localNodes
            elif elementNumber in tail and printTail:
              print elementNumber,localNodes

            elements.NodesSet(elementNumber,localNodes)

print "number of muscle Elements = %d\nnumber of right radius elements = %d\ntotal number of elements = %d\n"%(len(muscleElements)
,len(leftRadiusElements),len(muscleElements)+len(leftRadiusElements)+len(leftHumerusElements)+len(leftUlnaElements))

print "Elapsed time after reading ele file is: ", time.time()-t
raw_input("Press Enter to continue...")


elements.CreateFinish()
mesh.CreateFinish()

#================================================================================================================================
#  Decomposition
#================================================================================================================================
    

#!    !Create a decomposition
#!    CALL cmfe_Decomposition_Initialise(Decomposition,Err)
#!    CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
#!    !Set the decomposition to be a general decomposition with the specified number of domains
#!    CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
#!    CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
#!    CALL cmfe_Decomposition_CalculateFacesSet(decomposition,.TRUE.,err)    
#!    !Finish the decomposition
#!    CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

if (ProgressDiagnostics):
    print " == >> MESH DECOMPOSITION << == "

# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.calculateFaces = True  
decomposition.CreateFinish()

#================================================================================================================================
#  Geometric Field
#================================================================================================================================


#!  !Start to create a default (geometric) field on the region
#!  CALL cmfe_Field_Initialise(GeometricField,Err)
#!  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
#!  !Set the decomposition to use
#!  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
#!  !Set the domain to be used by the field components.
#!  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
#!  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
#!  IF(NUMBER_GLOBAL_Z_ELEMENTS/=0) THEN
#!    CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
#!  ENDIF
#!  !Finish creating the field
#!  CALL cmfe_Field_CreateFinish(GeometricField,Err)

       
#!    !Update the geometric field parameters
#!    CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)
#!!  ENDIF

#!  !IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD)) GEOMETRIC_FIELD=>REGION%FIELDS%FIELDS(1)%PTR

if (ProgressDiagnostics):
    print " == >> GEOMETRIC FIELD << == "
     
# Create a field for the geometry
#geometricField = iron.Field()
#geometricField.CreateStart(geometricFieldUserNumber,region)
#geometricField.meshDecomposition = decomposition
#geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
#geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
#if (numberGlobalZElements != 0):
#	geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
#geometricField.CreateFinish()

# Set geometry from the generated mesh
#generatedMesh.GeometricParametersCalculate(geometricField)



geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
#if (numberGlobalZElements != 0):
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)

geometricField.CreateFinish()

# Set geometry from the generated mesh
#generatedMesh.GeometricParametersCalculate(geometricField)

# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes
print numberOfNodes
# Get or calculate geometric Parameters

print "Elapsed time before reading node file is: ", time.time()-t

FileName_node = "MaxVol1000/correctedBoundaryMarker.node"

printHead = True
printTail = True
head = range(5)
tail = range(numberOfNodes-4,numberOfNodes+1)

#nodePositions = False
boundaryMarker=0
nodeNumber = 0
#numberOfNodes = 1
#frontBoundary  = []
#backBoundary   = []
#bottomBoundary = []
#rightBoundary  = []
#topBoundary    = []
#leftBoundary   = []
boundary = []
skinMarker = 2

with open(FileName_node, "r") as f:
    target=f.readlines()
    for lineNum,line in enumerate(target):
        target[lineNum] = target[lineNum].rstrip('\n\r').replace("\t"," ").replace(","," ").split()
    for lineNum,line in enumerate(target):
        if lineNum !=0:
          nodeNumber = int(target[lineNum][0])+offset
          x = float(target[lineNum][1])
          y = float(target[lineNum][2])
          z = float(target[lineNum][3])
          boundaryMarker = int(target[lineNum][4])
          if boundaryMarker == skinMarker:
            boundary.append(nodeNumber)

          if nodeNumber in head and printHead:
            print nodeNumber,x,y,z,boundaryMarker
          elif nodeNumber in tail and printTail:
            print nodeNumber,x,y,z,boundaryMarker

          geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,1,x)
          geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,2,y)
          geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,nodeNumber,3,z)
#        if 'POINTS' in line:
#          nodePositions = True
#          numberOfNodes = int(target[lineNum][1])
print "number of boundary nodes = %d"%len(boundary)

print "Elapsed time after reading node file is: ", time.time()-t

raw_input("Press Enter to continue...")

# Update the geometric field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

#================================================================================================================================
#  Equations Sets
#================================================================================================================================

  
#!  !Create the equations_set
#!  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
#!    CALL cmfe_Field_Initialise(EquationsSetField,Err)
#!CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
#!  & CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE], &
#!  & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
#!  !Set the equations set to be a standard Laplace problem
  
#!  !Finish creating the equations set
#!  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

if (ProgressDiagnostics):
    print " == >> EQUATIONS SET << == "

# Create standard Diffusion equations set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.CLASSICAL_FIELD,
        iron.EquationsSetTypes.DIFFUSION_EQUATION,
        iron.EquationsSetSubtypes.CONSTANT_SOURCE_DIFFUSION]
equationsSet.CreateStart(equationsSetUserNumber,region,geometricField,
        equationsSetSpecification,equationsSetFieldUserNumber,equationsSetField)
equationsSet.CreateFinish()

#================================================================================================================================
#  Dependent Field
#================================================================================================================================


#!  !Create the equations set dependent field variables
#!  CALL cmfe_Field_Initialise(DependentField,Err)
#!  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
#!  !Finish the equations set dependent field variables
#!  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)
  
#!  !Initialise the field with an initial guess
#!  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP, &
#!    & Err)  
#!SELECT CASE(1)
#!CASE(1)
#!  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP, &
#!    & Err)  
#!CASE(2)
#!  j = (NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
#!  k = NUMBER_GLOBAL_Y_ELEMENTS/2*(NUMBER_GLOBAL_X_ELEMENTS+1)+NUMBER_GLOBAL_X_ELEMENTS/2+1

#!  Do i = 1,NUMBER_GLOBAL_Z_ELEMENTS/2 
#!    CALL  cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,k+(i-1)*j,1,37.0_CMISSRP,Err)
#!  END DO
#!  k = k+(NUMBER_GLOBAL_Z_ELEMENTS/2-1)*j
#!  DO i = 1,NUMBER_GLOBAL_Z_ELEMENTS/2
#!    CALL  cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,k+i*(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+2),1,37.0_CMISSRP,Err)
#!    CALL  cmfe_Field_ParameterSetUpdateNode(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,k+i*(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS),1,37.0_CMISSRP,Err)   
#!  END DO
#!END SELECT

if (ProgressDiagnostics):
    print " == >> DEPENDENT FIELD << == "
    
# Create dependent field
dependentField = iron.Field()
equationsSet.DependentCreateStart(dependentFieldUserNumber,dependentField)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.U,iron.FieldDOFOrderTypes.SEPARATED)
dependentField.DOFOrderTypeSet(iron.FieldVariableTypes.DELUDELN,iron.FieldDOFOrderTypes.SEPARATED)
equationsSet.DependentCreateFinish()

# Initialise dependent field
dependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)

#================================================================================================================================
#  Materials Field
#================================================================================================================================

#!  !Create the equations set material field variables
#!  CALL cmfe_Field_Initialise(MaterialsField,Err)
#!  CALL cmfe_EquationsSet_MaterialsCreateStart#!(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
#!  !equation is a(x)*delTdelt + k(x)*grad(T) + c =0
#!!  MaterialsField%field%variables(1)%number_of_components = 3
#!  !Finish the equations set dependent field variables
#!!  MaterialsField%field%variables(1)%parameter_sets%parameter_sets(1)%ptr%parameters%cmiss%data_dp(1) =9009.99_DP
#!  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)
  
#!    !Initialise the field with an initial guess
#!    K_Conductivity  = 100.0_CMISSRP
#!!    Rho_Density     = 1050.0_CMISSRP
#!!    Cp_SpecificHeat = 3800.0_CMISSRP
#!    K_Conductivity  = 1.0_CMISSRP
#!    Rho_Density     = 1.0_CMISSRP
#!    Cp_SpecificHeat = 1.0_CMISSRP
#!    alpha = K_Conductivity/(Rho_Density*Cp_SpecificHeat) !Diffusivity

    
#!  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,alpha, &
#!    & Err) 
#!  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,alpha, &
#!    & Err)     

if (ProgressDiagnostics):
    print " == >> MATERIALS FIELD << == "

#!  !Create the equations set material field variables
materialsField = iron.Field()
equationsSet.MaterialsCreateStart(materialsFieldUserNumber,materialsField)
equationsSet.MaterialsCreateFinish()

K_Conductivity  = 1.0
Rho_Density     = 1.0
Cp_SpecificHeat = 1.0
alpha = K_Conductivity/(Rho_Density*Cp_SpecificHeat) #Diffusivity
#materialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,alpha) 
#materialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,alpha) 

#materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, userElementNumber,componentNumber, value)
#componentNumber = 1 conductivity
#componentNumber = 2 conductivity
#componentNumber = 3 conductivity
#componentNumber = 4 rho*cp
for elementNumber in muscleElements:
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,1, 0.42*1e-3)
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,2, 0.42*1e-3)
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,3, 0.42*1e-3)
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,4, 4.0e6*1e-9)

for elementNumber in leftRadiusElements:
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,1, 0.75*1e-3)
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,2, 0.75*1e-3)
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,3, 0.75*1e-3)
    materialsField.ParameterSetUpdateElement(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, 
elementNumber,4, 2.3e6*1e-9)
#================================================================================================================================
#  Source Field
#================================================================================================================================
  
#!  !Create the equations set source field variables
#!  CALL cmfe_Field_Initialise(SourceField,Err)
#!  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
#!  !Finish the equations set dependent field variables
#!  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)


#!  CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
#!    & 1,50.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)
    
#!Select CASE(1)
#!CASE(1)
#!  CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, & 
#!    & 1,50.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)
#!CASE(2)    
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,56,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,57,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)  
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,58,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)  
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,59,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)          
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,60,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)  
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,72,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,84,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)          
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,96,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)  
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,108,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)                                                                 
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,50,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,40,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)          
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,30,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)  
#!  CALL cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,20,1,48000000.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)       
#!CASE(3)
#!! ARM 3D. provisional goal
#!  j = (NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
#!  k = NUMBER_GLOBAL_Y_ELEMENTS/2*(NUMBER_GLOBAL_X_ELEMENTS+1)+NUMBER_GLOBAL_X_ELEMENTS/2+1

#!  Do i = 1,NUMBER_GLOBAL_Z_ELEMENTS/2 
#!    CALL  cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,k+(i-1)*j,1,50.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)    
#!  END DO
#!  k = k+(NUMBER_GLOBAL_Z_ELEMENTS/2-1)*j
#!  DO i = 1,NUMBER_GLOBAL_Z_ELEMENTS/2
#!    CALL  cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,k+i*(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+2), &
#!        & 1,50.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)
#!    CALL  cmfe_Field_ParameterSetUpdateNode(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, & 
#!        & CMFE_NO_GLOBAL_DERIV,k+i*(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS), &
#!        & 1,50.0_CMISSRP/(Rho_Density*Cp_SpecificHeat),Err)   
#!  END DO
#!END SELECT
#!  !Create the equations set analytic field variables
#!!  CALL cmfe_Field_Initialise(AnalyticField,Err)
#!!  CALL cmfe_EquationsSet_AnalyticCreateStart(EquationsSet,CMFE_EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1, &
#!!    & AnalyticFieldUserNumber, &
#!!    & AnalyticField,Err)
#!  !Finish the equations set analytic field variables
#!!  CALL cmfe_EquationsSet_AnalyticCreateFinish(EquationsSet,Err)

if (ProgressDiagnostics):
    print " == >> SOURCE FIELD << == "

sourceField = iron.Field()
equationsSet.SourceCreateStart(sourceFieldUserNumber,sourceField)
equationsSet.SourceCreateFinish()  

#sourceField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,50.0/(Rho_Density*Cp_SpecificHeat))
sourceField.ComponentValuesInitialise(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0.0)
#================================================================================================================================
#  Equations
#================================================================================================================================

#!  !Create the equations set equations
#!  CALL cmfe_Equations_Initialise(Equations,Err)
#!  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
#!  !Set the equations matrices sparsity type
#!  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
#!  !Set the equations set output
#!  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
#!  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_TIMING_OUTPUT,Err)
#!  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_MATRIX_OUTPUT,Err)
#!  !CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT,Err)
#!  !Finish the equations set equations
#!  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)


  
#!  !Create the equations set boundary conditions
#!!   CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
#!!   CALL cmfe_EquationsSetBoundaryConditionsCreateStart(EquationsSet,BoundaryConditions,Err)
#!!   !Set the first node to 0.0 and the last node to 1.0
#!! !   FirstNodeNumber=1
#!! !   IF(NUMBER_GLOBAL_Z_ELEMENTS==0) THEN
#!! !     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)
#!! !   ELSE
#!! !     LastNodeNumber=(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
#!! !   ENDIF
#!! !   CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_U_VARIABLE_TYPE,1,FirstNodeNumber,1, &
#!! !     & CMFE_BOUNDARY_CONDITION_FIXED,0.0_CMISSRP,Err)
#!! !   CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,LastNodeNumber,1, &
#!! !     & CMFE_BOUNDARY_CONDITION_FIXED,1.0_CMISSRP,Err)
#!!   !Finish the creation of the equations set boundary conditions
#!!   CALL cmfe_EquationsSetBoundaryConditionsCreateFinish(EquationsSet,Err)

if (ProgressDiagnostics):
    print " == >> EQUATIONS << == "

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.sparsityType = iron.EquationsSparsityTypes.SPARSE
equations.outputType = iron.EquationsOutputTypes.NONE
equationsSet.EquationsCreateFinish()

#================================================================================================================================
#  Problems
#================================================================================================================================

#!  !Create the problem
#!  CALL cmfe_Problem_Initialise(Problem,Err)
#!  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS,CMFE_PROBLEM_DIFFUSION_EQUATION_TYPE, &
#!    & CMFE_PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE],Problem,Err)
#!  !Finish the creation of a problem.
#!  CALL cmfe_Problem_CreateFinish(Problem,Err)

if (ProgressDiagnostics):
    print " == >> PROBLEM << == "

# Create Laplace problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.CLASSICAL_FIELD,
        iron.ProblemTypes.DIFFUSION_EQUATION,
        iron.ProblemSubtypes.LINEAR_SOURCE_DIFFUSION] #?????????
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()

#================================================================================================================================
#  Control Loops
#================================================================================================================================

#!  !Create the problem control
#!  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
#!  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
#!  !Get the control loop
#!  CALL cmfe_Problem_ControlLoopGet(Problem,ControlLoopNode,ControlLoop,Err)
#!  !Set the times
#!  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,0.8001_CMISSRP,0.001_CMISSRP,Err)
#!  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,10,Err)    
#!  !Finish creating the problem control loop
#!  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)

if (ProgressDiagnostics):
    print " == >> PROBLEM CONTROL LOOP << == "

# Create control loops
problem.ControlLoopCreateStart()
TimeLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],TimeLoop)
TimeLoop.LabelSet('Time Loop')
TimeLoop.TimesSet(startTime,stopTime,timeIncrement)
TimeLoop.TimeOutputSet(DYNAMIC_SOLVER_DIFFUSION_OUTPUT_FREQUENCY)
problem.ControlLoopCreateFinish()

#================================================================================================================================
#  Solvers
#================================================================================================================================



#!  !Start the creation of the problem solvers
#!!  
#!! !   !For the Direct Solver MUMPS, uncomment the below two lines and comment out the above five
#!! !   CALL SOLVER_LINEAR_TYPE_SET(LINEAR_SOLVER,SOLVER_LINEAR_DIRECT_SOLVE_TYPE,ERR,ERROR,*999)
#!! !   CALL SOLVER_LINEAR_DIRECT_TYPE_SET(LINEAR_SOLVER,SOLVER_DIRECT_MUMPS,ERR,ERROR,*999) 
#!! 



#!  !Start the creation of the problem solvers
#!  CALL cmfe_Solver_Initialise(Solver,Err)
#!  CALL cmfe_Solver_Initialise(LinearSolver,Err)
#!  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
#!  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
#!  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
#!  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
#!  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
#!  !CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
#!  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
#!  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
#!  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)
#!  !Finish the creation of the problem solver
#!  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

# Create problem solver
solver = iron.Solver()
LinearSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
#solver.outputType = iron.SolverOutputTypes.SOLVER
solver.outputType = iron.SolverOutputTypes.PROGRESS
solver.DynamicLinearSolverGet(LinearSolver)
#solver.linearType = iron.LinearSolverTypes.ITERATIVE
#solver.LinearIterativeMaximumIterationsSet(MAXIMUM_ITERATIONS)
#solver.linearIterativeAbsoluteTolerance = relativeTolerance
#solver.linearIterativeRelativeTolerance = absoluteTolerance
problem.SolversCreateFinish()

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

#!  !Create the problem solver equations
#!  CALL cmfe_Solver_Initialise(Solver,Err)
#!  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
#!  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
#!  !Get the solve equations
#!  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
#!  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
#!  !Set the solver equations sparsity
#!  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)
#!  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_FULL_MATRICES,Err)  
#!  !Add in the equations set
#!  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
#!  !Finish the creation of the problem solver equations
#!  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

if (ProgressDiagnostics):
    print " == >> SOLVER EQUATIONS << == "

# Create solver equations and add equations set to solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

  
#!  ! Grouping Elements
#!!  allocate( ElementUserNodes(8) )
#!!  CALL cmfe_MeshElements_NodesGet(regionUserNumber,meshUserNumber,1,64, &
#!!      & elementUserNodes,err)
#!!  Top = 
#!!  CALL cmfe_MeshElements_AdjacentElementGet(regionUserNumber,meshUserNumber,1,64, & 
#!!      & 2,adjacentElement,err)
#!!  if (adjacentElemnt == 0) then
#!!    adjacentElent is top element
#!!  end if
#!!  print *, GeneratedMesh%generatedmesh%mesh%topology(1)%ptr%nodes%nodes(125)%boundarynode

#!!  MeshElements_NodesGet
  

  
#!  !Start the creation of the equations set boundary conditions
#!  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
#!  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)




#!case(4)
#! !Dirichlet symmetric 3D 8*4*4 elements 4*2*2 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,90.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,2,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,75.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,3,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,75.0_CMISSRP,Err)                 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,4,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,75.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,5,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,6,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,7,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                     
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)
#!!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,8,1, &
#!!      & CMFE_BOUNDARY_CONDITION_ROBIN,22.0_CMISSRP,Err)      
#!!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,9,1, &
#!!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,9,1, &
#!      & CMFE_BOUNDARY_CONDITION_NEUMANN_INTEGRATED,16.0_CMISSRP,Err)
                    
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,10,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,11,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,12,1, &
#!!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                     
#!!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,13,1, &
#!!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,12,1, &
#!		    & CMFE_BOUNDARY_CONDITION_NEUMANN_POINT,19.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,13,1, &
#!		    & CMFE_BOUNDARY_CONDITION_NEUMANN_INTEGRATED_ONLY,21.0_CMISSRP,Err)		    
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,14,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,15,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)              
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,16,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,17,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                     
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,18,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)  
#!!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,18,1, &
#!!      & CMFE_BOUNDARY_CONDITION_DIRICHLET,27.0_CMISSRP,Err)        
#!!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,19,1, &
#!!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)               
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,19,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,20,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,21,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,22,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,150.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,23,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,150.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,24,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,150.0_CMISSRP,Err)   
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,25,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)                
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,26,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,27,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,28,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,29,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,30,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,31,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,35,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,36,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,40,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
           
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,41,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,45,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,46,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,47,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,48,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,49,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,50,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err)                  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,51,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,52,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,53,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,54,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,55,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,56,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,60,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,61,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,65,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,66,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,70,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,71,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,72,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,73,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,74,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,75,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,76,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,77,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,78,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,79,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,80,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,81,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,85,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,86,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,90,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,91,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,95,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,96,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,97,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,98,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,99,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,100,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,101,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,102,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,103,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,104,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,105,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,106,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,110,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,111,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,115,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,116,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,120,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,121,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,122,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,123,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,124,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,125,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,126,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,127,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,128,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,129,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,130,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,131,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,135,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,136,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,140,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,141,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,145,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,146,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,147,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,148,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,149,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,150,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,151,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,152,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,153,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,154,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,155,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,156,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,160,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,161,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,165,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,166,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,170,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,171,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,172,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,173,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,174,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,175,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err) 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,176,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,85.0_CMISSRP,Err)      
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,177,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,178,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)            
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,179,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,50.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,180,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,115.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,181,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,185,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,186,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,190,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,191,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,120.0_CMISSRP,Err)                         
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,195,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,180.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,196,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,197,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,198,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,199,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,200.0_CMISSRP,Err)             
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,200,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,190.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,201,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,90.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,202,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,75.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,203,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,75.0_CMISSRP,Err)                 
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,204,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,75.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,205,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,206,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,207,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                     
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,208,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,209,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)              
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,210,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,211,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,212,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                     
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,213,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,214,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)              
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,215,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)              
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,216,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,110.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,217,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)                     
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,218,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)  
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,219,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,100.0_CMISSRP,Err)              
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,220,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)        
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,221,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,140.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,222,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,150.0_CMISSRP,Err)       
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,223,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,150.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,224,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,150.0_CMISSRP,Err)
#!  CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,225,1, &
#!      & CMFE_BOUNDARY_CONDITION_FIXED,160.0_CMISSRP,Err) 

if (ProgressDiagnostics):
    print " == >> BOUNDARY CONDITIONS << == "

boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

#firstNodeNumber=1
#nodes = iron.Nodes()
#region.NodesGet(nodes)
#lastNodeNumber = nodes.numberOfNodes
#firstNodeDomain = decomposition.NodeDomainGet(firstNodeNumber,1)
#lastNodeDomain = decomposition.NodeDomainGet(lastNodeNumber,1)
#if firstNodeDomain == computationalNodeNumber:
#    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,firstNodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
#if lastNodeDomain == computationalNodeNumber:
#    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,lastNodeNumber,1,iron.BoundaryConditionsTypes.FIXED,1.0)
             
#boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,versioIdx,derivIdx,nodeNumber,componentIdx,iron.BoundaryConditionsTypes.FIXED,90.0)

#boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,1,1,iron.BoundaryConditionsTypes.ROBIN,[90.0,100.0,110.0])
#boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,1,1,iron.BoundaryConditionsTypes.FIXED,[90.0])
#for i in range (2,(numberGlobalXElements+1)):
#    boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,i,1,iron.BoundaryConditionsTypes.FIXED,[75.0]) 
#boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,numberGlobalXElements+1,1,iron.BoundaryConditionsTypes.NEUMANN_POINT,[110.0])
#boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,5,1,iron.BoundaryConditionsTypes.FIXED,[110.0])

#boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.DELUDELN,1,1,6,1,iron.BoundaryConditionsTypes.ROBIN,[110.0,120.0,110.0])
for nodeNumber in boundary:
  boundaryConditions.SetNode(dependentField,iron.FieldVariableTypes.U,1,1,nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,[100.0])


#!  !Finish the creation of the equations set boundary conditions
#!  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)	

solverEquations.BoundaryConditionsCreateFinish()

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#!  !Solve the problem
#!  CALL cmfe_Problem_Solve(Problem,Err)
print "Solving problem..."
start = time.time()
# Solve the problem
problem.Solve()
end = time.time()
elapsed = end - start
print "Total Number of Elements = %d " %totalNumberOfElements
print "Total Number of nodes = %d " %numberOfNodes
print "Calculation Time = %3.4f" %elapsed
print "Problem solved!"
print "#"
print "number of muscle Elements = %d\nnumber of right radius elements = %d\ntotal number of elements = %d\n"%(len(muscleElements)
,len(leftRadiusElements),len(muscleElements)+len(leftRadiusElements))
print "number of boundary nodes = %d"%len(boundary)
print "Elapsed time: ", time.time()-t
#!# Export results
#!baseName = "laplace"
#!dataFormat = "PLAIN_TEXT"
#!fml = iron.FieldMLIO()
#!fml.OutputCreate(mesh, "", baseName, dataFormat)
#!fml.OutputAddFieldNoType(baseName+".geometric", dataFormat, geometricField,
#!    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#!fml.OutputAddFieldNoType(baseName+".phi", dataFormat, dependentField,
#!    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
#!fml.OutputWrite("LaplaceExample.xml")
#!fml.Finalise()


#!  !Output Analytic analysis
#!!  CALL cmfe_AnalyticAnalysis_Output(DependentField,"DiffusionAnalytics_x10_y10_z10_L_T1",Err)
#exportField = True
#if (exportField): 
#    fields = iron.Fields()
#    fields.Create(Region)
#    fields.Finalise()
#!  EXPORT_FIELD=.TRUE.
#!  IF(EXPORT_FIELD) THEN
#!    CALL cmfe_Fields_Initialise(Fields,Err)
#!    CALL cmfe_Fields_Create(Region,Fields,Err)
#!!    CALL cmfe_Fields_NodesExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#!!    CALL cmfe_Fields_ElementsExport(Fields,"DiffusionConstantSourceAnalytic_x10_y10_z10_L_T1","FORTRAN",Err)
#!    CALL cmfe_Fields_Finalise(Fields,Err)

#!  ENDIF

#!  !CALL cmfe_Finalise(Err)
#!  WRITE(*,'(A)') "Program successfully completed."
  
#!  STOP
iron.Finalise() 
#!END PROGRAM DIFFUSIONCONSTANTSOURCEEXAMPLE
