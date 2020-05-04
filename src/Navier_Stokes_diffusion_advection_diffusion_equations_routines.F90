!> \file
!> \authors Elias Soltani
!> \brief This module handles all routines pertaining to diffusion coupled to diffusion.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand, the University of Oxford, Oxford, United
!> Kingdom and King's College, London, United Kingdom. Portions created
!> by the University of Auckland, the University of Oxford and King's
!> College, London are Copyright (C) 2007-2010 by the University of
!> Auckland, the University of Oxford and King's College, London.
!> All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!>TThis module handles all routines pertaining to diffusion coupled to diffusion.


MODULE NAVIER_STOKES_DIFFUSION_ADVECTION_DIFFUSION_ROUTINES

  USE ADVECTION_DIFFUSION_EQUATION_ROUTINES
  USE ADVECTION_EQUATION_ROUTINES
  USE BaseRoutines
  USE BasisRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CHARACTERISTIC_EQUATION_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE ComputationEnvironment
  USE CmissMPI
  USE DIFFUSION_EQUATION_ROUTINES
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsSetConstants
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
#ifndef NOMPIMOD
  USE MPI
#endif
  USE MatrixVector
  USE MESH_ROUTINES
  USE MeshAccessRoutines
  USE NAVIER_STOKES_EQUATIONS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverMappingAccessRoutines
  USE SolverMatricesAccessRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"


  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  TYPE :: SCHEDULE_TYPE
    REAL(DP) :: start
    INTEGER(INTG) :: timeIncrement
    INTEGER(INTG) :: timeSteps
    INTEGER(INTG) :: outputFrequency
    REAL(DP) :: convectionCoeff
    REAL(DP) :: heatFlux
    REAL(DP) :: Tair
    REAL(DP) :: Pv
  END TYPE SCHEDULE_TYPE

  PUBLIC NavierStokesDiffAdvDiff_EquationsSetSetup
  PUBLIC NavierStokesDiffAdvDiff_EquationsSetSpecSet
  PUBLIC NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet

  PUBLIC NavierStokesDiffAdvDiff_ProblemSetup
  PUBLIC NavierStokesDiffAdvDiff_ProblemSpecificationSet

  PUBLIC NavierStokesDiffAdvDiff_FiniteElementCalculate

  PUBLIC NavierStokesDiffAdvDiff_PreSolve
  PUBLIC NavierStokesDiffAdvDiff_PostSolve

  PUBLIC NavierStokesDiffAdvDiff_ControlLoopPreLoop
  PUBLIC NavierStokesDiffAdvDiff_ControlLoopPostLoop

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a coupled navier stokes diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "diffusion and advection-diffusion type equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a diffusion & advection-diffusion equation type of a multi physics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet",ERR,ERROR)
    EXITS("NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes & diffusion & advection-diffusion coupled equation.
  SUBROUTINE NavierStokesDiffAdvDiff_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables


    ENTERS("NavierStokesDiffAdvDiff_EquationsSetSetup",ERR,ERROR,*999)

    CALL FlagError("Not implemented.",ERR,ERROR,*999)

    EXITS("NavierStokesDiffAdvDiff_EquationsSetSetup")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_EquationsSetSetup",ERR,ERROR)
    EXITS("NavierStokesDiffAdvDiff_EquationsSetSetup")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a coupled Navier Stokes & diffusion & advection-diffusion equation finite element equations set.
  SUBROUTINE NavierStokesDiffAdvDiff_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("NavierStokesDiffAdvDiff_FiniteElementCalculate",ERR,ERROR,*999)

    CALL FlagError("Not implemented.",ERR,ERROR,*999)

    EXITS("NavierStokesDiffAdvDiff_FiniteElementCalculate")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_FiniteElementCalculate",ERR,ERROR)
    EXITS("NavierStokesDiffAdvDiff_FiniteElementCalculate")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE NavierStokesDiffAdvDiff_EquationsSetSpecSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("NavierStokesDiffAdvDiff_EquationsSetSpecSet",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)

    EXITS("NavierStokesDiffAdvDiff_EquationsSetSpecSet")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_EquationsSetSpecSet",err,error)
    EXITS("NavierStokesDiffAdvDiff_EquationsSetSpecSet")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_EquationsSetSpecSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a coupled Navier-Stokes & diffusion & advection-diffusion problem.
  SUBROUTINE NavierStokesDiffAdvDiff_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("NavierStokesDiffAdvDiff_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
          & PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a coupled Navier-Stokes & diffusion & advection-diffusion type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_NAVIER_STOKES_DIFFUSION_ADVECTION_DIFFUSION_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Navier-Stokes diffusion advection-diffusion problem specification must have 3 entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokesDiffAdvDiff_ProblemSpecificationSet")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_ProblemSpecificationSet",err,error)
    EXITS("NavierStokesDiffAdvDiff_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the coupled Navier-Stokes diffusion-diffusion equations problem.
  SUBROUTINE NavierStokesDiffAdvDiff_ProblemSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT,simpleLoopCellML1,whileLoopConditional,simpleLoop2, &
      & iterativeWhileLoop,timeLoopFlow
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION, SOLVER_ADVECTION_DIFFUSION,cellMLSolver,solverAdvectionDiffusion, &
      & solverDiffusion, solver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,solverEquations
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS,cellMLEquations
    TYPE(VARYING_STRING) :: LOCAL_ERROR,localError

    ENTERS("NavierStokesDiffAdvDiff_ProblemSetup",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_DIFFUSION)
    NULLIFY(SOLVER_ADVECTION_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a coupled Navier-Stokes & diffusion & advection-diffusion problem.", &
          & err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   coupled  source diffusion--advection-diffusion
      !--------------------------------------------------------------------
      CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,3,ERR,ERROR,*999)
            !
            !Set the first solver to be a linear solver for the advection-diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !
            !Set the second solver to be a linear solver for the diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_DIFFUSION,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !

!!!-- CELLML EVALUATOR SOLVER --!!!
            NULLIFY(cellMLSolver)
            CALL SOLVERS_SOLVER_GET(solvers,1,cellMLSolver,err,error,*999)
            CALL SOLVER_TYPE_SET(cellMLSolver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(cellMLSolver,"Evaluator Solver",err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !
            !Get the advection-diffusion solver and create the advection-diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION, &
              & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            !Get the diffusion solver and create the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_DIFFUSION,SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_DIFFUSION,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_DIFFUSION, &
              & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_DIFFUSION,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !

          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the advection-diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            !
            !Finish the creation of the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_DIFFUSION,SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            !

          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Elias \*
          !Create the CELLML solver equations
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
              ! NULLIFY(iterativeWhileLoop)
              ! CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              ! NULLIFY(simpleLoop)
              ! CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              ! CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
            ! ELSE
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            END IF
             ! NULLIFY(solver_cellML)
             NULLIFY(cellMLSolver)
             NULLIFY(CELLML_EQUATIONS)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,cellMLSolver,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_START(cellMLSolver,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_STATIC,err,error,*999)
              !Set the linearity
              CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_LINEAR,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for cellML equations setup Navier-Stokes equation type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            !   NULLIFY(iterativeWhileLoop)
            !   CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
            !   NULLIFY(simpleLoop)
            !   CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            !   CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
            ! ELSE
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            END IF
            NULLIFY(cellMLSolver)
            NULLIFY(CELLML_EQUATIONS)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,cellMLSolver,err,error,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(cellMLSolver,CELLML_EQUATIONS,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "// &
                & TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for cellML equations setup Navier-Stokes fluid mechanics problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a CellML setup for a 1D Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT    !Elias */
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a coupled diffusion & advection-diffusion equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT


      !--------------------------------------------------------------------
      !   coupled Navier-Stokes & source diffusion--advection-diffusion
      !--------------------------------------------------------------------
      CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled Navier-Stokes diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up the root time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,3,err,error,*999)
            CALL CONTROL_LOOP_LABEL_SET(CONTROL_LOOP,"main time loop",err,error,*999)
            !Setup the simple loop for shivering cellML sovler
            NULLIFY(simpleLoopCellML1)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,simpleLoopCellML1,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(simpleLoopCellML1,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL CONTROL_LOOP_LABEL_SET(simpleLoopCellML1,"Shivering CellML solver Loop",err,error,*999)
            !Setup a control loop for solving flow whenever the flag continue loop is on
            NULLIFY(whileLoopConditional)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,whileLoopConditional,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(whileLoopConditional,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(whileLoopConditional,1000,err,error,*999) !I do not need this probably
            CALL ControlLoop_AbsoluteToleranceSet(whileLoopConditional,0.1_DP,err,error,*999)   !I do not need this probalby
            CALL CONTROL_LOOP_LABEL_SET(whileLoopConditional,"Conditional while loop",err,error,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(whileLoopConditional,1,err,error,*999)
            !Setup the simple loop for shivering cellML sovler
            NULLIFY(simpleLoop2)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(simpleLoop2,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL CONTROL_LOOP_LABEL_SET(simpleLoop2,"Diffusion advection-diffusion solver Loop",err,error,*999)
            !Set up a time loop for flow solver
            NULLIFY(timeLoopFlow)
            CALL CONTROL_LOOP_SUB_LOOP_GET(whileLoopConditional,1,timeLoopFlow,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(timeLoopFlow,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(timeLoopFlow,1,err,error,*999)
            CALL CONTROL_LOOP_LABEL_SET(timeLoopFlow,"Flow time loop",err,error,*999)
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop. TODO: Add another subtype for RCR or 1D0D subtype and create the control loop and solvers for it.
            NULLIFY(iterativeWhileLoop)
            CALL CONTROL_LOOP_SUB_LOOP_GET(timeLoopFlow,1,iterativeWhileLoop,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E3_DP,err,error,*999)
            CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled Navier-Stokes & diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !Create the solvers
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !We need 5 solvers: 1 CELLML evalueator solver for shivering, 2 solvers for solving flow i.e. characteristic nonlinear solver and 1D NS solver, and 1
            !solver for solving energy equation in blood vessels and another one to solve transient diffusion equation for bioheat equation.

            !Start the solvers creation
            NULLIFY(solvers)
            NULLIFY(cellMLSolver)
            NULLIFY(simpleLoopCellML1)
            !Time loop has 3 sub loops, simple loop, while loop and another simple loop

            !Simple loop for cellML solver
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,simpleLoopCellML1,err,error,*999)
            CALL SOLVERS_CREATE_START(simpleLoopCellML1,solvers,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- CELLML EVALUATOR SOLVER --!!!
            CALL SOLVERS_SOLVER_GET(solvers,1,cellMLSolver,err,error,*999)
            CALL SOLVER_TYPE_SET(cellMLSolver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(cellMLSolver,"Evaluator Solver",err,error,*999)

            !Conditional while loop
            NULLIFY(whileLoopConditional)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,whileLoopConditional,err,error,*999)
            !Flow time loop
            NULLIFY(timeLoopFlow)
            CALL CONTROL_LOOP_SUB_LOOP_GET(whileLoopConditional,1,timeLoopFlow,err,error,*999)
            !Iterative while loop for Flow
            NULLIFY(iterativeWhileLoop)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SUB_LOOP_GET(timeLoopFlow,1,iterativeWhileLoop,err,error,*999)
            CALL SOLVERS_CREATE_START(iterativeWhileLoop,solvers,err,error,*999)
            CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
            NULLIFY(SOLVER)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
            CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)

            !Simple loop for solving vessels energy equation and bioheat equation in tissues
            NULLIFY(simpleLoop2)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL SOLVERS_CREATE_START(simpleLoop2,solvers,err,error,*999)
            CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- E N E R G Y  E Q U A T I O N  V E S S E L S --!!!
            !Set the first solver to be a linear solver for the advection-diffusion problem
            NULLIFY(solverAdvectionDiffusion)
            CALL SOLVERS_SOLVER_GET(solvers,1,solverAdvectionDiffusion,err,error,*999)
            CALL SOLVER_TYPE_SET(solverAdvectionDiffusion,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(solverAdvectionDiffusion,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(solverAdvectionDiffusion,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(solverAdvectionDiffusion,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(solverAdvectionDiffusion,SOLVER_CMISS_LIBRARY,err,error,*999)
!!!-- B I O H E A T --!!!
            !Set the second solver to be a linear solver for the diffusion problem
            NULLIFY(solverDiffusion)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,solverDiffusion,err,error,*999)
            CALL SOLVER_TYPE_SET(solverDiffusion,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(solverDiffusion,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(solverDiffusion,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(solverDiffusion,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(solverDiffusion,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            NULLIFY(solvers)
            NULLIFY(cellMLSolver)
            NULLIFY(simpleLoopCellML1)
            !Time loop has 3 sub loops, simple loop, while loop and another simple loop

            !Simple loop for cellML solver
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,simpleLoopCellML1,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(simpleLoopCellML1,solvers,err,error,*999)
            CALL SOLVERS_CREATE_FINISH(solvers,err,error,*999)

            !Conditional while loop
            NULLIFY(whileLoopConditional)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,whileLoopConditional,err,error,*999)
            !Flow time loop
            NULLIFY(timeLoopFlow)
            CALL CONTROL_LOOP_SUB_LOOP_GET(whileLoopConditional,1,timeLoopFlow,err,error,*999)
            !Iterative while loop for Flow
            NULLIFY(iterativeWhileLoop)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SUB_LOOP_GET(timeLoopFlow,1,iterativeWhileLoop,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,solvers,err,error,*999)
            CALL SOLVERS_CREATE_FINISH(solvers,err,error,*999)

            !Simple loop for solving vessels energy equation and bioheat equation in tissues
            NULLIFY(simpleLoop2)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop2,solvers,err,error,*999)
            CALL SOLVERS_CREATE_FINISH(solvers,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a coupled Navier-Stokes & diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        !Create the solver equations
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Control loops have 3 solvers and solvers have 4 solvers in total plus one cellML.
            !The first simple loop has a cellML evaluator. see CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE).
            !Get the solver for each solvers of control loops and create solver equations for them

            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Flow solvers
            NULLIFY(whileLoopConditional)
            NULLIFY(iterativeWhileLoop)
            NULLIFY(timeLoopFlow)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,whileLoopConditional,err,error,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(whileLoopConditional,1,timeLoopFlow,err,error,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(timeLoopFlow,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,solvers,err,error,*999)
  !!!-- C H A R A C T E R I S T I C --!!!
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            NULLIFY(solver)
            NULLIFY(solverEquations)
  !!!-- N A V I E R   S T O K E S --!!!
            CALL SOLVERS_SOLVER_GET(solvers,2,solver,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)

            NULLIFY(simpleLoop2)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop2,solvers,err,error,*999)
            !It has two solvers, one for solvine energy equation for vessels and the other one two solve bioheat equation in tissues
!!!-- E N E R G Y  E Q U A T I O N  V E S S E L S --!!!
            !Get the advection-diffusion solver and create the advection-diffusion solver equations
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
!!!-- B I O H E A T --!!!
            !Get the diffusion solver and create the diffusion solver equations
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,2,solver,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Flow solvers
            NULLIFY(whileLoopConditional)
            NULLIFY(iterativeWhileLoop)
            NULLIFY(timeLoopFlow)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,whileLoopConditional,err,error,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(whileLoopConditional,1,timeLoopFlow,err,error,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(timeLoopFlow,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,solvers,err,error,*999)
  !!!-- C H A R A C T E R I S T I C --!!!
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(solverEquations,err,error,*999)
  !!!-- N A V I E R   S T O K E S --!!!
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,2,solver,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(solverEquations,err,error,*999)

            NULLIFY(simpleLoop2)
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop2,solvers,err,error,*999)
            !It has two solvers, one for solvine energy equation for vessels and the other one two solve bioheat equation in tissues
  !!!-- E N E R G Y  E Q U A T I O N  V E S S E L S --!!!
            !Get the advection-diffusion solver and create the advection-diffusion solver equations
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(solverEquations,err,error,*999)
  !!!-- B I O H E A T --!!!
            !Get the diffusion solver and create the diffusion solver equations
            NULLIFY(solver)
            NULLIFY(solverEquations)
            CALL SOLVERS_SOLVER_GET(solvers,2,solver,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(solver,solverEquations,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(solverEquations,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled Navier-Stokes & diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Elias \*
          !Create the CELLML solver equations
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE) THEN
              NULLIFY(simpleLoopCellML1)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,simpleLoopCellML1,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoopCellML1,solvers,err,error,*999)
            END IF
            ! create cellML equations for controling equations
            NULLIFY(solver)
            NULLIFY(cellMLEquations)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_START(solver,cellMLEquations,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_STATIC,err,error,*999)
              !Set the linearity
              CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_LINEAR,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for cellML equations setup of controling equations of Navier-Stokes diffusion advection-diffusion."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE) THEN
              NULLIFY(simpleLoopCellML1)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,simpleLoopCellML1,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoopCellML1,solvers,err,error,*999)
            END IF
            NULLIFY(solver)
            NULLIFY(cellMLEquations)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(solver,cellMLEquations,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_FINISH(cellMLEquations,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "// &
                & TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for cellML equations setup of controling equations of Navier-Stokes diffusion advection-diffusion."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a CellML setup of controling equations of Navier-Stokes diffusion advection-diffusion."
            CALL FlagError(localError,err,error,*999)
          END SELECT    !Elias */
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a coupled Navier-Stokes & diffusion & advection-diffusion equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   c a s e   d e f a u l t
      !-----------------------------------------------------------------
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a coupled Navier-Stokes source diffusion & advection-diffusion equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_ProblemSetup")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_ProblemSetup",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets up the coupled Navier-Stokes & diffusion & advection-diffusion problem pre-solve.
  SUBROUTINE NavierStokesDiffAdvDiff_PreSolve(controlLoop,solver,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: localError


    ENTERS("NavierStokesDiffAdvDiff_PreSolve",ERR,ERROR,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(controlLoop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(controlLoop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes & diffusion advection-diffusion &
              & problem.",err,error,*999)
          END IF
          SELECT CASE(controlLoop%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
            !Update Tskin, Tcore and controlling equations
            SELECT CASE(solver%solve_type)
            ! This switch takes advantage of solve types. However, we type for bioheat solver and NS solver is the same.
            CASE(SOLVER_CELLML_EVALUATOR_TYPE)
              ! Update Tskin, Tcore and then update control parameters such as perfusion and heat generation
! --- C E L L M L  E v a l u a t o r --- !
              CALL NavierStokesDiffAdvDiff_UpdateControlParameters(controlLoop,solver,err,error,*999)
            CASE(SOLVER_NONLINEAR_TYPE)
! --- C h a r a c t e r i s t i c   S o l v e r --- !
              ! Do extrapolation for the first iteration. Copy parameter sets.
              CALL NavierStokesDiffAdvDiff_CharacteristicPreSolve(controlLoop,solver,err,error,*999)
            CASE(SOLVER_DYNAMIC_TYPE)
              SELECT CASE(controlLoop%loop_type)
              CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
! --- 1 D   N a v i e r - S t o k e s   S o l v e r ---
                !Copy parameter sets.
                CALL NavierStokesDiffAdvDiff_1DNSPreSolve(controlLoop,solver,err,error,*999)
              CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
! --- B i o h e a t  S o l v e r --- !
              CASE DEFAULT
                localError="The control type of "//TRIM(NumberToVString(solver%solvers%control_loop%loop_type,"*",err,error))// &
                  & " is invalid for a coupled Navier-Stokes & diffusion & advection-diffusion proplem type."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The solve type of "//TRIM(NumberToVString(SOLVER%SOLVE_TYPE,"*",err,error))// &
                & " is invalid for a coupled Navier-Stokes & diffusion & diffusion-advection problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            ! CALL NavierStokesDiffAdvDiff_UpdateParameters(controlLoop,solver,err,error,*999) !Update Shivering and T_skin
          CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
            !Do nothing
            ! CALL NavierStokesGeometric_Intersection(controlLoop,err,error,*999) !Elias
            CALL NavierStokesDiffAdvDiff_UpdateParameters(controlLoop,solver,err,error,*999) !Update Shivering and T_skin
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a coupled Navier-Stokes & diffusion & advection-diffusion type of a multi physics problem class."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_PreSolve")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_PreSolve",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_PreSolve
  !
  !================================================================================================================================
  !

  !>Sets up the characteristics pre-solve.
  SUBROUTINE NavierStokesDiffAdvDiff_CharacteristicPreSolve(controlLoop,solver,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: localError
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG) :: solver_matrix_idx,iteration,equationsSetIdx
    TYPE(FIELD_TYPE), POINTER :: dependentField
    REAL(DP) :: timeIncrement,currentTime
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable

    ENTERS("NavierStokesDiffAdvDiff_CharacteristicPreSolve",ERR,ERROR,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(controlLoop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(controlLoop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes & diffusion advection-diffusion &
              & problem.",err,error,*999)
          END IF
          SELECT CASE(SOLVER%SOLVE_TYPE)
          ! This switch takes advantage of the uniqueness of the solver types to do pre-solve operations
          ! for each of solvers in the various possible 1D subloops

          ! --- C h a r a c t e r i s t i c   S o l v e r ---
          CASE(SOLVER_NONLINEAR_TYPE)
            CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,currentTime,timeIncrement,err,error,*999)
            iteration = controlLoop%WHILE_LOOP%ITERATION_NUMBER
            EQUATIONS_SET=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
            dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
            ! Characteristic solver effectively solves for the mass/momentum conserving fluxes at the
            ! *NEXT* timestep by extrapolating current field values and then solving a system of nonlinear
            ! equations: cons mass, continuity of pressure, and the characteristics.
            NULLIFY(fieldVariable)
            CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
            IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_INPUT_DATA1_SET_TYPE)%ptr)) THEN
              CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
               & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
              CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
               & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
            END IF
            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
             & FIELD_INPUT_DATA1_SET_TYPE,1.0_DP,err,error,*999)
            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_RESIDUAL_SET_TYPE, &
             & FIELD_INPUT_DATA2_SET_TYPE,1.0_DP,err,error,*999)

            IF(iteration == 1) THEN
              NULLIFY(fieldVariable)
              CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
              IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
              END IF
              ! Extrapolate new W from Q,A if this is the first timestep (otherwise will be calculated based on Navier-Stokes
              ! values)
              CALL Characteristic_Extrapolate(SOLVER,currentTime,timeIncrement,err,error,*999)
            END IF
          CASE DEFAULT
            localError="The solve type of "//TRIM(NumberToVString(SOLVER%SOLVE_TYPE,"*",err,error))// &
              & " is invalid for a 1D Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_CharacteristicPreSolve")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_CharacteristicPreSolve",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_CharacteristicPreSolve


  !
  !================================================================================================================================
  !

  !>Sets up the 1D Navier-Stokes pre-solve.
  SUBROUTINE NavierStokesDiffAdvDiff_1DNSPreSolve(controlLoop,solver,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: localError
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    INTEGER(INTG) :: solver_matrix_idx,iteration,equationsSetIdx
    TYPE(FIELD_TYPE), POINTER :: dependentField
    REAL(DP) :: timeIncrement,currentTime
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX

    ENTERS("NavierStokesDiffAdvDiff_1DNSPreSolve",ERR,ERROR,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(controlLoop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(controlLoop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes & diffusion advection-diffusion &
              & problem.",err,error,*999)
          END IF
          SELECT CASE(SOLVER%SOLVE_TYPE)
          ! This switch takes advantage of the uniqueness of the solver types to do pre-solve operations
          ! for each of solvers in the various possible 1D subloops

      ! --- 1 D   N a v i e r - S t o k e s   S o l v e r ---
          CASE(SOLVER_DYNAMIC_TYPE)
            IF(SOLVER%global_number==2) THEN
              ! update solver matrix
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  SOLVER_MATRICES=>SOLVER_equations%SOLVER_MATRICES
                  IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                    DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                      SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%ptr
                      IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                        SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                      ELSE
                        CALL FlagError("Solver Matrix is not associated.",err,error,*999)
                      END IF
                    END DO
                  ELSE
                    CALL FlagError("Solver Matrices is not associated.",err,error,*999)
                  END IF
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                    IF(ASSOCIATED(dependentField)) THEN
                      CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
                       & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                      CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE, &
                       & FIELD_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
                    ELSE
                      CALL FlagError("Dependent field is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations is not associated.",err,error,*999)
              END IF
            ELSE
              ! --- A d v e c t i o n   S o l v e r ---
              CALL Advection_PreSolve(solver,err,error,*999)
            END IF
            ! Update boundary conditions
            CALL NavierStokesDiffAdvDiff_PreSolveBoundaryConditions(SOLVER,err,error,*999)
          CASE DEFAULT
            localError="The solve type of "//TRIM(NumberToVString(SOLVER%SOLVE_TYPE,"*",err,error))// &
              & " is invalid for a 1D Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_1DNSPreSolve")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_1DNSPreSolve",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_1DNSPreSolve


  !
  !================================================================================================================================
  !

  !>Update the parameters muscle volume, Tskin, Tcore.
  SUBROUTINE NavierStokesDiffAdvDiff_UpdateControlParameters(controlLoop,solver,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    INTEGER(INTG) :: faceIdx,gaussIdx,MPI_IERROR,elemIdx,faceNumber,ms,globalDof,nodeIdx,derivIdx,nodeNumber,counter, &
      & myComputationalNodeNumber,nodeDomain
    !
    REAL(DP) :: area,Tn,phim,muscleVolume,organType,sourceValue ! Because I do not want to define twe independentField variables (problem with CellML) I define organType as Real instead of integer
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    REAL(DP) :: T_skin,T_core,T_shiv,Qshiv_max,Qshiv
    TYPE(SOLVER_TYPE), POINTER :: solverDiffusion  !<A pointer to the solvers !Elias
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquationsDiffusion
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMappingDiffusion
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetDiffusion
    TYPE(FIELD_TYPE), POINTER :: geometricField,dependentField,independentField,sourceField
    LOGICAL :: dependentGeometry
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP,simpleLoop2
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    REAL(DP) :: currentTime,timeIncrement

    !
    ENTERS("NavierStokesDiffAdvDiff_UpdateControlParameters",ERR,ERROR,*999)


    !Get the equations sets, dependent and geometric fields.

    !Get the bioheat solver
    NULLIFY(CONTROL_LOOP_ROOT)
    CONTROL_LOOP_ROOT=>controlLoop%PROBLEM%CONTROL_LOOP
    NULLIFY(CONTROL_LOOP)
    CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
    NULLIFY(simpleLoop2)
    CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
    NULLIFY(solvers)
    CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop2,solvers,err,error,*999)
    ! NULLIFY(solver)
    NULLIFY(solverDiffusion)
    CALL SOLVERS_SOLVER_GET(solvers,2,solverDiffusion,err,error,*999)

    NULLIFY(solverEquationsDiffusion)
    CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)

    NULLIFY(solverMappingDiffusion)
    CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
    NULLIFY(equationsSetDiffusion)

    CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,1,equationsSetDiffusion,err,error,*999)
    !Get dependent, independent, source and geometric fields
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(sourceField)

    CALL EquationsSet_DependentFieldGet(equationsSetDiffusion,dependentField,err,error,*999)
    CALL EquationsSet_SourceFieldGet(equationsSetDiffusion,sourceField,err,error,*999)
    !Get the independent field where T_skin and other parameters are stored
    CALL EquationsSet_IndependentFieldGet(equationsSetDiffusion,independentField,err,error,*999)
    ! get the geometry
    NULLIFY(geometricField)
    CALL Field_GeometricGeneralFieldGet(dependentField,geometricField,dependentGeometry,err,error,*999)


    IF(.NOT.ASSOCIATED(geometricField)) CALL FlagError("Field is not associated.",err,error,*999)
    IF(.NOT.geometricField%FIELD_FINISHED) THEN
      localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//" has not been finished."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(geometricField%TYPE/=FIELD_GEOMETRIC_TYPE) THEN
      localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//" is not a geometric field."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(.NOT.ASSOCIATED(geometricField%GEOMETRIC_FIELD_PARAMETERS)) THEN
      localError="Geometric parameters are not associated for field number "// &
        & TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    ! Calculate skin temperature by averaging Tt over boundary which is skin.
    NULLIFY(coordinateSystem)
    CALL Field_CoordinateSystemGet(geometricField,coordinateSystem,err,error,*999)
    IF(coordinateSystem%NUMBER_OF_DIMENSIONS==3) THEN !only calculate Skin temperature if the body is in 3D
      !Get basis type for the first component of the mesh defined with this geometric field
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(geometricField,decomposition,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)

      !==================== Calculate average skin temperature ===================
      area=geometricField%GEOMETRIC_FIELD_PARAMETERS%surfaceArea

      CALL field_VariableSurfaceIntegral(dependentField,FIELD_U_VARIABLE_TYPE,1,T_skin,err,error,*999)
      T_skin = T_skin/area
      ! CALL Field_ParameterSetUpdateElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,T_skin,ERR,ERROR,*999)
      CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,T_skin,err,error,*999)
      ! ===============================================================================
      ! Calculate total muscle volume (for each partition first)
      ! CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,currentTime,timeIncrement,err,error,*999)
      ! IF (currentTime<timeIncrement) THEN ! Just calculate this for the first time step
      muscleVolume = 0.0_DP
      DO elemIdx=1,decompositionTopology%ELEMENTS%NUMBER_OF_ELEMENTS
        !CALL Field_ParameterSetGetElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elemIdx,3, &
          ! & organType,err,error,*999)
        CALL Field_ParameterSetGetLocalElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elemIdx,3,organType,err,error,*999)
        IF(organType==1) THEN ! TODO change 1 to muscle type
          muscleVolume=muscleVolume+geometricField%GEOMETRIC_FIELD_PARAMETERS%VOLUMES(elemIdx)
        ENDIF
      ENDDO


      ! Collect and sum the total muscle volumes of each partition to obtain total muscle volume. vol=SUM(vol_partition)
      IF(computationalEnvironment%numberOfComputationalNodes>1) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,muscleVolume, &
        & 1,MPI_REAL8,MPI_SUM,computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
      END IF

      ! Update the muscle total volume in independent field.
      CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & 4,muscleVolume,err,error,*999)

      ! ===================== obtain T_core =======================
      ! TODO Should it be rectal temperature only??

      ! CALL Field_VariableVolumeIntegral(dependentField,FIELD_U_VARIABLE_TYPE,1,T_core,err,error,*999)
      ! T_core=T_core/muscleVolume
      !Temporarily use node 28589 temperature for Tcore.
      !get the node domain.
      T_core=0.0_DP
      myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(err,error)
      open(21,file='./input/bioheat/core',status='unknown')
      read(21,*) nodeNumber
      REWIND(21)


      CALL DECOMPOSITION_NODE_DOMAIN_GET(decomposition,nodeNumber,1,nodeDomain,err,error,*999)
      IF (nodeDomain == myComputationalNodeNumber) THEN
        !get the temperature for the node
        CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
          & nodeNumber,1,T_core,err,error,*999)
        close(21)
      END IF

      IF(computationalEnvironment%numberOfComputationalNodes>1) THEN
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,T_core, &
        & 1,MPI_REAL8,MPI_SUM,computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
      END IF
      CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,T_core,err,error,*999)


      CALL Field_ParameterSetUpdateStart(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

      ! TODO: You need to print some values on terminal to check them. Also, you need to exprot them into a file.


      IF(myComputationalNodeNumber==0)THEN
        !get Tcore for a user element number of 1
        CALL Field_ParameterSetGetLocalElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,2, &
          & T_core,err,error,*999)
        CALL Field_ParameterSetGetLocalElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & 1,1,sourceValue,err,error,*999)
        print*, "Tskin,Tcore,muscleVolume,area=",T_skin,T_core, muscleVolume,area
        open(1,file='data',status='old')
        write(1,*) T_skin,T_core,muscleVolume,sourceValue
      ENDIF

      ! DO elementIdx=1,decompositionElements%NUMBER_OF_ELEMENTS
      !   CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
      !     & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
      !   elementVolume=0.0_DP
      !   DO gaussPointIdx=1,numberOfGaussPoints
      !     CALL Field_InterpolateXi(FIRST_PART_DERIV,gaussPoints(1:3,gaussPointIdx),interpolatedPoint(FIELD_U_VARIABLE_TYPE)%PTR, &
      !       & err,error,*999)
      !     CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE, &
      !       & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
      !     elementVolume=elementVolume+InterpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian*gaussWeights(gaussPointIdx)
      !   ENDDO !gaussPointIdx
      !   field%GEOMETRIC_FIELD_PARAMETERS%volumes(elementIdx)=elementVolume
      ! ENDDO !elementIdx
      !==================== Update T_sweating =====================
      ! IF(T_skin<=33.0_DP) THEN
      !   T_swe=42.084_DP-0.15833_DP*T_skin
      ! ELSE IF(T_skin>33.0_DP) THEN
      !   T_swe=36.85_DP
      ! END IF
      !
      ! IF(T_core>T_swe) THEN
      !   mdot_swe=(45.8_DP+739.4_DP*(T_core-T_swe))/(3.6D6)
      ! END IF
      ! !TODO find P_out and ....
      ! w=0.06_DP+mdot_swe*(1.0_DP-0.06_DP)/0.000193_DP !TODO strange relation
      ! q_swe=w*(P_skin-P_out)/(Rswe_cl+1.0_DP/(f_cl*h_swe)) ! W/m2 TODO I need to change this also be carefull about rhoc

      !========================= q_breathing for lung =======================
      ! q_bre=1.0_DP/V_lung*(0.0014_DP*Qm_glob*(34.0_0-T_out)+0.0173_DP*Qm_glob*(5.87_DP-P_out)) ! W/cm3 !TODO change the units divide by rhoc and find T_out and ....

    ENDIF

    EXITS("NavierStokesDiffAdvDiff_UpdateControlParameters")
    RETURN
  999 ERRORSEXITS("NavierStokesDiffAdvDiff_UpdateControlParameters",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_UpdateControlParameters
  !
  !================================================================================================================================
  !

  !>Update the parameters
  SUBROUTINE NavierStokesDiffAdvDiff_UpdateParameters(controlLoop,solver,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    INTEGER(INTG) :: faceIdx,gaussIdx,MPI_IERROR,elemIdx,faceNumber,ms,globalDof,nodeIdx,derivIdx,nodeNumber,counter,myComputationalNodeNumber
    !
    REAL(DP) :: area,Tn,phim,muscleVolume,organType,sourceValue ! Because I do not want to define twe independentField variables (problem with CellML) I define organType as Real instead of integer
    ! TYPE(BASIS_TYPE), POINTER:: basis
    TYPE(COORDINATE_SYSTEM_TYPE), POINTER :: coordinateSystem
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    REAL(DP) :: T_skin,T_core,T_shiv,Qshiv_max,Qshiv
    TYPE(SOLVER_TYPE), POINTER :: solverDiffusion  !<A pointer to the solvers !Elias
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquationsDiffusion
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMappingDiffusion
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetDiffusion
    TYPE(FIELD_TYPE), POINTER :: geometricField,dependentField,independentField,sourceField
    LOGICAL :: dependentGeometry
    TYPE(VARYING_STRING) :: dummyError,localError
    !
    ENTERS("NavierStokesDiffAdvDiff_UpdateParameters",ERR,ERROR,*999)

    IF(SOLVER%GLOBAL_NUMBER==1) THEN
    !
      !Update skin temperature
      NULLIFY(solverDiffusion)
      NULLIFY(solverEquationsDiffusion)

      NULLIFY(solverMappingDiffusion)

      !Get the equations sets, dependent and geometric fields.
      CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,3,solverDiffusion,err,error,*999)
      CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)

      CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
      NULLIFY(equationsSetDiffusion)

      CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,1,equationsSetDiffusion,err,error,*999)
      NULLIFY(dependentField)
      NULLIFY(independentField)
      NULLIFY(sourceField)
      !Get blood temperature
      CALL EquationsSet_DependentFieldGet(equationsSetDiffusion,dependentField,err,error,*999)
      CALL EquationsSet_SourceFieldGet(equationsSetDiffusion,sourceField,err,error,*999)
      !Get the independent field where T_skin and other parameters are stored
      CALL EquationsSet_IndependentFieldGet(equationsSetDiffusion,independentField,err,error,*999)
      NULLIFY(geometricField)
      CALL Field_GeometricGeneralFieldGet(dependentField,geometricField,dependentGeometry,err,error,*999)


      IF(.NOT.ASSOCIATED(geometricField)) CALL FlagError("Field is not associated.",err,error,*999)
      IF(.NOT.geometricField%FIELD_FINISHED) THEN
        localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//" has not been finished."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(geometricField%TYPE/=FIELD_GEOMETRIC_TYPE) THEN
        localError="Field number "//TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//" is not a geometric field."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      IF(.NOT.ASSOCIATED(geometricField%GEOMETRIC_FIELD_PARAMETERS)) THEN
        localError="Geometric parameters are not associated for field number "// &
          & TRIM(NumberToVString(geometricField%USER_NUMBER,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      NULLIFY(coordinateSystem)
      CALL Field_CoordinateSystemGet(geometricField,coordinateSystem,err,error,*999)
      IF(coordinateSystem%NUMBER_OF_DIMENSIONS==3) THEN !only calculate Skin temperature if the body is in 3D
        !Get basis type for the first component of the mesh defined with this geometric field
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(geometricField,decomposition,err,error,*999)
        NULLIFY(decompositionTopology)
        CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)


        !==================== Calculate average skin temperature ===================
        area=geometricField%GEOMETRIC_FIELD_PARAMETERS%surfaceArea        !Allocate Gauss points

        CALL field_VariableSurfaceIntegral(dependentField,FIELD_U_VARIABLE_TYPE,1,T_skin,err,error,*999)
        T_skin = T_skin/area
        ! CALL Field_ParameterSetUpdateElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,T_skin,ERR,ERROR,*999)
        CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,T_skin,err,error,*999)

        ! ===============================================================================
        ! Calculate total muscle volume (for each partition first)
        muscleVolume = 0.0_DP
        DO elemIdx=1,decompositionTopology%ELEMENTS%NUMBER_OF_ELEMENTS
          !CALL Field_ParameterSetGetElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elemIdx,3, &
            ! & organType,err,error,*999)
          CALL Field_ParameterSetGetLocalElement(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
           & elemIdx,3,organType,err,error,*999)
          IF(organType==1) THEN ! TODO change 1 to muscle type
            muscleVolume=muscleVolume+geometricField%GEOMETRIC_FIELD_PARAMETERS%VOLUMES(elemIdx)
          ENDIF
        ENDDO

        ! Collect and sum the total muscle volumes of each partition to obtain total muscle volume. vol=SUM(vol_partition)
        IF(computationalEnvironment%numberOfComputationalNodes>1) THEN
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,muscleVolume, &
          & 1,MPI_REAL8,MPI_SUM,computationalEnvironment%mpiCommunicator,MPI_IERROR)
          CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        END IF

        ! Update the muscle total volume in independent field.
        CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 4,muscleVolume,err,error,*999)

        ! ===================== obtain T_core =======================
        ! TODO Should it be rectal temperature only??

        CALL Field_VariableVolumeIntegral(dependentField,FIELD_U_VARIABLE_TYPE,1,T_core,err,error,*999)
        T_core=T_core/muscleVolume

        CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,T_core,err,error,*999)

        CALL Field_ParameterSetUpdateStart(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

        myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(err,error)
        IF(myComputationalNodeNumber==0)THEN
          CALL Field_ParameterSetGetLocalElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
           & 1,1,sourceValue,err,error,*999)
          print*, "Tskin,Tcore=",T_skin,T_core
          open(1,file='data',status='old')
          write(1,*) T_skin,T_core,sourceValue
        ENDIF

        ! DO elementIdx=1,decompositionElements%NUMBER_OF_ELEMENTS
        !   CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
        !     & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        !   elementVolume=0.0_DP
        !   DO gaussPointIdx=1,numberOfGaussPoints
        !     CALL Field_InterpolateXi(FIRST_PART_DERIV,gaussPoints(1:3,gaussPointIdx),interpolatedPoint(FIELD_U_VARIABLE_TYPE)%PTR, &
        !       & err,error,*999)
        !     CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE, &
        !       & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
        !     elementVolume=elementVolume+InterpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian*gaussWeights(gaussPointIdx)
        !   ENDDO !gaussPointIdx
        !   field%GEOMETRIC_FIELD_PARAMETERS%volumes(elementIdx)=elementVolume
        ! ENDDO !elementIdx
        !==================== Update T_sweating =====================
        ! IF(T_skin<=33.0_DP) THEN
        !   T_swe=42.084_DP-0.15833_DP*T_skin
        ! ELSE IF(T_skin>33.0_DP) THEN
        !   T_swe=36.85_DP
        ! END IF
        !
        ! IF(T_core>T_swe) THEN
        !   mdot_swe=(45.8_DP+739.4_DP*(T_core-T_swe))/(3.6D6)
        ! END IF
        ! !TODO find P_out and ....
        ! w=0.06_DP+mdot_swe*(1.0_DP-0.06_DP)/0.000193_DP !TODO strange relation
        ! q_swe=w*(P_skin-P_out)/(Rswe_cl+1.0_DP/(f_cl*h_swe)) ! W/m2 TODO I need to change this also be carefull about rhoc

        !========================= q_breathing for lung =======================
        ! q_bre=1.0_DP/V_lung*(0.0014_DP*Qm_glob*(34.0_0-T_out)+0.0173_DP*Qm_glob*(5.87_DP-P_out)) ! W/cm3 !TODO change the units divide by rhoc and find T_out and ....



      ENDIF

    END IF


    EXITS("NavierStokesDiffAdvDiff_UpdateParameters")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_UpdateParameters",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_UpdateParameters
  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Navier-Stokes flow pre solve
  SUBROUTINE NavierStokesDiffAdvDiff_PreSolveBoundaryConditions(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET,SOLID_EQUATIONS_SET,FLUID_EQUATIONS_SET
    TYPE(EQUATIONS_SET_DEPENDENT_TYPE), POINTER :: SOLID_DEPENDENT
    TYPE(EQUATIONS_SET_GEOMETRY_TYPE), POINTER :: FLUID_GEOMETRIC
    TYPE(EquationsType), POINTER :: EQUATIONS,SOLID_EQUATIONS,FLUID_EQUATIONS
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,dependentField,geometricField,materialsField
    TYPE(FIELD_TYPE), POINTER :: independentField,SOLID_dependentField,FLUID_GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentFieldVariable,independentFieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,SOLID_SOLVER_EQUATIONS,FLUID_SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING,SOLID_SOLVER_MAPPING,FLUID_SOLVER_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: Solver2
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeIdx,derivativeIdx,versionIdx,variableIdx,numberOfSourceTimesteps,timeIdx,componentIdx
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE,GLOBAL_DERIV_INDEX,node_idx,variable_type
    INTEGER(INTG) :: variable_idx,local_ny,ANALYTIC_FUNCTION_TYPE,component_idx,deriv_idx,dim_idx,version_idx
    INTEGER(INTG) :: element_idx,en_idx,I,J,K,number_of_nodes_xic(3),search_idx,localDof,globalDof,componentBC,previousNodeNumber
    INTEGER(INTG) :: componentNumberVelocity,numberOfDimensions,numberOfNodes,numberOfGlobalNodes
    INTEGER(INTG) :: dependentVariableType,independentVariableType,dependentDof,independentDof,userNodeNumber,localNodeNumber
    INTEGER(INTG) :: EquationsSetIndex,SolidNodeNumber,FluidNodeNumber,equationsSetIdx
    INTEGER(INTG) :: currentTimeLoopIteration,outputIterationNumber,numberOfFittedNodes,computationalNode
    INTEGER(INTG), ALLOCATABLE :: InletNodes(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,DISPLACEMENT_VALUE,VALUE,XI_COORDINATES(3),timeData,QP,QPP,componentValues(3)
    REAL(DP) :: T_COORDINATES(20,3),MU_PARAM,RHO_PARAM,X(3),FluidGFValue,SolidDFValue,NewLaplaceBoundaryValue,Lref,Tref,Mref
    REAL(DP) :: startTime,stopTime,currentTime,timeIncrement
    REAL(DP), POINTER :: MESH_VELOCITY_VALUES(:), GEOMETRIC_PARAMETERS(:), BOUNDARY_VALUES(:)
    REAL(DP), POINTER :: TANGENTS(:,:),NORMAL(:),TIME,ANALYTIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    REAL(DP), POINTER :: independentParameters(:),dependentParameters(:)
    REAL(DP), ALLOCATABLE :: nodeData(:,:),qSpline(:),qValues(:),tValues(:),BoundaryValues(:),fittedNodes(:)
    LOGICAL :: ghostNode,nodeExists,importDataFromFile,ALENavierStokesEquationsSetFound=.FALSE.
    LOGICAL :: SolidEquationsSetFound=.FALSE.,SolidNodeFound=.FALSE.,FluidEquationsSetFound=.FALSE.,parameterSetCreated
    CHARACTER(70) :: inputFile,tempString

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(ANALYTIC_FIELD)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(materialsField)
    NULLIFY(independentField)
    NULLIFY(ANALYTIC_VARIABLE)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(GEOMETRIC_VARIABLE)
    NULLIFY(MATERIALS_VARIABLE)
    NULLIFY(DOMAIN)
    NULLIFY(DOMAIN_NODES)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(MESH_VELOCITY_VALUES)
    NULLIFY(GEOMETRIC_PARAMETERS)
    NULLIFY(BOUNDARY_VALUES)
    NULLIFY(TANGENTS)
    NULLIFY(NORMAL)
    NULLIFY(TIME)
    NULLIFY(ANALYTIC_PARAMETERS)
    NULLIFY(MATERIALS_PARAMETERS)
    NULLIFY(independentParameters)
    NULLIFY(dependentParameters)

    ENTERS("NavierStokesDiffAdvDiff_PreSolveBoundaryConditions",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,startTime,stopTime,CURRENT_TIME,timeIncrement, &
          & currentTimeLoopIteration,outputIterationNumber,ERR,ERROR,*999)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF

          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))


          CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              !If analytic flow waveform, calculate and update
              SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                  IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                        & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        ! Calculate analytic values
                        CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*999)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                        ! Perform spline interpolation of values from a file
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                        DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
                          dependentVariableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
                          NULLIFY(dependentFieldVariable)
                          CALL Field_VariableGet(dependentField,dependentVariableType,dependentFieldVariable,err,error,*999)
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
                            & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                            IF(ASSOCIATED(dependentFieldVariable)) THEN
                              DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE== &
                                  & FIELD_NODE_BASED_INTERPOLATION) THEN
                                  domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
                                  IF(ASSOCIATED(domain)) THEN
                                    IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                      DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                        ! Create the analytic field values type on the dependent field if it does not exist
                                        CALL FIELD_PARAMETER_SET_CREATED(dependentField,dependentVariableType, &
                                          & FIELD_ANALYTIC_VALUES_SET_TYPE,parameterSetCreated,ERR,ERROR,*999)
                                        IF (.NOT. parameterSetCreated) THEN
                                          CALL FIELD_PARAMETER_SET_CREATE(dependentField,dependentVariableType, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                        END IF
                                        !Loop over the local nodes excluding the ghosts.
                                        DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                          userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                          DO derivativeIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                            DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                & numberOfVersions
                                              dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                & VERSIONS(versionIdx)
                                              ! Update dependent field value if this is a splint BC
                                              BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                & CONDITION_TYPES(dependentDof)
                                              IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                                !Update analytic field if file exists and dependent field if boundary condition set
                                                inputFile = './input/interpolatedData/1D/'
                                                IF(dependentVariableType == FIELD_U_VARIABLE_TYPE) THEN
                                                  inputFile = TRIM(inputFile) // 'U/component'
                                                END IF
                                                WRITE(tempString,"(I1.1)") componentIdx
                                                inputFile = TRIM(inputFile) // tempString(1:1) // '/derivative'
                                                WRITE(tempString,"(I1.1)") derivativeIdx
                                                inputFile = TRIM(inputFile) // tempString(1:1) // '/version'
                                                WRITE(tempString,"(I1.1)") versionIdx
                                                inputFile = TRIM(inputFile) // tempString(1:1) // '/'
                                                WRITE(tempString,"(I4.4)") userNodeNumber
                                                inputFile = TRIM(inputFile) // tempString(1:4) // '.dat'
                                                inputFile = TRIM(inputFile)
                                                INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                                                IF(importDataFromFile) THEN
                                                  !Read fitted data from input file (if exists)
                                                  OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')
                                                  ! Header timeData = numberOfTimesteps
                                                  READ(10,*) timeData
                                                  numberOfSourceTimesteps = INT(timeData)
                                                  ALLOCATE(nodeData(numberOfSourceTimesteps,2))
                                                  ALLOCATE(qValues(numberOfSourceTimesteps))
                                                  ALLOCATE(tValues(numberOfSourceTimesteps))
                                                  ALLOCATE(qSpline(numberOfSourceTimesteps))
                                                  nodeData = 0.0_DP
                                                  ! Read in time and dependent value
                                                  DO timeIdx=1,numberOfSourceTimesteps
                                                    READ(10,*) (nodeData(timeIdx,component_idx), component_idx=1,2)
                                                  END DO
                                                  CLOSE(UNIT=10)
                                                  tValues = nodeData(:,1)
                                                  qValues = nodeData(:,2)
                                                  CALL spline_cubic_set(numberOfSourceTimesteps,tValues,qValues, &
                                                    & 2,0.0_DP,2,0.0_DP,qSpline,err,error,*999)
                                                  CALL spline_cubic_val(numberOfSourceTimesteps,tValues,qValues,qSpline, &
                                                    & CURRENT_TIME,VALUE,QP,QPP,err,error,*999)
                                                  DEALLOCATE(nodeData)
                                                  DEALLOCATE(qSpline)
                                                  DEALLOCATE(qValues)
                                                  DEALLOCATE(tValues)
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                    & dependentVariableType,FIELD_VALUES_SET_TYPE,dependentDof, &
                                                    & VALUE,ERR,ERROR,*999)
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                    & dependentVariableType,FIELD_ANALYTIC_VALUES_SET_TYPE,dependentDof, &
                                                    & VALUE,ERR,ERROR,*999)
                                                END IF
                                              END IF ! check if import data file exists
                                            END DO !versionIdx
                                          END DO !derivativeIdx
                                        END DO !nodeIdx
                                        ! Update distributed field values
                                        CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                          & FIELD_VALUES_SET_TYPE,err,error,*999)
                                        CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                          & FIELD_VALUES_SET_TYPE,err,error,*999)
                                        CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                          & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                        CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                          & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                      ELSE
                                        CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain topology is not associated.",err,error,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Domain is not associated.",err,error,*999)
                                  END IF
                                ELSE
                                  CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                                END IF
                              END DO !componentIdx
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                            END IF
                          END IF
                        END DO !variableIdx
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART)
                        ! Using heart lumped parameter model for input
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                        CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,5, &
                          & Lref,err,error,*999)
                        CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,6, &
                          & Tref,err,error,*999)
                        CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,7, &
                          & Mref,err,error,*999)
                        DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
                          dependentVariableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
                          NULLIFY(dependentFieldVariable)
                          CALL Field_VariableGet(dependentField,dependentVariableType,dependentFieldVariable,err,error,*999)
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,dependentFieldVariable, &
                            & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                            IF(ASSOCIATED(dependentFieldVariable)) THEN
                              DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE== &
                                  & FIELD_NODE_BASED_INTERPOLATION) THEN
                                  domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
                                  IF(ASSOCIATED(domain)) THEN
                                    IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                      DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                        !Loop over the local nodes excluding the ghosts.
                                        DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                          userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                          DO derivativeIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                            DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                               & numberOfVersions
                                              dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                & VERSIONS(versionIdx)
                                              BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                & CONDITION_TYPES(dependentDof)
                                              IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE, &
                                                  & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,1,VALUE, &
                                                  & err,error,*999)
                                                ! Convert Q from ml/s to non-dimensionalised form.
                                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,dependentVariableType, &
                                                  & FIELD_VALUES_SET_TYPE,dependentDof,((Lref**3.0)/Tref)*VALUE,err,error,*999)
                                              END IF
                                            END DO !versionIdx
                                          END DO !derivativeIdx
                                        END DO !nodeIdx
                                        ! Update distributed field values
                                        CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                          & FIELD_VALUES_SET_TYPE,err,error,*999)
                                        CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                          & FIELD_VALUES_SET_TYPE,err,error,*999)
                                      ELSE
                                        CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain topology is not associated.",err,error,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Domain is not associated.",err,error,*999)
                                  END IF
                                ELSE
                                  CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                                END IF
                              END DO !componentIdx
                            ELSE
                              CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                            END IF
                          END IF
                        END DO !variableIdx
                      CASE DEFAULT
                        ! Do nothing (might have another use for analytic equations)
                      END SELECT
                    END IF ! Check for analytic equations
                    ELSE
                      CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",err,error,*999)
                END IF
              END IF ! solver equations associated
              ! Update any multiscale boundary values (coupled 0D or non-reflecting)
              CALL NavierStokes_UpdateMultiscaleBoundary(EQUATIONS_SET,BOUNDARY_CONDITIONS,TIME_INCREMENT,err,error,*999)

          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes equation fluid type of a Multi-physics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokesDiffAdvDiff_PreSolveBoundaryConditions")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_PreSolveBoundaryConditions",err,error)
    EXITS("NavierStokesDiffAdvDiff_PreSolveBoundaryConditions")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_PreSolveBoundaryConditions

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes diffusion advection-diffusion problem post solve.
  SUBROUTINE NavierStokesDiffAdvDiff_PostSolve(controlLoop,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: equationsSetIdx,elemIdx,arteryNode,MODEL,elemNum,nodeIdx,numberOfElementNodes,node
    TYPE(SOLVER_TYPE), POINTER :: solverAdvectionDiffusion, solverDiffusion  !<A pointer to the solvers !Elias
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetDiffusion,equationsSetAdvectionDiffusion
    TYPE(FIELD_TYPE), POINTER :: sourceField,dependentField,materialsField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquationsDiffusion,solverEquationsAdvectionDiffusion
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMappingDiffusion,solverMappingAdvectionDiffusion
    REAL(DP) :: Tb,c_param,Tb_old,sourceValue,Tt,Tt_old
    REAL(DP) :: radius,alpha,Tt_ave,Tb_ave
    REAL(DP) :: elemList(15)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(coupledElementsType) :: coupledElements
    INTEGER(INTG) :: node1,node2,arteryElemIdx,tissueElemIdx,arteryElement
    INTEGER(INTG), PARAMETER :: AVERAGE=1,ELEMENT_BASED=2,NO_COUPLED=3
    INTEGER(INTG) :: iteration

    ENTERS("NavierStokesDiffAdvDiff_PostSolve",ERR,ERROR,*999)


    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(controlLoop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(controlLoop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(controlLoop%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)

            SELECT CASE(SOLVER%SOLVE_TYPE)
            CASE(SOLVER_NONLINEAR_TYPE)
              ! Characteristic solver- copy branch Q,A values to new parameter set
              NULLIFY(dependentField)
              NULLIFY(fieldVariable)
              dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
              CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
              IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
              END IF
              iteration = controlLoop%WHILE_LOOP%ITERATION_NUMBER
              IF(iteration == 1) THEN
                CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                 & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)
              END IF
            CASE(SOLVER_DYNAMIC_TYPE)
              SELECT CASE(controlLoop%loop_type)
              CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
                ! Navier-Stokes solver: do nothing
              CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
! --- B i o h e a t  S o l v e r --- !
                CALL NavierStokesDiffAdvDiff_PostSolveOutputData(controlLoop,SOLVER,ERR,ERROR,*999)
              CASE DEFAULT
                localError="The control type of "//TRIM(NumberToVString(controlLoop%loop_type,"*",err,error))// &
                  & " is invalid for a coupled Navier-Stokes & diffusion & advection-diffusion proplem type."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(SOLVER_CELLML_EVALUATOR_TYPE)
! --- C e l l M L  E v a l u a t o r --- !
            CASE DEFAULT
              localError="The solver type of "//TRIM(NumberToVString(SOLVER%SOLVE_TYPE,"*",err,error))// &
                & " is invalid for a 1D Navier-Stokes problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
            MODEL=NO_COUPLED
            SELECT CASE(MODEL)
            CASE(NO_COUPLED)
              !Do nothing
              !Output results
              CALL NavierStokesDiffAdvDiff_PostSolveOutputData(controlLoop,SOLVER,ERR,ERROR,*999)
            CASE(AVERAGE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN

                !Update source field for energy equation of blood based on updated temperature difference between artery wall and blood temperature, Tb(n+1)-Tw(n) !Elias */
                NULLIFY(solverDiffusion)
                NULLIFY(solverEquationsDiffusion)
                NULLIFY(solverEquationsAdvectionDiffusion)
                NULLIFY(solverMappingDiffusion)
                NULLIFY(solverMappingAdvectionDiffusion)
                !Get the equations sets, dependent and source fields.
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,3,solverDiffusion,err,error,*999)
                CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)
                CALL Solver_SolverEquationsGet(solver,solverEquationsAdvectionDiffusion,err,error,*999)
                CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
                CALL SolverEquations_SolverMappingGet(solverEquationsAdvectionDiffusion, &
                  & solverMappingAdvectionDiffusion,err,error,*999)
                !
                DO equationsSetIdx=1,solverMappingDiffusion%NUMBER_OF_EQUATIONS_SETS
                  NULLIFY(equationsSetDiffusion)
                  NULLIFY(equationsSetAdvectionDiffusion)
                  CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,equationsSetIdx,equationsSetDiffusion,err,error,*999)
                  CALL SolverMapping_EquationsSetGet(solverMappingAdvectionDiffusion,equationsSetIdx, &
                    & equationsSetAdvectionDiffusion,err,error,*999)
                  NULLIFY(dependentField)
                  !Get blood temperature
                  CALL EquationsSet_DependentFieldGet(equationsSetAdvectionDiffusion,dependentField,err,error,*999)
                  CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                    & 95,1,Tb,err,error,*999)
                  CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,1, &
                    & 95,1,Tb_old,err,error,*999)
                  IF(equationsSetDiffusion%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSetDiffusion%specification(2)==EQUATIONS_SET_DIFFUSION_EQUATION_TYPE) THEN
                    IF(equationsSetDiffusion%specification(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                      !Get tissue source field
                      NULLIFY(sourceField)
                      NULLIFY(materialsField)
                      CALL EquationsSet_SourceFieldGet(equationsSetDiffusion,sourceField,err,error,*999)
                      CALL EquationsSet_MaterialsFieldGet(equationsSetDiffusion,materialsField,err,error,*999)
                      CALL Field_ParameterSetGetElement(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,3, &
                        & c_param,err,error,*999)
                      CALL Field_ParameterSetGetElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                        & sourceValue,err,error,*999)
                      !Add change in source value to previous source value.
                      sourceValue= sourceValue+c_param*(Tb-Tb_old)
  !                      CALL FIELD_COMPONENT_VALUES_INITIALISE(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
  !                        & sourceValue,ERR,ERROR,*999)
                      !Update tissue source field
                      DO elemIdx=1,80
                        CALL Field_ParameterSetUpdateElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elemIdx,1, &
                          & sourceValue,err,error,*999)
                      END DO
  !                      !Update U1 variable
  !                      CALL Field_ParameterSetsCopy(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
  !                        & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                    END IF
                  END IF
                END DO !equationsSetIdx !Elias /*

                !Output results
  !                CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                CALL NavierStokesDiffAdvDiff_PostSolveOutputData(controlLoop,SOLVER,ERR,ERROR,*999) !Elias
  !                CALL Advection_PostSolve(solver,err,error,*999)
              ELSE IF(SOLVER%GLOBAL_NUMBER==3) THEN

                !Update the source fields if we have converged !Elias */
                NULLIFY(solverAdvectionDiffusion)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,solverAdvectionDiffusion,err,error,*999)
                NULLIFY(solverEquationsAdvectionDiffusion)
                NULLIFY(solverEquationsDiffusion)
                CALL Solver_SolverEquationsGet(solverAdvectionDiffusion,solverEquationsAdvectionDiffusion,err,error,*999)
                CALL Solver_SolverEquationsGet(solver,solverEquationsDiffusion,err,error,*999)
                NULLIFY(solverMappingDiffusion)
                NULLIFY(solverMappingAdvectionDiffusion)
                CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
                CALL SolverEquations_SolverMappingGet(solverEquationsAdvectionDiffusion, &
                  & solverMappingAdvectionDiffusion,err,error,*999)
                DO equationsSetIdx=1,solverMappingAdvectionDiffusion%NUMBER_OF_EQUATIONS_SETS
                  NULLIFY(equationsSetDiffusion)
                  NULLIFY(equationsSetAdvectionDiffusion)
                  CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,equationsSetIdx,equationsSetDiffusion,err,error,*999)
                  CALL SolverMapping_EquationsSetGet(solverMappingAdvectionDiffusion,equationsSetIdx, &
                    & equationsSetAdvectionDiffusion,err,error,*999)
                  NULLIFY(dependentField)
                  CALL EquationsSet_DependentFieldGet(equationsSetDiffusion,dependentField,err,error,*999)
                  IF(equationsSetAdvectionDiffusion%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSetAdvectionDiffusion%specification(2)==EQUATIONS_SET_ADVECTION_EQUATION_TYPE) THEN
                    IF(equationsSetAdvectionDiffusion%specification(3)==EQUATIONS_SET_ADVECTION_DIFFUSION_SUBTYPE) THEN
                      !Update new source value
                      NULLIFY(sourceField)
                      NULLIFY(materialsField)
                      CALL EquationsSet_SourceFieldGet(equationsSetAdvectionDiffusion,sourceField,err,error,*999)
                      CALL EquationsSet_MaterialsFieldGet(equationsSetAdvectionDiffusion,materialsField,err,error,*999)
                      CALL Field_ParameterSetGetElement(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,10,2, &
                        & c_param,err,error,*999)
                      arteryNode=80
                      DO elemIdx=1,80
                        arteryNode=arteryNode+1
                        CALL Field_ParameterSetGetElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,arteryNode-1,1, &
                          & sourceValue,err,error,*999)
                        CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          & elemIdx,1,Tt,err,error,*999)
                        CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,1, &
                          & elemIdx,1,Tt_old,err,error,*999)
                        sourceValue= sourceValue+c_param*(Tt-Tt_old)
                        CALL Field_ParameterSetUpdateElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & arteryNode-1,1,sourceValue,err,error,*999)
                      END DO
  !                      !Update U1 variable
  !                      CALL Field_ParameterSetsCopy(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
  !                        & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                    END IF
                  END IF
                END DO !equationsSetIdx !Elias /*

                !Output results
  !                CALL DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999) !Elias
                CALL NavierStokesDiffAdvDiff_PostSolveOutputData(controlLoop,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE(ELEMENT_BASED)
              !Update Tissue source field for the elements with artery passing through
              !Find these elements and loop over them using an input file
              ! CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, &
              !  & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,1.0_DP, &
              !  & err,error,*999)
              !DO elementIdx in arteryElements
                !Get element source parameters
                !Get the vessel diameter in the element
                !Get the vessel area in the element
                !Get the element volume
                !Get the temperature difference
              !END DO
              !================= Temporarily I store coupled elements this way =================
              !elemList=(/159515,159517,159516,159520,90930,164237,111939,52299,198091,118535,38793,143279,120911,28174, &
              !  & 105491/)
              coupledElements%numberOfArteryElements=24
              ! elemIdx=1
              ALLOCATE(coupledElements%arteryElementNumbers(coupledElements%numberOfArteryElements))

              coupledElements%arteryElementNumbers(1)%arteryElementNumber=2
              coupledElements%arteryElementNumbers(2)%arteryElementNumber=3
              coupledElements%arteryElementNumbers(3)%arteryElementNumber=4
              coupledElements%arteryElementNumbers(4)%arteryElementNumber=5
              coupledElements%arteryElementNumbers(5)%arteryElementNumber=6
              coupledElements%arteryElementNumbers(6)%arteryElementNumber=7
              coupledElements%arteryElementNumbers(7)%arteryElementNumber=8
              coupledElements%arteryElementNumbers(8)%arteryElementNumber=9
              coupledElements%arteryElementNumbers(9)%arteryElementNumber=10
              coupledElements%arteryElementNumbers(10)%arteryElementNumber=11
              coupledElements%arteryElementNumbers(11)%arteryElementNumber=12
              coupledElements%arteryElementNumbers(12)%arteryElementNumber=13
              coupledElements%arteryElementNumbers(13)%arteryElementNumber=14
              coupledElements%arteryElementNumbers(14)%arteryElementNumber=15
              coupledElements%arteryElementNumbers(15)%arteryElementNumber=16
              coupledElements%arteryElementNumbers(16)%arteryElementNumber=17
              coupledElements%arteryElementNumbers(17)%arteryElementNumber=18
              coupledElements%arteryElementNumbers(18)%arteryElementNumber=19
              coupledElements%arteryElementNumbers(19)%arteryElementNumber=20
              coupledElements%arteryElementNumbers(20)%arteryElementNumber=21
              coupledElements%arteryElementNumbers(21)%arteryElementNumber=22
              coupledElements%arteryElementNumbers(22)%arteryElementNumber=23
              coupledElements%arteryElementNumbers(23)%arteryElementNumber=24
              coupledElements%arteryElementNumbers(24)%arteryElementNumber=25

              coupledElements%arteryElementNumbers(1)%numberOfTissueElements=8
              coupledElements%arteryElementNumbers(2)%numberOfTissueElements=8
              coupledElements%arteryElementNumbers(3)%numberOfTissueElements=6
              coupledElements%arteryElementNumbers(4)%numberOfTissueElements=15
              coupledElements%arteryElementNumbers(5)%numberOfTissueElements=3
              coupledElements%arteryElementNumbers(6)%numberOfTissueElements=8
              coupledElements%arteryElementNumbers(7)%numberOfTissueElements=6
              coupledElements%arteryElementNumbers(8)%numberOfTissueElements=6
              coupledElements%arteryElementNumbers(9)%numberOfTissueElements=8
              coupledElements%arteryElementNumbers(10)%numberOfTissueElements=9
              coupledElements%arteryElementNumbers(11)%numberOfTissueElements=11
              coupledElements%arteryElementNumbers(12)%numberOfTissueElements=7
              coupledElements%arteryElementNumbers(13)%numberOfTissueElements=7
              coupledElements%arteryElementNumbers(14)%numberOfTissueElements=9
              coupledElements%arteryElementNumbers(15)%numberOfTissueElements=3
              coupledElements%arteryElementNumbers(16)%numberOfTissueElements=2
              coupledElements%arteryElementNumbers(17)%numberOfTissueElements=4
              coupledElements%arteryElementNumbers(18)%numberOfTissueElements=4
              coupledElements%arteryElementNumbers(19)%numberOfTissueElements=7
              coupledElements%arteryElementNumbers(20)%numberOfTissueElements=7
              coupledElements%arteryElementNumbers(21)%numberOfTissueElements=4
              coupledElements%arteryElementNumbers(22)%numberOfTissueElements=7
              coupledElements%arteryElementNumbers(23)%numberOfTissueElements=7
              coupledElements%arteryElementNumbers(24)%numberOfTissueElements=5

              DO tissueElemIdx=1,coupledElements%numberOfArteryElements
              ALLOCATE(coupledElements%arteryElementNumbers(tissueElemIdx)%tissueElementNumbers(coupledElements% &
                & arteryElementNumbers(tissueElemIdx)%numberOfTissueElements))
              END DO


              coupledElements%arteryElementNumbers(1)%tissueElementNumbers=(/105451,	40713,	76212,	61454, &
                &	142995,	80342,	78436,	121888/)
              coupledElements%arteryElementNumbers(2)%tissueElementNumbers=(/185218,	204426,	36025,	182714, &
                &	124903,	30781,	76633,	13440/)
              coupledElements%arteryElementNumbers(3)%tissueElementNumbers=(/43829,	38359,	6893,	45834,	79633,	159527/)
              coupledElements%arteryElementNumbers(4)%tissueElementNumbers=(/159515,159517,159516,159520, &
                & 90930,164237,111939,52299,198091,118535,38793,143279,120911,28174, 105491/)
              coupledElements%arteryElementNumbers(5)%tissueElementNumbers=(/130595,135972,37742/)
              coupledElements%arteryElementNumbers(6)%tissueElementNumbers=(/135973,90543,99402,176834,214543,168651,212321,107730/)
              coupledElements%arteryElementNumbers(7)%tissueElementNumbers=(/205562,88201,121339,123420,3190,146481/)
              coupledElements%arteryElementNumbers(8)%tissueElementNumbers=(/146486,92427,40865,170495,69683,16321/)
              coupledElements%arteryElementNumbers(9)%tissueElementNumbers=(/22769,43503,66769,78691,76067,27618,172114,158460/)
              coupledElements%arteryElementNumbers(10)%tissueElementNumbers=(/107886,176950,184144,120630,167660,32173,233071, &
                & 218992,8439/)
              coupledElements%arteryElementNumbers(11)%tissueElementNumbers=(/197253,164566,223200,207751,200070,146570,230474, &
                & 218327,39087,155572,171909/)
              coupledElements%arteryElementNumbers(12)%tissueElementNumbers=(/209281,155242,27089,10600,193086,155226,113906/)
              coupledElements%arteryElementNumbers(13)%tissueElementNumbers=(/89737,92674,233534,81544,172854,14619,223750/)
              coupledElements%arteryElementNumbers(14)%tissueElementNumbers=(/85126,5016,217711,39886,222841,84374,144024, &
                & 129163,6674/)
              coupledElements%arteryElementNumbers(15)%tissueElementNumbers=(/127334,127336,1489/)
              coupledElements%arteryElementNumbers(16)%tissueElementNumbers=(/115283,48508/)
              coupledElements%arteryElementNumbers(17)%tissueElementNumbers=(/91763,38220,114435,204146/)
              coupledElements%arteryElementNumbers(18)%tissueElementNumbers=(/65240,39150,180796,180803/)
              coupledElements%arteryElementNumbers(19)%tissueElementNumbers=(/102487,169995,185606,162436,26230,124705,153332/)
              coupledElements%arteryElementNumbers(20)%tissueElementNumbers=(/122755,52417,57010,39758,184885,61083,118837/)
              coupledElements%arteryElementNumbers(21)%tissueElementNumbers=(/70081,62178,210722,202215/)
              coupledElements%arteryElementNumbers(22)%tissueElementNumbers=(/26369,223548,79043,123687,88734,20746,369/)
              coupledElements%arteryElementNumbers(23)%tissueElementNumbers=(/162853,222664,95158,202132,222658,55675,35076/)
              coupledElements%arteryElementNumbers(24)%tissueElementNumbers=(/20939,136314,162936,32088,28476/)

              coupledElements%arteryElementNumbers(1)%arteryElementRadius=3.00235_DP   !mm
              coupledElements%arteryElementNumbers(2)%arteryElementRadius=2.7313_DP    !mm
              coupledElements%arteryElementNumbers(3)%arteryElementRadius=2.44155_DP   !mm
              coupledElements%arteryElementNumbers(4)%arteryElementRadius=2.0967_DP    !mm
              coupledElements%arteryElementNumbers(5)%arteryElementRadius=1.8592_DP    !mm
              coupledElements%arteryElementNumbers(6)%arteryElementRadius=1.7905_DP    !mm
              coupledElements%arteryElementNumbers(7)%arteryElementRadius=1.7414_DP    !mm
              coupledElements%arteryElementNumbers(8)%arteryElementRadius=1.68845_DP   !mm
              coupledElements%arteryElementNumbers(9)%arteryElementRadius=1.63325_DP   !mm
              coupledElements%arteryElementNumbers(10)%arteryElementRadius=1.57155_DP  !mm
              coupledElements%arteryElementNumbers(11)%arteryElementRadius=1.5016_DP   !mm
              coupledElements%arteryElementNumbers(12)%arteryElementRadius=1.43075_DP  !mm
              coupledElements%arteryElementNumbers(13)%arteryElementRadius=1.3496_DP   !mm
              coupledElements%arteryElementNumbers(14)%arteryElementRadius=1.121755_DP !mm
              coupledElements%arteryElementNumbers(15)%arteryElementRadius=1.86625_DP  !mm
              coupledElements%arteryElementNumbers(16)%arteryElementRadius=1.81345_DP  !mm
              coupledElements%arteryElementNumbers(17)%arteryElementRadius=1.7777_DP   !mm
              coupledElements%arteryElementNumbers(18)%arteryElementRadius=1.72575_DP  !mm
              coupledElements%arteryElementNumbers(19)%arteryElementRadius=1.6521_DP   !mm
              coupledElements%arteryElementNumbers(20)%arteryElementRadius=1.5616_DP   !mm
              coupledElements%arteryElementNumbers(21)%arteryElementRadius=1.47835_DP  !mm
              coupledElements%arteryElementNumbers(22)%arteryElementRadius=1.40775_DP  !mm
              coupledElements%arteryElementNumbers(23)%arteryElementRadius=1.3397_DP   !mm
              coupledElements%arteryElementNumbers(24)%arteryElementRadius=1.23915_DP  !mm

              !=============
              IF(SOLVER%GLOBAL_NUMBER==2) THEN

                !Update source field for energy equation of blood based on updated temperature difference between artery wall and blood temperature, Tb(n+1)-Tw(n) !Elias */
                NULLIFY(solverDiffusion)
                NULLIFY(solverEquationsDiffusion)
                NULLIFY(solverEquationsAdvectionDiffusion)
                NULLIFY(solverMappingDiffusion)
                NULLIFY(solverMappingAdvectionDiffusion)
                !Get the equations sets, dependent and source fields.
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,3,solverDiffusion,err,error,*999)
                CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)
                CALL Solver_SolverEquationsGet(solver,solverEquationsAdvectionDiffusion,err,error,*999)
                CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
                CALL SolverEquations_SolverMappingGet(solverEquationsAdvectionDiffusion, &
                  & solverMappingAdvectionDiffusion,err,error,*999)
                !
                DO equationsSetIdx=1,solverMappingDiffusion%NUMBER_OF_EQUATIONS_SETS
                  NULLIFY(equationsSetDiffusion)
                  NULLIFY(equationsSetAdvectionDiffusion)
                  CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,equationsSetIdx,equationsSetDiffusion,err,error,*999)
                  CALL SolverMapping_EquationsSetGet(solverMappingAdvectionDiffusion,equationsSetIdx, &
                    & equationsSetAdvectionDiffusion,err,error,*999)
                  NULLIFY(dependentField)
                  !Get blood temperature
                  CALL EquationsSet_DependentFieldGet(equationsSetAdvectionDiffusion,dependentField,err,error,*999)


                  IF(equationsSetDiffusion%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSetDiffusion%specification(2)==EQUATIONS_SET_DIFFUSION_EQUATION_TYPE) THEN
                    IF(equationsSetDiffusion%specification(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
                      !Get tissue source field
                      NULLIFY(sourceField)
                      NULLIFY(materialsField)
                      CALL EquationsSet_SourceFieldGet(equationsSetDiffusion,sourceField,err,error,*999)

                      DO arteryElemIdx=1,coupledElements%numberOfArteryElements
                        !Calculate average blood temperature in each  artery element for its tissue elements
                        Tb_ave=0
                        node1=coupledElements%arteryElementNumbers(arteryElemIdx)%arteryElementNumber
                        node2=node1+1 !TODO We have only one exception for bifurcation. So the right way to do it is to get node numbers for each artery element
                        !TODO What if we have more than 2 nodes in the element or for the quadratic interpolation?
                        CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          & node1,1,Tb,err,error,*999)
                        Tb_ave=Tb_ave+Tb
                        CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          & node2,1,Tb,err,error,*999)
                        Tb_ave=Tb_ave+Tb
                        ! CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,1, &
                          ! & 95,1,Tb_old,err,error,*999)
                        Tb_ave=Tb_ave/2
                        ! CALL EquationsSet_MaterialsFieldGet(equationsSetDiffusion,materialsField,err,error,*999)
                        ! CALL Field_ParameterSetGetElement(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,3, &
                          ! & c_param,err,error,*999)
                        ! CALL Field_ParameterSetGetElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                          ! & sourceValue,err,error,*999)
                        !Add change in source value to previous source value.
                        ! sourceValue= sourceValue+c_param*(Tb-Tb_old)
                        sourceValue=700.0e-9/(1085e-9*3768)+0.0004771995_DP*Tb_ave
    !                      CALL FIELD_COMPONENT_VALUES_INITIALISE(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
    !                        & sourceValue,ERR,ERROR,*999)
                        !Update tissue source field
                        DO tissueElemIdx=1,coupledElements%arteryElementNumbers(arteryElemIdx)%numberOfTissueElements
                          elemNum=coupledElements%arteryElementNumbers(arteryElemIdx)%tissueElementNumbers(tissueElemIdx)
                          CALL Field_ParameterSetUpdateElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elemNum,1, &
                            & sourceValue,err,error,*999)
                        END DO
                      END DO
  !                      !Update U1 variable
  !                      CALL Field_ParameterSetsCopy(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
  !                        & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                    END IF
                  END IF
                END DO !equationsSetIdx !Elias /*                NULLIFY(sourceField)


              ELSE IF(SOLVER%GLOBAL_NUMBER==3) THEN
                alpha=0.12814434_DP
                !Update the source fields if we have converged !Elias */
                NULLIFY(solverAdvectionDiffusion)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,solverAdvectionDiffusion,err,error,*999)
                NULLIFY(solverEquationsAdvectionDiffusion)
                NULLIFY(solverEquationsDiffusion)
                CALL Solver_SolverEquationsGet(solverAdvectionDiffusion,solverEquationsAdvectionDiffusion,err,error,*999)
                CALL Solver_SolverEquationsGet(solver,solverEquationsDiffusion,err,error,*999)
                NULLIFY(solverMappingDiffusion)
                NULLIFY(solverMappingAdvectionDiffusion)
                CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
                CALL SolverEquations_SolverMappingGet(solverEquationsAdvectionDiffusion, &
                  & solverMappingAdvectionDiffusion,err,error,*999)
                DO equationsSetIdx=1,solverMappingAdvectionDiffusion%NUMBER_OF_EQUATIONS_SETS
                  NULLIFY(equationsSetDiffusion)
                  NULLIFY(equationsSetAdvectionDiffusion)
                  CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,equationsSetIdx,equationsSetDiffusion,err,error,*999)
                  CALL SolverMapping_EquationsSetGet(solverMappingAdvectionDiffusion,equationsSetIdx, &
                    & equationsSetAdvectionDiffusion,err,error,*999)
                  NULLIFY(dependentField)
                  CALL EquationsSet_DependentFieldGet(equationsSetDiffusion,dependentField,err,error,*999)
                  IF(equationsSetAdvectionDiffusion%specification(1)==EQUATIONS_SET_CLASSICAL_FIELD_CLASS.AND. &
                    & equationsSetAdvectionDiffusion%specification(2)==EQUATIONS_SET_ADVECTION_EQUATION_TYPE) THEN
                    IF(equationsSetAdvectionDiffusion%specification(3)==EQUATIONS_SET_ADVECTION_DIFFUSION_SUBTYPE) THEN
                      !Update new source value
                      NULLIFY(sourceField)
                      NULLIFY(materialsField)
                      CALL EquationsSet_SourceFieldGet(equationsSetAdvectionDiffusion,sourceField,err,error,*999)
                      ! CALL EquationsSet_MaterialsFieldGet(equationsSetAdvectionDiffusion,materialsField,err,error,*999)
                      ! CALL Field_ParameterSetGetElement(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,10,2, &
                        ! & c_param,err,error,*999)
                      ! arteryNode=80
                      ! DO elemIdx=1,80
                      ! arteryNode=arteryNode+1
                      ! CALL Field_ParameterSetGetElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,arteryNode-1,1, &
                        ! & sourceValue,err,error,*999)
                      NULLIFY(fieldVariable)
                      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                      topology=>fieldVariable%COMPONENTS(1)%DOMAIN%TOPOLOGY

                      DO arteryElemIdx=1,coupledElements%numberOfArteryElements
                        !Calculate average tissue temperature along the artery element
                        Tt_ave=0
                        DO tissueElemIdx=1,coupledElements%arteryElementNumbers(arteryElemIdx)%numberOfTissueElements
                          elemNum=coupledElements%arteryElementNumbers(arteryElemIdx)%tissueElementNumbers(tissueElemIdx)
                          !TODO Probably does not work for parallel.
                          numberOfElementNodes=4 !TODO get the number of nodes
                          DO nodeIdx=1,numberOfElementNodes
                            node=topology%elements%elements(elemNum)%element_nodes(nodeIdx)
                            CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
                              & node,1,Tt,err,error,*999)
                            Tt_ave=Tt_ave+Tt
                          END DO
                        END DO
                        Tt_ave=Tt_ave/(numberOfElementNodes*coupledElements%arteryElementNumbers(arteryElemIdx)% &
                          & numberOfTissueElements)
                        ! CALL Field_ParameterSetGetNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,1, &
                          ! & elemIdx,1,Tt_old,err,error,*999)
                        ! sourceValue= sourceValue+c_param*(Tt-Tt_old)
                        radius=coupledElements%arteryElementNumbers(arteryElemIdx)%arteryElementRadius
                        sourceValue=4*alpha/radius**2*Tt_ave
                        arteryElement=coupledElements%arteryElementNumbers(arteryElemIdx)%arteryElementNumber
                        CALL Field_ParameterSetUpdateElement(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & arteryElement,1,sourceValue,err,error,*999)
                        ! END DO
                      END DO
    !                      !Update U1 variable
    !                      CALL Field_ParameterSetsCopy(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_NEXT_VALUES_SET_TYPE, &
    !                        & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                    END IF
                  END IF
                END DO !equationsSetIdx !Elias /*

              END IF
              !Output results
              CALL NavierStokesDiffAdvDiff_PostSolveOutputData(controlLoop,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              CALL FlagError(""//TRIM(NUMBER_TO_VSTRING(MODEL,"*",ERR,ERROR))//" model is not implemented",ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a Navier-Stokes & diffusion & advection-diffusion type of a multi physics problem class."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_PostSolve")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_PostSolve",ERR,ERROR)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_PostSolve

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE NavierStokesDiffAdvDiff_ControlLoopPreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: subloop,subloop2,subloop3,iterativeWhileLoop2,iterativeWhileLoop3,whileLoopConditional
    TYPE(SOLVER_TYPE), POINTER :: navierStokesSolver,navierStokesSolver3D,navierStokesSolver1D,solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,solverMapping2
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,equationsSet2,coupledEquationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: numberOfSolvers,solverIdx,solverIdx2,equationsSetIdx,equationsSetIdx2
    INTEGER(INTG) :: subloopIdx,subloopIdx2,subloopIdx3,iteration3D1D
    REAL(DP) :: absolute3D0DTolerance,relative3D0DTolerance
    LOGICAL :: convergedFlag,conditionalWhileLoopContinue
    character(70) :: label

    ENTERS("NavierStokesDiffAdvDiff_ControlLoopPreLoop",err,error,*999)

    NULLIFY(equationsSet)
    NULLIFY(coupledEquationsSet)
    NULLIFY(dependentField)
    NULLIFY(fieldVariable)
    conditionalWhileLoopContinue = .FALSE.   !True for solving flow
    convergedFlag = .FALSE.
    absolute3D0DTolerance = 0.0_DP
    relative3D0DTolerance = 0.0_DP

    IF(ASSOCIATED(controlLoop)) THEN
      SELECT CASE(controlLoop%PROBLEM%specification(3))

      CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          ! Do nothing
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          IF(controlLoop%CONTROL_LOOP_LEVEL==1) THEN
            ! set the condition for solving flow problem.
            NULLIFY(whileLoopConditional)
            CALL CONTROL_LOOP_SUB_LOOP_GET(controlLoop,2,whileLoopConditional,err,error,*999)
            whileLoopConditional%while_loop%CONTINUE_LOOP=conditionalWhileLoopContinue
          END IF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          ! Do nothing
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NumberToVString(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokesDiffAdvDiff_ControlLoopPreLoop")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_ControlLoopPreLoop",err,error)
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_ControlLoopPreLoop
  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE NavierStokesDiffAdvDiff_ControlLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: subloop,subloop2,subloop3,iterativeWhileLoop2,iterativeWhileLoop3,timeLoopFlow
    TYPE(SOLVER_TYPE), POINTER :: navierStokesSolver,navierStokesSolver3D,navierStokesSolver1D,solverDiffusion
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,solverMapping2
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,equationsSet2,coupledEquationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsRobinType), POINTER :: robinConditions
    INTEGER(INTG) :: numberOfSolvers,solverIdx,solverIdx2,equationsSetIdx,equationsSetIdx2
    INTEGER(INTG) :: subloopIdx,subloopIdx2,subloopIdx3,iteration3D1D
    REAL(DP) :: absolute3D0DTolerance,relative3D0DTolerance
    LOGICAL :: convergedFlag
    character(70) :: label
    REAL(DP) :: currentTime,timeIncrement

    TYPE(SCHEDULE_TYPE) :: schedule(6)
    CHARACTER(len=18) :: header
    INTEGER(INTG) :: phases,phaseIdx,phaseNumber
    REAL(DP) :: stopTime

    ENTERS("NavierStokesDiffAdvDiff_ControlLoopPostLoop",err,error,*999)

    NULLIFY(equationsSet)
    NULLIFY(coupledEquationsSet)
    NULLIFY(dependentField)
    NULLIFY(fieldVariable)
    convergedFlag = .FALSE.
    absolute3D0DTolerance = 0.0_DP
    relative3D0DTolerance = 0.0_DP

    IF(ASSOCIATED(controlLoop)) THEN
      SELECT CASE(controlLoop%PROBLEM%specification(3))

      CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
        CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
        open(22,file='./input/bioheat/schedule',status='unknown')
        read(22,*) header
        read(22,*) phases
        DO phaseIdx=1,phases
          read(22,*) schedule(phaseIdx)
          stopTime=schedule(phaseIdx)%timeIncrement*schedule(phaseIdx)%timeSteps
          IF(currentTime<stopTime .AND. currentTime>=schedule(phaseIdx)%start) &
            & phaseNumber=phaseIdx
        END DO
        REWIND(22)
        close(22)

        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          ! For bioheat equaiton in tissue, we do not need to update all the matrices
          !Only two time steps allowed to update the matrices
          IF(controlLoop%control_loop_level==2 .AND. controlLoop%sub_loop_index==3) THEN
            ! Get the equations vector and matrices
            NULLIFY(solvers)
            CALL CONTROL_LOOP_SOLVERS_GET(controlLoop,solvers,err,error,*999)
            NULLIFY(solverDiffusion)
            CALL SOLVERS_SOLVER_GET(solvers,2,solverDiffusion,err,error,*999)
            NULLIFY(solverEquations)
            NULLIFY(solverMapping)
            NULLIFY(equationsSet)
            NULLIFY(equations)
            NULLIFY(vectorEquations)
            NULLIFY(vectorMatrices)
            CALL Solver_SolverEquationsGet(solverDiffusion,solverEquations,err,error,*999)
            CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
            NULLIFY(boundaryConditions)
            CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
            robinConditions=>boundaryConditions%robinBoundaryConditions
            NULLIFY(solverMatrices)
            NULLIFY(solverMatrix)
            CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
            CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
            IF(currentTime>=2*timeIncrement+schedule(phaseNumber)%start) THEN
              ! Turn off the flags for updating matrices that do not vary with time.
              vectorMatrices%dynamicMatrices%matrices(1)%ptr%updateMatrix=.FALSE.
              vectorMatrices%dynamicMatrices%matrices(2)%ptr%updateMatrix=.FALSE.
              solverMatrix%update_matrix=.FALSE.
              IF(ASSOCIATED(robinConditions)) THEN
                robinConditions%updateMatrix=.FALSE.
              END IF
            ELSE
              vectorMatrices%dynamicMatrices%matrices(1)%ptr%updateMatrix=.TRUE.
              vectorMatrices%dynamicMatrices%matrices(2)%ptr%updateMatrix=.TRUE.
              solverMatrix%update_matrix=.TRUE.
              IF(ASSOCIATED(robinConditions)) THEN
                robinConditions%updateMatrix=.TRUE.
              END IF
            END IF
          END IF
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          IF(controlLoop%CONTROL_LOOP_LEVEL/=1) THEN
            ! inner time loop - export data?
            navierStokesSolver=>controlLoop%SUB_LOOPS(1)%ptr%SOLVERS%SOLVERS(2)%ptr
            CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(navierStokesSolver,err,error,*999)
            CALL NavierStokesDiffAdvDiff_CoupleFlow(controlLoop,navierStokesSolver,err,error,*999)
          END IF
          !If we have different phases, modify paramerters
          IF(controlLoop%CONTROL_LOOP_LEVEL==1) THEN
            !Update time increment
            !controlLoop%TIME_LOOP%TIME_INCREMENT=schedule(phaseNumber)%timeIncrement
            IF(currentTime>=schedule(phaseNumber)%start .AND. currentTime<schedule(phaseNumber)%start+timeIncrement) THEN
              CALL NavierStokesDiffAdvDiff_UpdateBoundary(controlLoop,schedule(phaseNumber),err,error,*999)
            END IF
          END IF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          IF(controlLoop%CONTROL_LOOP_LEVEL>2) THEN
            navierStokesSolver=>controlLoop%SOLVERS%SOLVERS(2)%ptr
            CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,err,error,*999)
          END IF
          IF(controlLoop%CONTROL_LOOP_LEVEL==2) THEN
            !Get flow time loop
            NULLIFY(timeLoopFlow)
            CALL CONTROL_LOOP_SUB_LOOP_GET(controlLoop,1,timeLoopFlow,err,error,*999)
            ! Stop executing flow solver for the next time step
            timeLoopFlow%TIME_LOOP%STOP_TIME=0.0_DP
            timeLoopFlow%TIME_LOOP%NUMBER_OF_ITERATIONS=0
          END IF
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NumberToVString(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokesDiffAdvDiff_ControlLoopPostLoop")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_ControlLoopPostLoop",err,error)
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_ControlLoopPostLoop
  !
  !================================================================================================================================
  !

  !>Sets up the diffuion-diffusion problem post solve output data.
  SUBROUTINE NavierStokesDiffAdvDiff_PostSolveOutputData(controlLoop,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError, METHOD, FILENAME
    INTEGER(INTG) :: energy,bioheat

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    LOGICAL, SAVE :: firstCallArtery = .TRUE.
    LOGICAL, SAVE :: firstCallTissue = .TRUE.
    LOGICAL :: EXPORT_FIELD
    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    REAL(DP) :: start,finish,stopTime
    INTEGER:: cr,c1,c2,n,stopIteration
    CHARACTER(10)::c1f
    REAL::t1,t2,rate
    CHARACTER(4)::EN
    CHARACTER(10)::GRGR,KW,GB,YR
    CHARACTER(7)::BL,YE,WH,STR1,STR2
    CHARACTER(15)::STR3
    CHARACTER(50)::BAR
    ! Elias. Color codes
    GRGR=char(27)//'[1;32;40m'
    KW=char(27)//'[1;30;47m'
    GB=char(27)//'[1;32;44m'
    YR=char(27)//'[1;33;41m'
    BL=char(27)//'[1;34m'
    YE=char(27)//'[1;33m'
    WH=char(27)//'[1;37m'
    EN=char(27)//'[0m'
    BAR='                                       '
    n=0


    ENTERS("NavierStokesDiffAdvDiff_PostSolveOutputData",ERR,ERROR,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(controlLoop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(controlLoop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(controlLoop%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
            & PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(controlLoop,CURRENT_TIME,TIME_INCREMENT,err,error,*999) !Elias */
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    !In case of thermoregulation energy is solver 2 and bioheat is solver 3 but in coupled version
                    ! energy is solver 1 and bioheat is solver2
                    SELECT CASE(controlLoop%PROBLEM%SPECIFICATION(3))
                    CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                      energy=2
                      bioheat=3
                    CASE(PROBLEM_COUPLED_BIOHEAT_NAVIERSTOKES_DIFF_ADV_DIFF_SUBTYPE)
                      energy=1
                      bioheat=2
                      controlLoop=>controlLoop%parent_loop
                    CASE DEFAULT
                      LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%SPECIFICATION(3),"*", &
                        & ERR,ERROR))//" is not valid for a coupled Navier-Stokes & diffusion & advection-diffusion type."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    END SELECT

                    CURRENT_LOOP_ITERATION=controlLoop%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=controlLoop%TIME_LOOP%OUTPUT_NUMBER !FREQUENCY

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(controlLoop%TIME_LOOP%CURRENT_TIME<=controlLoop%TIME_LOOP%STOP_TIME) THEN
                        IF(CURRENT_LOOP_ITERATION<10) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        DEPENDENT_REGION=>EQUATIONS_SET%REGION
                        FILE=OUTPUT_FILE
  !          FILE="TRANSIENT_OUTPUT"
!!!!!!!!ADAPT THIS TO WORK WITH DIFFUSION AND NOT JUST FLUID MECHANICS
                         METHOD="FORTRAN"

                        IF(SOLVER%GLOBAL_NUMBER==energy) THEN
                          FILENAME = "./outputArtery/"//"MainTime_"//TRIM(NUMBER_TO_VSTRING(CURRENT_LOOP_ITERATION,"*",ERR,ERROR))
                        ELSE IF(SOLVER%GLOBAL_NUMBER==bioheat) THEN
                          FILENAME = "./outputDiffusion/"//"MainTime_"//TRIM(NUMBER_TO_VSTRING(CURRENT_LOOP_ITERATION,"*",ERR,ERROR))
                        END IF
                         EXPORT_FIELD=.TRUE.
                         IF(EXPORT_FIELD) THEN
                           IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
!                             CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
!                               & err,error,*999)
                             CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                             IF(SOLVER%GLOBAL_NUMBER==energy)THEN
                               IF(firstCallArtery)THEN

                                 CALL FIELD_IO_ELEMENTS_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                                 firstCALLArtery=.FALSE.
                               END IF
                             ELSE IF(SOLVER%GLOBAL_NUMBER==bioheat)THEN
                               IF(firstCallTissue)THEN

                                 CALL FIELD_IO_ELEMENTS_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                                 firstCallTissue=.FALSE.
                               END IF
                             END IF
                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)

                             ! Output iteration information. Elias
                             CALL cpu_time(start)
                             WRITE(STR1,'(f4.2)') TIME_INCREMENT
                             WRITE(STR2,'(f6.1)') start
                             stopTime=controlLoop%TIME_LOOP%STOP_TIME
                             stopIteration=(stopTime/TIME_INCREMENT)
                             n=REAL(CURRENT_LOOP_ITERATION)/stopIteration*50
                             WRITE(STR3,*) CURRENT_LOOP_ITERATION
                             CALL system_clock(count_rate=cr)
                             rate = REAL(cr)
                             CALL system_clock (c1)
                             WRITE(c1f,'(f10.2)')c1/rate

                             WRITE(*,'(a41,a29,a33,a38)')GRGR//"   time step(sec)   "//EN,&
                             & GB//"   cpu_time   "//EN,KW//"   # time steps   "//EN,&
                             & YR//"   Elapsed Time(sec)   "//EN
                             WRITE(*,'(a32,a31,a29,a34)')GRGR//STR1//EN,&
                             & BL//STR2//EN,YE//STR3//EN,&
                             & WH//c1f//EN
                             print*,' '
                             WRITE(*,'(a20,a60,a13)')char(27)//'[1;;41m'//' '//EN,char(27)//'[;47m'//BAR(1:n)//EN//BAR(n+1:50),&
                               & char(27)//'[1;;41m'//' '//EN

                           ENDIF
                         ENDIF

!                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
!                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1 .OR. &
!                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
!                           & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1) THEN
!                            CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,ERR,ERROR,*999)
!                          ENDIF
!                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF !Elias /*
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a coupled Navier-Stokes & diffusion & advection-diffusion type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff_PostSolveOutputData")
    RETURN
999 ERRORS("NavierStokesDiffAdvDiff_PostSolveOutputData",ERR,ERROR)
    EXITS("NavierStokesDiffAdvDiff_PostSolveOutputData")
    RETURN 1

  END SUBROUTINE NavierStokesDiffAdvDiff_PostSolveOutputData

  !
  !================================================================================================================================
  !
  !>This subroutine couples flow and energy for 1D Navier Stokes:
  SUBROUTINE NavierStokesDiffAdvDiff_CoupleFlow(controlLoop,navierStokesSolver, &
    & err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: navierStokesSolver  !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err              !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error     !<The error string

    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoopRoot,CONTROL_LOOP,simpleLoop2
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: solverEnergy
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentFieldNavierStokes,independentFieldEnergy
    REAL(DP) :: currentTime,timeIncrement,storingTimeStep
    INTEGER(INTG) :: compIdx,compIdx2,timeStep

    ENTERS("NavierStokesDiffAdvDiff_CoupleFlow",err,error,*999)

    !Get the dependent field for navier stokes solver
    NULLIFY(solverEquations)
    NULLIFY(solverMapping)
    NULLIFY(equationsSet)
    NULLIFY(dependentFieldNavierStokes)
    CALL Solver_SolverEquationsGet(navierStokesSolver,solverEquations,err,error,*999)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentFieldNavierStokes,err,error,*999)

    !Get the energy solver
    IF(controlLoop%CONTROL_LOOP_LEVEL>2) THEN
      NULLIFY(controlLoopRoot)
      controlLoopRoot=>controlLoop%problem%control_loop
      NULLIFY(CONTROL_LOOP)
      CALL CONTROL_LOOP_GET(controlLoopRoot,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
      !Simple loop for solving vessels energy equation and bioheat equation in tissues
      NULLIFY(simpleLoop2)
      CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
      NULLIFY(solvers)
      CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop2,solvers,err,error,*999)
      NULLIFY(solverEnergy)
      CALL SOLVERS_SOLVER_GET(solvers,1,solverEnergy,err,error,*999)
      ! Get the current time. The unit could be ms or s
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    ELSE
      CALL FlagError("control loop level is not greater than two.",err,error,*999)
    END IF

    ! Get the energy solver independent field
    NULLIFY(solverEquations)
    NULLIFY(solverMapping)
    NULLIFY(equationsSet)
    NULLIFY(independentFieldEnergy)
    CALL Solver_SolverEquationsGet(solverEnergy,solverEquations,err,error,*999)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentFieldEnergy,err,error,*999)

    ! Copy the current time value parameters set from flow dependent field to energy independent field
    currentTime=currentTime/1000.0_DP ! convert to s
    timeIncrement=timeIncrement/1000.0_DP
    storingTimeStep=0.1_DP
    timeStep=NINT(currentTime/storingTimeStep)
    IF(ABS(currentTime-timeStep*storingTimeStep)<(timeIncrement/2.0_DP) .AND. currentTime<1.09) THEN
      Do compIdx=1,dependentFieldNavierStokes%variables(1)%NUMBER_OF_COMPONENTS
        compIdx2=2*timeStep+compIdx
        CALL Field_ParametersToFieldParametersCopy(dependentFieldNavierStokes, &
          & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,compIdx,independentFieldEnergy, &
          & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,compIdx2,err,error,*999)
      END DO
    END IF

    EXITS("NavierStokesDiffAdvDiff_CoupleFlow")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_CoupleFlow",err,error)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_CoupleFlow

  !
  !================================================================================================================================
  !
  !Elias
  !>Calculates the intersection between a given tetrahedron and a segment.
  SUBROUTINE NavierStokesGeometric_Intersection(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: x,y,z
    REAL(DP), ALLOCATABLE :: verticesCoordinates(:)
    INTEGER(INTG) :: numberOfVertices,nodeIdx,node,nic
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: solverDiffusion  !<A pointer to the solvers !Elias
    INTEGER(INTG) :: equationsSetIdx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetDiffusion
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquationsDiffusion
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMappingDiffusion,solverMappingAdvectionDiffusion
    ! REAL(DP) :: elemList(15)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: topology
    TYPE(VerticesType) :: vertices
    REAL(DP) :: xP(3),xQ(3),point(3),point1(3),point2(3)
    TYPE(CrossedType) :: crossedElements
    REAL(DP) :: xLp(4),xLq(4),xLpoint(4)
    REAL(DP) :: tParam
    INTEGER :: counter,splineIdx,acceptedAdjacentElement,foundElementIdx
    INTEGER(INTG) :: adjacentElemIdx(4),adjacentElement(4)
    LOGICAL :: found



    ENTERS("NavierStokesGeometric_Intersection",err,error,*999)

    NULLIFY(solverDiffusion)
    crossedElements%numberOfCrossedElements=0
    ALLOCATE(crossedElements%crossedElements(crossedElements%numberOfCrossedElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate crossed elements array.",err,error,*999)

    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    CALL SOLVERS_SOLVER_GET(solvers,3,solverDiffusion,err,error,*999)


    NULLIFY(solverEquationsDiffusion)
    CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)
    NULLIFY(solverMappingDiffusion)
    CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)



    DO equationsSetIdx=1,solverMappingDiffusion%NUMBER_OF_EQUATIONS_SETS
      NULLIFY(equationsSetDiffusion)
      CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,equationsSetIdx,equationsSetDiffusion,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSetDiffusion,geometricField,err,error,*999)

      NULLIFY(fieldVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
      topology=>fieldVariable%COMPONENTS(1)%DOMAIN%TOPOLOGY

      numberOfVertices=4
      vertices%numberOfVertices=numberOfVertices
      ALLOCATE(vertices%vertices(numberOfVertices))


      ! ALLOCATE(vertices%vertices(1:4)%verticesCoordinates(3))
      DO nodeIdx=1,numberOfVertices
        node=topology%elements%elements(6674)%element_nodes(nodeIdx)
        CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
          & node,1,x,err,error,*999)
        CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
          & node,2,y,err,error,*999)
        CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
          & node,3,z,err,error,*999)
        vertices%vertices(nodeIdx)%verticesCoordinates=(/x,y,z/)
      END DO
    END DO

    xP=(/278.41,-119.14,792.14/)
    xQ=(/277.919,-120.915,797.265/)



    CALL NavierStokesGeometric_BarycentricCoordinates(xP,vertices,xLp,err,error,*999)
    CALL NavierStokesGeometric_BarycentricCoordinates(xQ,vertices,xLq,err,error,*999)

    counter=0
    adjacentElemIdx=(/0,0,0,0/)
    DO nic=1,4
      IF (xLq(nic)<0) THEN
        tParam=-xLp(nic)/(xLq(nic)-xLp(nic))
        IF (tParam>0 .AND. tParam<1) THEN !TODO if the curve passes close enough to an edge of the element then we have 2 tParam that passes the conditions but they both give us almost the same point so no problem but the artery will pass two tissue elements from here.
          point(1)=xP(1)+tParam*(xQ(1)-xP(1))
          point(2)=xP(2)+tParam*(xQ(2)-xP(2))
          point(3)=xP(3)+tParam*(xQ(3)-xP(3))
          counter=counter+1
          adjacentElemIdx(counter)=nic
        END IF
      END IF
    END DO

    CALL NavierStokesIntersection_AddCrossedElement(crossedElements,6674,xP,point,err,error,*999)
    adjacentElement=(/0,0,0,0/)
    adjacentElement(1)=topology%domain%decomposition%topology%elements%elements(6674)%adjacent_elements(adjacentElemIdx(1))% &
      & adjacent_elements(1)

    xP=point

    DO nodeIdx=1,numberOfVertices
      node=topology%elements%elements(adjacentElement(1))%element_nodes(nodeIdx)
      CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
        & node,1,x,err,error,*999)
      CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
        & node,2,y,err,error,*999)
      CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
        & node,3,z,err,error,*999)
      vertices%vertices(nodeIdx)%verticesCoordinates=(/x,y,z/)
    END DO


    CALL NavierStokesGeometric_BarycentricCoordinates(xP,vertices,xLp,err,error,*999)
    CALL NavierStokesGeometric_BarycentricCoordinates(xQ,vertices,xLq,err,error,*999)

    counter=0
    DO nic=1,4
      IF (xLq(nic)<0) THEN
        tParam=-xLp(nic)/(xLq(nic)-xLp(nic))
        IF (tParam>0 .AND. tParam<1) THEN
          point(1)=xP(1)+tParam*(xQ(1)-xP(1))
          point(2)=xP(2)+tParam*(xQ(2)-xP(2))
          point(3)=xP(3)+tParam*(xQ(3)-xP(3))
          counter=counter+1
          adjacentElemIdx(1)=nic
        END IF
      END IF
    END DO

    CALL NavierStokesIntersection_AddCrossedElement(crossedElements,adjacentElement(1),xP,point,err,error,*999)

    ! adjacentElemIdx=2 !MINLOC(ABS(xLpoint), DIM=1,MASK=(xLpoint > 0))

    ! adjacentElement(1)=topology%domain%decomposition%topology%elements%elements(6647)%adjacent_elements(2)% &
    !   & adjacent_elements(1)
    adjacentElement(1)=topology%domain%decomposition%topology%elements%elements(adjacentElement(1))% &
      & adjacent_elements(2)%adjacent_elements(1)

    counter=1
    xP=point
    found=.FALSE.
    OPEN(unit=18, file='arterySpline.csv' , status='old', &
      &  access ='sequential',form='formatted')!,recl=71781*10)
    READ(18,*)
    splineLoop: DO splineIdx=1,51
      READ(18,*) x,y,z
      IF (z<799.975) CYCLE splineLoop

      xQ=(/x,y,z/)


      DO foundElementIdx=1,counter
        !Find vertices of the tetrahedron.
        IF (.NOT. found) THEN
          DO nodeIdx=1,numberOfVertices
            node=topology%elements%elements(adjacentElement(foundElementIdx))%element_nodes(nodeIdx)
            CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
              & node,1,x,err,error,*999)
            CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
              & node,2,y,err,error,*999)
            CALL Field_ParameterSetGetNode(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1, &
              & node,3,z,err,error,*999)
            vertices%vertices(nodeIdx)%verticesCoordinates=(/x,y,z/)
            acceptedAdjacentElement=adjacentElement(foundElementIdx)
          END DO
          CALL NavierStokesGeometric_BarycentricCoordinates(xP,vertices,xLp,err,error,*999)
        END IF


        CALL NavierStokesGeometric_BarycentricCoordinates(xQ,vertices,xLq,err,error,*999)

        IF (ALL(xLq>0 .AND. xLq<1)) THEN
          found=.TRUE.
          CYCLE splineLoop
        END IF

        counter=0
        adjacentElemIdx=(/0,0,0,0/)
        DO nic=1,4
          IF (xLq(nic)<0) THEN
            tParam=-xLp(nic)/(xLq(nic)-xLp(nic))
            IF (tParam>0 .AND. tParam<1) THEN
              point(1)=xP(1)+tParam*(xQ(1)-xP(1))
              point(2)=xP(2)+tParam*(xQ(2)-xP(2))
              point(3)=xP(3)+tParam*(xQ(3)-xP(3))
              counter=counter+1
              adjacentElemIdx(counter)=nic
              IF (counter==1)  THEN
                point1=point
              ELSE IF (counter==2) THEN
                point2=point
              END IF
            ! tParam=tParam+0.1
            END IF
          END IF
        END DO
        IF (counter/=0)  THEN
          point=point1
          EXIT
        ELSE IF (counter==0) THEN
          xP=point2
        END IF

      END DO


      CALL NavierStokesIntersection_AddCrossedElement(crossedElements,acceptedAdjacentElement,xP,point,err,error,*999)
      DO foundElementIdx=1,counter
        adjacentElement(foundElementIdx)=topology%domain%decomposition%topology%elements%elements(acceptedAdjacentElement)% &
          & adjacent_elements(adjacentElemIdx(foundElementIdx))%adjacent_elements(1)
      END DO
      xP=point1
      found=.FALSE.
      ! write(*,*) point
    END DO splineLoop



    EXITS("NavierStokesGeometric_Intersection")
    RETURN
  999 ERRORS("NavierStokesGeometric_Intersection",err,error)
    EXITS("NavierStokesGeometric_Intersection")
    RETURN 1

  END SUBROUTINE NavierStokesGeometric_Intersection

  !
  !================================================================================================================================
  !
  !Elias
  !>Calculates the barycentric ccordinated for a point with respect to a given tetrahedron. See https://en.wikipedia.org/wiki/Barycentric_coordinate_system
  SUBROUTINE NavierStokesGeometric_BarycentricCoordinates(point,vertices,xL,err,error,*)

    !Argument variables
    REAL(DP), INTENT(IN) :: point(3) !<The point coordinates
    TYPE(VerticesType), INTENT(IN) :: vertices !<A pointer to the vertices of the tetrahedron.
    REAL(DP), INTENT(OUT) :: xL(4) !<The barrycentric coordinates for the point.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: T(3,3),A(3,3)
    REAL(DP) :: b(3),c(3)
    REAL(DP) :: det

    ENTERS("NavierStokesGeometric_BarycentricCoordinates",err,error,*999)

    !Find the linear transformation matrix
    T(:,1)=vertices%vertices(1)%verticesCoordinates(:)-vertices%vertices(4)%verticesCoordinates(:)
    T(:,2)=vertices%vertices(2)%verticesCoordinates(:)-vertices%vertices(4)%verticesCoordinates(:)
    T(:,3)=vertices%vertices(3)%verticesCoordinates(:)-vertices%vertices(4)%verticesCoordinates(:)

    CALL INVERT(T,A,det,ERR,ERROR,*999)
    !T.L is invariant and we have L=inv(T)*(r-r4) where r is the point, r4 is the forth vertex, L is the barycentric coordinates.
    b=point-vertices%vertices(4)%verticesCoordinates(:) !r-r4
    CALL MatrixVectorProduct(A,b,c,err,error,*999)

    xL(1:3)=c
    xL(4)=1-xL(1)-xL(2)-xL(3) !Sum of L coordinates is 1. L represents volume to total volume ratio which is a ratio of a segment from point to vertex  to the face.


    EXITS("NavierStokesGeometric_BarycentricCoordinates")
    RETURN
999 ERRORS("NavierStokesGeometric_BarycentricCoordinates",err,error)
    EXITS("NavierStokesGeometric_BarycentricCoordinates")
    RETURN 1

  END SUBROUTINE NavierStokesGeometric_BarycentricCoordinates


  !
  !================================================================================================================================
  !
  !Elias
  !>Adds a new crossed elements to the list.
  SUBROUTINE NavierStokesIntersection_AddCrossedElement(crossedElements,elementNumber,point1,point2,err,error,*)

    !Argument variables
    TYPE(CrossedType), INTENT(INOUT) :: crossedElements !<The elements that are crossed by a given curve
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number
    REAL(DP), INTENT(IN) :: point1(3) !The coordinates of the first intersection
    REAL(DP), INTENT(IN) :: point2(3) !The coordinates of the second intersection
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elemIdx
    TYPE(CrossedElementsType), ALLOCATABLE :: newCrossedElements(:)


    ENTERS("NavierStokesIntersection_AddCrossedElement",err,error,*999)


    ALLOCATE(newCrossedElements(crossedElements%numberOfCrossedElements+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new crossed elements array.",err,error,*999)
    IF(ALLOCATED(crossedElements%crossedElements)) THEN
      DO elemIdx=1,crossedElements%numberOfCrossedElements
        newCrossedElements(elemIdx)=crossedElements%crossedElements(elemIdx)
      ENDDO
    ENDIF

    newCrossedElements(crossedElements%numberOFCrossedElements+1)%elementNumber=elementNumber
    newCrossedElements(crossedElements%numberOFCrossedElements+1)%intersections(1)%intersection=point1
    newCrossedElements(crossedElements%numberOFCrossedElements+1)%intersections(2)%intersection=point2
    CALL MOVE_ALLOC(newCrossedElements,crossedElements%crossedElements)
    crossedElements%numberOFCrossedElements=crossedElements%numberOFCrossedElements+1


    EXITS("NavierStokesIntersection_AddCrossedElement")
    RETURN
999 ERRORS("NavierStokesIntersection_AddCrossedElement",err,error)
    EXITS("NavierStokesIntersection_AddCrossedElement")
    RETURN 1

  END SUBROUTINE NavierStokesIntersection_AddCrossedElement


  !
  !================================================================================================================================
  !
  !> Reads input data from a file
!   SUBROUTINE FLUID_MECHANICS_IO_READ_DATA(SOLVER_TYPE,INPUT_VALUES,NUMBER_OF_DIMENSIONS,INPUT_TYPE, &
!     & INPUT_OPTION,TIME_STEP,LENGTH_SCALE,ERR,ERROR,*)
!
!     INTEGER(INTG):: SOLVER_TYPE,I,INPUT_OPTION,CHECK
!     INTEGER(INTG) :: ERR
!     INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
!     TYPE(VARYING_STRING):: ERROR
!     REAL(DP), POINTER :: INPUT_VALUES(:)
!     INTEGER(INTG):: INPUT_TYPE
!     REAL(DP) :: LENGTH_SCALE
!
!     INTEGER(INTG):: ENDI,TIME_STEP
!     CHARACTER(35) :: INPUT_FILE
!     CHARACTER(29) :: UVEL_FILE
!
!     ENTERS("FLUID_MECHANICS_IO_READ_DATA",ERR,ERROR,*999)
!
!     I=SOLVER_TYPE
!     NUMBER_OF_DIMENSIONS=42
!     ENDI=42
!
!
! !     IF(SOLVER_TYPE==1) THEN !LINEAR
!       IF(INPUT_TYPE==1)THEN !POISSON VECTOR SOURCE TEMPORARY
!         ENDI=SIZE(INPUT_VALUES)
!
! ! WRITE(*,*) "TIME_STEP", TIME_STEP
!
!         IF(INPUT_OPTION==1) THEN
!           IF(TIME_STEP<10) THEN
!             WRITE(UVEL_FILE,'("./input/data/VEL_DATA_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(UVEL_FILE,'("./input/data/VEL_DATA_",I0,".dat")') TIME_STEP
!           ENDIF
!         ELSE IF(INPUT_OPTION==2) THEN
!           IF(TIME_STEP<=10) THEN
!             WRITE(UVEL_FILE,'("./input/data/VEL_DATA_0",I0,".dat")') TIME_STEP-1
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(UVEL_FILE,'("./input/data/VEL_DATA_",I0,".dat")') TIME_STEP-1
!           ENDIF
!         ELSE IF(INPUT_OPTION==3) THEN
!           IF(TIME_STEP<10) THEN
!             WRITE(UVEL_FILE,'("./input/data/ORI_DATA_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(UVEL_FILE,'("./input/data/ORI_DATA_",I0,".dat")') TIME_STEP
!           ENDIF
!         ELSE IF(INPUT_OPTION==4) THEN
!           IF(TIME_STEP<10) THEN
!             WRITE(UVEL_FILE,'("./input/data/U_DATA_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(UVEL_FILE,'("./input/data/U_DATA_",I0,".dat")') TIME_STEP
!           ENDIF
!         ELSE IF(INPUT_OPTION==5) THEN
!           IF(TIME_STEP<10) THEN
!             WRITE(UVEL_FILE,'("./input/data/V_DATA_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(UVEL_FILE,'("./input/data/V_DATA_",I0,".dat")') TIME_STEP
!           ENDIF
!         ELSE IF(INPUT_OPTION==6) THEN
!           IF(TIME_STEP<10) THEN
!             WRITE(UVEL_FILE,'("./input/data/W_DATA_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(UVEL_FILE,'("./input/data/W_DATA_",I0,".dat")') TIME_STEP
!           ENDIF
!         ENDIF
!           CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,UVEL_FILE,ERR,ERROR,*999)
!           OPEN(UNIT=42, FILE=UVEL_FILE,STATUS='unknown')
!           READ(42,*) CHECK
!           IF(CHECK/=ENDI) THEN
!             STOP 'Error during data input - probably wrong Lagrangian/Hermite input file!'
!           ENDIF
!           DO I=1,ENDI
!             READ(42,*) INPUT_VALUES(I)
!           ENDDO
!           CLOSE(42)
!
!       ELSE IF(INPUT_TYPE==42) THEN
! ! do nothing for now
!         IF(INPUT_OPTION==0) THEN
!           !do nothing (default)
!         ELSE IF(INPUT_OPTION==1) THEN
!           ENDI=SIZE(INPUT_VALUES)
!           IF(TIME_STEP<=10) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_00",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<1000) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
!           ENDIF
!           OPEN(UNIT=TIME_STEP, FILE=INPUT_FILE,STATUS='unknown')
!           READ(TIME_STEP,*) CHECK
!           IF(CHECK/=ENDI) THEN
!             STOP 'Error during data input - probably wrong Lagrangian/Hermite input file!'
!           ENDIF
!           DO I=1,ENDI
!             READ(TIME_STEP,*) INPUT_VALUES(I)
!           ENDDO
! ! ! TESTETSTEST
!           INPUT_VALUES=INPUT_VALUES/LENGTH_SCALE
!           WRITE(*,*)'1! INPUT_VALUES=INPUT_VALUES/LENGTH_SCALE'
!           CLOSE(TIME_STEP)
!         ELSE IF(INPUT_OPTION==2) THEN ! For Darcy, invoke the length scale (consistent with reading in the geometry data)
!           ENDI=SIZE(INPUT_VALUES)
!           IF(TIME_STEP<=10) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_00",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<1000) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
!           ENDIF
!           OPEN(UNIT=TIME_STEP, FILE=INPUT_FILE,STATUS='unknown')
!           DO I=1,ENDI
!             READ(TIME_STEP,*) INPUT_VALUES(I)
!             INPUT_VALUES(I) = LENGTH_SCALE * INPUT_VALUES(I)
!           ENDDO
!           CLOSE(TIME_STEP)
!         ELSE IF(INPUT_OPTION==3) THEN
!           ENDI=SIZE(INPUT_VALUES)
!           IF(TIME_STEP<=10) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_00",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<100) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_0",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<500) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') TIME_STEP
!           ELSE IF(TIME_STEP<1000) THEN
!             WRITE(INPUT_FILE,'("./input/motion/DISPLACEMENT_",I0,".dat")') 1000-TIME_STEP
!           ENDIF
!           OPEN(UNIT=TIME_STEP, FILE=INPUT_FILE,STATUS='unknown')
!           DO I=1,ENDI
!             READ(TIME_STEP,*) INPUT_VALUES(I)
!             INPUT_VALUES(I) = LENGTH_SCALE * INPUT_VALUES(I)
!           ENDDO
!           CLOSE(TIME_STEP)
!         ELSE
!           STOP 'Error during data input'
!         ENDIF
!       ENDIF
!
!     EXITS("FLUID_MECHANICS_IO_READ_DATA")
!     RETURN
! 999 ERRORSEXITS("FLUID_MECHANICS_IO_READ_DATA",ERR,ERROR)
!     RETURN 1
!
!   END SUBROUTINE FLUID_MECHANICS_IO_READ_DATA

  !
  !================================================================================================================================
  !
   !Elias
 !>Update ambient boundary conditions
  SUBROUTINE NavierStokesDiffAdvDiff_UpdateBoundary(controlLoop,schedule, &
    & err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop.
    TYPE(SCHEDULE_TYPE), INTENT(IN) :: schedule !<stores parameters in different phases of the simulation
    INTEGER(INTG), INTENT(OUT) :: err              !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error     !<The error string

    !Local variables
    TYPE(VARYING_STRING) :: localError

    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP_ROOT,CONTROL_LOOP,simpleLoop2
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_TYPE), POINTER :: solverDiffusion
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquationsDiffusion
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMappingDiffusion
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSetDiffusion
    TYPE(FIELD_TYPE), POINTER :: independentField
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsRobinType), POINTER :: robinConditions
    REAL(DP) :: value

    ENTERS("NavierStokesDiffAdvDiff_UpdateBoundary",err,error,*999)

    !Get the bioheat solver
    NULLIFY(CONTROL_LOOP_ROOT)
    CONTROL_LOOP_ROOT=>controlLoop%PROBLEM%CONTROL_LOOP
    NULLIFY(CONTROL_LOOP)
    CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
    NULLIFY(simpleLoop2)
    CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,3,simpleLoop2,err,error,*999)
    NULLIFY(solvers)
    CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop2,solvers,err,error,*999)
    NULLIFY(solverDiffusion)
    CALL SOLVERS_SOLVER_GET(solvers,2,solverDiffusion,err,error,*999)
    NULLIFY(solverEquationsDiffusion)
    CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)
    NULLIFY(solverMappingDiffusion)
    CALL SolverEquations_SolverMappingGet(solverEquationsDiffusion,solverMappingDiffusion,err,error,*999)
    NULLIFY(equationsSetDiffusion)
    CALL SolverMapping_EquationsSetGet(solverMappingDiffusion,1,equationsSetDiffusion,err,error,*999)
    !Get dependent, independent, source and geometric fields
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSetDiffusion,independentField,err,error,*999)

   !get robin conditions
    NULLIFY(boundaryConditions)
    boundaryConditions=>solverEquationsDiffusion%BOUNDARY_CONDITIONS
    NULLIFY(robinConditions)
    IF(ASSOCIATED(boundaryConditions)) THEN
        robinConditions=>boundaryConditions%robinBoundaryConditions
    ELSE
      CALL FlagError("boundaryConditions is not associated.",err,error,*999)
    END IF
    ! Update robin conditions
    CALL DistributedVector_AllValuesSet(robinConditions%convectionCoeff,schedule%convectionCoeff,err,error,*999)
    CALL DistributedVector_AllValuesSet(robinConditions%heatFlux,schedule%heatFlux,err,error,*999)

    ! Update Tair and Pv
    CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & 5,schedule%Tair,err,error,*999)
    CALL Field_ComponentValuesInitialise(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & 6,schedule%Pv,err,error,*999)

    CALL Field_ParameterSetUpdateStart(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL Field_ParameterSetUpdateFinish(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

    EXITS("NavierStokesDiffAdvDiff_UpdateBoundary")
    RETURN
999 ERRORSEXITS("NavierStokesDiffAdvDiff_UpdateBoundary",err,error)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff_UpdateBoundary
  !
  !================================================================================================================================
  !
  !Elias
  !>Calculates the volume integral of a component of field variable over the entire volume.
  Subroutine NavierStokesDiffAdvDiff__TcoreCalculate(field,variableType,componentNumber,T_core,err,error,*)

    !Argument variables
    TYPE(FIELD_TYPE), POINTER :: field !<A pointer to the field to get the volume integral for
    INTEGER(INTG),  INTENT(IN) :: variableType !<The field variable type of the field variable component to set \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The component number of the field variable to calculate the integral
    REAL(DP), INTENT(OUT) :: T_core !<The value for volume integral of variable component
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(BASIS_TYPE), POINTER:: basis
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoint(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: interpolatedPointMetrics(:)
    INTEGER(INTG) :: gaussIdx,nodeNumber,globalDof,ms,elemIdx,derivIdx,nodeIdx,MPI_IERROR
    TYPE(VARYING_STRING) :: localError
    REAL(DP) :: phim,Tn,value,coreVolume,integralValue
    LOGICAL :: dependentGeometry,core_element
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable

    ENTERS("NavierStokesDiffAdvDiff__TcoreCalculate",err,error,*999)

    IF(ASSOCIATED(field)) THEN
      IF(field%FIELD_FINISHED) THEN
        IF(variableType>=1.AND.variableType<=FIELD_NUMBER_OF_VARIABLE_TYPES) THEN
          fieldVariable=>FIELD%VARIABLE_TYPE_MAP(variableType)%PTR
          IF(ASSOCIATED(fieldVariable)) THEN
            IF(componentNumber>=1.AND.componentNumber<=fieldVariable%NUMBER_OF_COMPONENTS) THEN
              ! get the geometry
              NULLIFY(geometricField)
              CALL Field_GeometricGeneralFieldGet(field,geometricField,dependentGeometry,err,error,*999)
              IF(.NOT.ASSOCIATED(geometricField)) CALL FlagError("Field is not associated.",err,error,*999)
              !Initialise interpolated point metrics
              NULLIFY(interpolationParameters)
              NULLIFY(interpolatedPoint)
              NULLIFY(interpolatedPointMetrics)
              CALL Field_InterpolationParametersInitialise(geometricField,interpolationParameters,err,error,*999)
              CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoint,err,error,*999)
              CALL Field_InterpolatedPointsMetricsInitialise(interpolatedPoint,interpolatedPointMetrics,err,error,*999)

              NULLIFY(decomposition)
              CALL Field_DecompositionGet(geometricField,decomposition,err,error,*999)
              NULLIFY(decompositionTopology)
              CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
              NULLIFY(domain)
              CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainElements)
              CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)

              SELECT CASE(fieldVariable%COMPONENTS(componentNumber)%INTERPOLATION_TYPE)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                localError="Not implemented for the constant interpolation "
                CALL FLAG_error(localError,err,error,*999)
              CASE(FIELD_NODE_BASED_INTERPOLATION)


                coreVolume=0.0_DP
                integralValue=0.0_DP
                !Loop over domain elements
                DO elemIdx=1,decompositionTopology%ELEMENTS%NUMBER_OF_ELEMENTS
                  IF (core_element) THEN
                    coreVolume=coreVolume+geometricField%GEOMETRIC_FIELD_PARAMETERS%VOLUMES(elemIdx)
                    !caclculate integral (T*dV) for core elements
                    basis=>domainElements%elements(elemIdx)%basis
                    IF(.NOT.ASSOCIATED(basis)) THEN
                      CALL FlagError("basis is not associated.",err,error,*999)
                    END IF
                    quadratureScheme=>basis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
                    IF(.NOT.ASSOCIATED(quadratureScheme)) THEN
                      CALL FlagError("Element basis default quadrature scheme is not associated.",err,error,*999)
                    END IF

                    CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elemIdx, &
                      & interpolationParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                    !Loop over gauss points
                    DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                      CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                        & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)

                      CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE, &
                        & interpolatedPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
                      !Loop over the nodes of the element
                      DO nodeIdx=1,basis%NUMBER_OF_NODES
                        nodeNumber=domainElements%elements(elemIdx)%element_nodes(nodeIdx)
                        DO derivIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)

                          globalDof=domain%mappings%nodes%local_to_global_map(nodeNumber)
                          CALL Field_ParameterSetGetNode(field,variableType,FIELD_VALUES_SET_TYPE,1,derivIdx, &
                            & globalDof,componentNumber,value,err,error,*999)

                          ms=basis%ELEMENT_PARAMETER_INDEX(derivIdx,nodeIdx)
                          phim=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussIdx)
                          integralValue=integralValue+value*phim* &
                            & InterpolatedPointMetrics(variableType)%ptr%jacobian*quadratureScheme%GAUSS_WEIGHTS(gaussIdx)
                        END DO !derivIdx
                      END DO !nodeIdx
                    ENDDO !gaussIdx
                  END IF
                ENDDO !elementIdx

                IF(computationalEnvironment%numberOfComputationalNodes>1) THEN
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,coreVolume, &
                  & 1,MPI_REAL8,MPI_SUM,computationalEnvironment%mpiCommunicator,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                END IF

                IF(computationalEnvironment%numberOfComputationalNodes>1) THEN
                  CALL MPI_ALLREDUCE(MPI_IN_PLACE,integralValue, &
                  & 1,MPI_REAL8,MPI_SUM,computationalEnvironment%mpiCommunicator,MPI_IERROR)
                  CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
                END IF

                T_core=integralValue/coreVolume

              CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                localError="Not implemented for the element based interpolation "
                CALL FLAG_error(localError,err,error,*999)
              CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                localError="Not implemented for the grid point based interpolation "
                CALL FLAG_error(localError,err,error,*999)
              CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                localError="Not implemented for the Gauss point based interpolation "
                CALL FLAG_error(localError,err,error,*999)
              CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                localError="Not implemented for the data point based interpolation "
                CALL FLAG_error(localError,err,error,*999)
              CASE DEFAULT
                localError="The field component interpolation type of "//TRIM(NumberToVString(fieldVariable% &
                  & COMPONENTS(componentNumber)%INTERPOLATION_TYPE,"*",err,error))// &
                  & " is invalid for component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
                  & " of variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
                  & " of field number "//TRIM(NumberToVString(FIELD%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              localError="Component number "//TRIM(NumberToVString(componentNumber,"*",err,error))// &
                & " is invalid for variable type "//TRIM(NumberToVString(variableType,"*",err,error))// &
                & " of field number "//TRIM(NumberToVString(FIELD%USER_NUMBER,"*",err,error))//" which has "// &
                & TRIM(NumberToVString(fieldVariable%NUMBER_OF_COMPONENTS,"*",err,error))// &
                & " components."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The specified field variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
              & " has not been defined on field number "//TRIM(NumberToVString(FIELD%USER_NUMBER,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="The specified variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
            & " is invalid. The variable type must be between 1 and  "// &
            & TRIM(NumberToVString(FIELD_NUMBER_OF_VARIABLE_TYPES,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="Field number "//TRIM(NumberToVString(FIELD%USER_NUMBER,"*",err,error))// &
          & " has not been finished."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Field is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokesDiffAdvDiff__TcoreCalculate")
    RETURN
  999 ERRORSEXITS("NavierStokesDiffAdvDiff__TcoreCalculate",err,error)
    RETURN 1
  END SUBROUTINE NavierStokesDiffAdvDiff__TcoreCalculate
  !
  !================================================================================================================================
  !

END MODULE NAVIER_STOKES_DIFFUSION_ADVECTION_DIFFUSION_ROUTINES
