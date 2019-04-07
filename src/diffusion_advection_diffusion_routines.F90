!> \file
!> \authors Andrew Cookson
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
!> Contributor(s): Chris Bradley, Elias Ghadam Soltani
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


MODULE DIFFUSION_ADVECTION_DIFFUSION_ROUTINES

  USE ADVECTION_DIFFUSION_EQUATION_ROUTINES
  USE ADVECTION_EQUATION_ROUTINES
  USE BaseRoutines
  USE BasisRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CONSTANTS
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
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
  USE MatrixVector
  USE MESH_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverMappingAccessRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"


  IMPLICIT NONE

  PUBLIC DiffusionAdvectionDiffusion_EquationsSetSetup
  PUBLIC DiffusionAdvectionDiffusion_EquationsSetSpecSet
  PUBLIC DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet

  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP
  PUBLIC DiffusionAdvectionDiffusion_ProblemSpecificationSet

  PUBLIC DiffusionAdvectionDiffusion_FiniteElementCalculate

  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE
  PUBLIC DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE



CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet",ERR,ERROR,*999)

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

    EXITS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet",ERR,ERROR)
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion & advection-diffusion coupled equation.
  SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables


    ENTERS("DiffusionAdvectionDiffusion_EquationsSetSetup",ERR,ERROR,*999)

    CALL FlagError("Not implemented.",ERR,ERROR,*999)

    EXITS("DiffusionAdvectionDiffusion_EquationsSetSetup")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_EquationsSetSetup",ERR,ERROR)
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSetup")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a coupled diffusion & advection-diffusion equation finite element equations set.
  SUBROUTINE DiffusionAdvectionDiffusion_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("DiffusionAdvectionDiffusion_FiniteElementCalculate",ERR,ERROR,*999)

    CALL FlagError("Not implemented.",ERR,ERROR,*999)

    EXITS("DiffusionAdvectionDiffusion_FiniteElementCalculate")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_FiniteElementCalculate",ERR,ERROR)
    EXITS("DiffusionAdvectionDiffusion_FiniteElementCalculate")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSpecSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DiffusionAdvectionDiffusion_EquationsSetSpecSet",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)

    EXITS("DiffusionAdvectionDiffusion_EquationsSetSpecSet")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_EquationsSetSpecSet",err,error)
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSpecSet")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSpecSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a coupled diffusion & advection-diffusion problem.
  SUBROUTINE DiffusionAdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("DiffusionAdvectionDiffusion_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
          & PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a coupled diffusion & advection-diffusion type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Diffusion advection-diffusion transport problem specification must have 3 entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("DiffusionAdvectionDiffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_ProblemSpecificationSet",err,error)
    EXITS("DiffusionAdvectionDiffusion_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the coupled diffusion-diffusion equations problem.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION, SOLVER_ADVECTION_DIFFUSION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION, SOLVER_EQUATIONS_ADVECTION_DIFFUSION
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP",ERR,ERROR,*999)

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
        CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
          & err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   coupled source diffusion--advection-diffusion
      !--------------------------------------------------------------------
    CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
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
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,ERR,ERROR,*999)
            !
            !Set the first solver to be a linear solver for the advection-diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_ADVECTION_DIFFUSION,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !
            !Set the second solver to be a linear solver for the diffusion problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER_DIFFUSION,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER_DIFFUSION,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            !
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
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION, &
              & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            !
            !Get the diffusion solver and create the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
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
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_ADVECTION_DIFFUSION,SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_ADVECTION_DIFFUSION,ERR,ERROR,*999)
            !
            !Finish the creation of the diffusion solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER_DIFFUSION,SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS_DIFFUSION,ERR,ERROR,*999)
            !

          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a coupled diffusion & advection-diffusion equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a coupled diffusion & advection-diffusion equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      !-----------------------------------------------------------------
      !   c a s e   d e f a u l t
      !-----------------------------------------------------------------
      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a coupled source diffusion & advection-diffusion equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion-diffusion problem pre-solve.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR


    ENTERS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)

              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !copy current value of concentration_one to another variable
                CALL AdvectionDiffusion_PreSolveStoreCurrentSoln(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !Set source term to be updated value of concentration_two
                !CALL ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)

                !Update indpendent data fields
!                 CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Read in vector data... ",ERR,ERROR,*999)
!                 CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_INPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)

              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                !copy current value of concentration_one to another variable
                CALL Diffusion_PreSolveStoreCurrentSolution(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !compute value of constant source term - evaluated from lamdba*(0.5*(c_1^{t+1}+c_1^{t}) - c_2^{t})
                !CALL DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
              ENDIF
            CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              !Do nothing
              CALL Geometric_Intersection(CONTROL_LOOP,err,error,*999) !Elias
            CASE DEFAULT
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
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

    EXITS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion-diffusion problem post solve.
  SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
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
    INTEGER(INTG), PARAMETER :: AVERAGE=1,ELEMENT_BASED=2

    ENTERS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
            IF(SOLVER%GLOBAL_NUMBER==1) THEN
              !Output results
            !  CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
              !Output results
            !  CALL DIFFUSION_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            END IF
          CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
            MODEL=ELEMENT_BASED
            SELECT CASE(MODEL)
            CASE(AVERAGE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN

                !Update source field for energy equation of blood based on updated temperature difference between artery wall and blood temperature, Tb(n+1)-Tw(n) !Elias */
                NULLIFY(solverDiffusion)
                NULLIFY(solverEquationsDiffusion)
                NULLIFY(solverEquationsAdvectionDiffusion)
                NULLIFY(solverMappingDiffusion)
                NULLIFY(solverMappingAdvectionDiffusion)
                !Get the equations sets, dependent and source fields.
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,solverDiffusion,err,error,*999)
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
                CALL DiffusionAdvectionDiffusion_PostSolveOutputData(CONTROL_LOOP,SOLVER,ERR,ERROR,*999) !Elias
  !                CALL Advection_PostSolve(solver,err,error,*999)
              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN

                !Update the source fields if we have converged !Elias */
                NULLIFY(solverAdvectionDiffusion)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,solverAdvectionDiffusion,err,error,*999)
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
                CALL DiffusionAdvectionDiffusion_PostSolveOutputData(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
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
              IF(SOLVER%GLOBAL_NUMBER==1) THEN

                !Update source field for energy equation of blood based on updated temperature difference between artery wall and blood temperature, Tb(n+1)-Tw(n) !Elias */
                NULLIFY(solverDiffusion)
                NULLIFY(solverEquationsDiffusion)
                NULLIFY(solverEquationsAdvectionDiffusion)
                NULLIFY(solverMappingDiffusion)
                NULLIFY(solverMappingAdvectionDiffusion)
                !Get the equations sets, dependent and source fields.
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,solverDiffusion,err,error,*999)
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


              ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
                alpha=0.12814434_DP
                !Update the source fields if we have converged !Elias */
                NULLIFY(solverAdvectionDiffusion)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,solverAdvectionDiffusion,err,error,*999)
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
              CALL DiffusionAdvectionDiffusion_PostSolveOutputData(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE DEFAULT
              CALL FlagError(""//TRIM(NUMBER_TO_VSTRING(MODEL,"*",ERR,ERROR))//" model is not implemented",ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
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

    EXITS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE")
    RETURN
999 ERRORSEXITS("DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE DIFFUSION_ADVECTION_DIFFUSION_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the diffuion-diffusion problem post solve output data.
  SUBROUTINE DiffusionAdvectionDiffusion_PostSolveOutputData(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
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


    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

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


    ENTERS("DiffusionAdvectionDiffusion_PostSolveOutputData",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(control_loop%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(control_loop%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a diffusion-advection diffusion problem.", &
              & err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                !CALL ADVECTION_DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                !CALL DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_THERMOREGULATION_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999) !Elias */
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr

                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER !FREQUENCY

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
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
                        IF(SOLVER%GLOBAL_NUMBER==1) THEN
                          FILENAME = "./outputArtery/"//"MainTime_"//TRIM(NUMBER_TO_VSTRING(CURRENT_LOOP_ITERATION,"*",ERR,ERROR))
                        ELSE IF(SOLVER%GLOBAL_NUMBER==2) THEN
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
!                             CALL FIELD_IO_ELEMENTS_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                             CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)

                             ! Output iteration information. Elias
                             CALL cpu_time(start)
                             WRITE(STR1,'(f4.2)') TIME_INCREMENT
                             WRITE(STR2,'(f6.1)') start
                             stopTime=CONTROL_LOOP%TIME_LOOP%STOP_TIME
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
              LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
                & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
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

    EXITS("DiffusionAdvectionDiffusion_PostSolveOutputData")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_PostSolveOutputData",ERR,ERROR)
    EXITS("DiffusionAdvectionDiffusion_PostSolveOutputData")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_PostSolveOutputData


  !
  !================================================================================================================================
  !
  !Elias
  !>Calculates the intersection between a given tetrahedron and a segment.
  SUBROUTINE Geometric_Intersection(controlLoop,err,error,*)

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



    ENTERS("Geometric_Intersection",err,error,*999)

    NULLIFY(solverDiffusion)
    crossedElements%numberOfCrossedElements=0
    ALLOCATE(crossedElements%crossedElements(crossedElements%numberOfCrossedElements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate crossed elements array.",err,error,*999)

    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    CALL SOLVERS_SOLVER_GET(solvers,2,solverDiffusion,err,error,*999)


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



    CALL Geometric_BarycentricCoordinates(xP,vertices,xLp,err,error,*999)
    CALL Geometric_BarycentricCoordinates(xQ,vertices,xLq,err,error,*999)

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

    CALL Intersection_AddCrossedElement(crossedElements,6674,xP,point,err,error,*999)
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


    CALL Geometric_BarycentricCoordinates(xP,vertices,xLp,err,error,*999)
    CALL Geometric_BarycentricCoordinates(xQ,vertices,xLq,err,error,*999)

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

    CALL Intersection_AddCrossedElement(crossedElements,adjacentElement(1),xP,point,err,error,*999)

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
          CALL Geometric_BarycentricCoordinates(xP,vertices,xLp,err,error,*999)
        END IF


        CALL Geometric_BarycentricCoordinates(xQ,vertices,xLq,err,error,*999)

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


      CALL Intersection_AddCrossedElement(crossedElements,acceptedAdjacentElement,xP,point,err,error,*999)
      DO foundElementIdx=1,counter
        adjacentElement(foundElementIdx)=topology%domain%decomposition%topology%elements%elements(acceptedAdjacentElement)% &
          & adjacent_elements(adjacentElemIdx(foundElementIdx))%adjacent_elements(1)
      END DO
      xP=point1
      found=.FALSE.
      ! write(*,*) point
    END DO splineLoop



    EXITS("Geometric_Intersection")
    RETURN
  999 ERRORS("Geometric_Intersection",err,error)
    EXITS("Geometric_Intersection")
    RETURN 1

  END SUBROUTINE Geometric_Intersection

  !
  !================================================================================================================================
  !
  !Elias
  !>Calculates the barycentric ccordinated for a point with respect to a given tetrahedron. See https://en.wikipedia.org/wiki/Barycentric_coordinate_system
  SUBROUTINE Geometric_BarycentricCoordinates(point,vertices,xL,err,error,*)

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

    ENTERS("Geometric_BarycentricCoordinates",err,error,*999)

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


    EXITS("Geometric_BarycentricCoordinates")
    RETURN
999 ERRORS("Geometric_BarycentricCoordinates",err,error)
    EXITS("Geometric_BarycentricCoordinates")
    RETURN 1

  END SUBROUTINE Geometric_BarycentricCoordinates


  !
  !================================================================================================================================
  !
  !Elias
  !>Adds a new crossed elements to the list.
  SUBROUTINE Intersection_AddCrossedElement(crossedElements,elementNumber,point1,point2,err,error,*)

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


    ENTERS("Intersection_AddCrossedElement",err,error,*999)


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


    EXITS("Intersection_AddCrossedElement")
    RETURN
999 ERRORS("Intersection_AddCrossedElement",err,error)
    EXITS("Intersection_AddCrossedElement")
    RETURN 1

  END SUBROUTINE Intersection_AddCrossedElement


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

  ! OK
  !================================================================================================================================
  !


END MODULE DIFFUSION_ADVECTION_DIFFUSION_ROUTINES
