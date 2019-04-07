  !>Assembles the solver matrices and rhs from the static equations.
  SUBROUTINE SOLVER_MATRICES_STATIC_ASSEMBLE(SOLVER,SELECTION_TYPE,ERR,ERROR,*)

    !Argument variable
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The type of matrix selection \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: dependent_variable_type,interface_variable_type,equations_column_number,equations_matrix_idx, &
      & equations_matrix_number,interface_row_number,equations_row_number,equations_row_number2,equations_set_idx, & 
      & interface_column_number,interface_condition_idx,interface_matrix_idx,LINEAR_VARIABLE_TYPE,rhs_boundary_condition, &
      & rhs_global_dof,equations_matrix_idx2,rhs_variable_dof,rhs_variable_type,variable_boundary_condition,solver_matrix_idx, &
      & solver_row_idx,solver_row_number,variable_dof,variable_global_dof,variable_idx,variable_type,&
      & dirichlet_idx,dirichlet_row,number_of_interface_matrices
    REAL(SP) :: systemElapsed,SYSTEM_TIME1(1),SYSTEM_TIME2(1),userElapsed,USER_TIME1(1),USER_TIME2(1)
    REAL(DP) :: DEPENDENT_VALUE,LINEAR_VALUE,LINEAR_VALUE_SUM,MATRIX_VALUE,RESIDUAL_VALUE,RHS_VALUE,row_coupling_coefficient, &
      & SOURCE_VALUE,VALUE,RHS_INTEGRATED_VALUE
    REAL(DP), POINTER :: RHS_PARAMETERS(:),CHECK_DATA(:),CHECK_DATA2(:),CHECK_DATA3(:),CHECK_DATA4(:)
    LOGICAL :: SUBTRACT_FIXED_BCS_FROM_RESIDUAL,HAS_INTEGRATED_VALUES,R_HAS_INTEGRATED_VALUES !Elias
    TYPE(REAL_DP_PTR_TYPE), ALLOCATABLE :: DEPENDENT_PARAMETERS(:)
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: DEPENDENT_BOUNDARY_CONDITIONS,RHS_BOUNDARY_CONDITIONS
    TYPE(DistributedMatrixType), POINTER :: PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(DistributedVectorType), POINTER :: LAGRANGE_VECTOR,DEPENDENT_VECTOR,DISTRIBUTED_SOURCE_VECTOR,EQUATIONS_RHS_VECTOR, &
      & LINEAR_TEMP_VECTOR,INTERFACE_TEMP_VECTOR,RESIDUAL_VECTOR,SOLVER_RESIDUAL_VECTOR,SOLVER_RHS_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: RHS_DOMAIN_MAPPING,VARIABLE_DOMAIN_MAPPING
    TYPE(EquationsJacobianType), POINTER :: JACOBIAN_MATRIX
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix,LINEAR_MATRIX
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_TO_SOLVER_MAPS_TYPE), POINTER :: EQUATIONS_TO_SOLVER_MAP
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: INTERFACE_VARIABLE,DEPENDENT_VARIABLE,LINEAR_VARIABLE,RHS_VARIABLE
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: INTERFACE_LAGRANGE
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: INTERFACE_RHS_MAPPING
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: INTERFACE_RHS_VECTOR
    TYPE(INTERFACE_TO_SOLVER_MAPS_TYPE), POINTER :: INTERFACE_TO_SOLVER_MAP
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES
  
    ENTERS("SOLVER_MATRICES_STATIC_ASSEMBLE",ERR,ERROR,*999)
  
    IF(ASSOCIATED(SOLVER)) THEN
      SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
      IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
        SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
        IF(ASSOCIATED(SOLVER_MAPPING)) THEN
          SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
          IF(ASSOCIATED(SOLVER_MATRICES)) THEN
            !Assemble the solver matrices
            NULLIFY(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
              !Assemble solver matrices
              IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPUTimer(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              !Loop over the solver matrices
              DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                SOLVER_MATRIX=>SOLVER_MATRICES%matrices(solver_matrix_idx)%ptr
                IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                  IF(SOLVER_MATRIX%UPDATE_MATRIX) THEN
                    SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                    IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN
                      !Initialise matrix to zero
                      CALL DistributedMatrix_AllValuesSet(SOLVER_DISTRIBUTED_MATRIX,0.0_DP,ERR,ERROR,*999)
                      !Loop over the equations sets
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        !First Loop over the linear equations matrices
                        DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                          & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_LINEAR_EQUATIONS_MATRICES
                          EQUATIONS_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                            & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%LINEAR_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & equations_matrix_idx)%ptr
                          IF(ASSOCIATED(EQUATIONS_TO_SOLVER_MAP)) THEN
                            equationsMatrix=>EQUATIONS_TO_SOLVER_MAP%EQUATIONS_MATRIX
                            IF(ASSOCIATED(equationsMatrix)) THEN

                              rhsMapping=>vectorMapping%rhsMapping !Elias /*
                              IF(ASSOCIATED(rhsMapping)) THEN
                                BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS 
                                IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                  RHS_VARIABLE=>rhsMapping%rhsVariable
                                  CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,RHS_VARIABLE, &
                                    & RHS_BOUNDARY_CONDITIONS,ERR,ERROR,*999) 
                                  IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN
                                    !Add K_R to the stiffness matrix=K_e+K_R. 
                                      !Later on this updated matrix is used in solver_rhs_vector as well.
                                    CALL BoundaryConditions_RobinDynamicIntegrate(RHS_BOUNDARY_CONDITIONS, &
                                      & equationsMatrix,ERR,ERROR,*999)   
                                  ELSE
                                    CALL FlagError("RHS boundary conditions variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Equations set boundary conditions is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                              ENDIF !Elias */

                              CALL SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx,1.0_DP,equationsMatrix, &
                                & ERR,ERROR,*999)
                            ELSE
                              CALL FlagError("The equations matrix is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("The equations matrix equations to solver map is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !equations_matrix_idx
                        IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                          & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                          & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
                          !Now set the values from the equations Jacobian
                          DO equations_matrix_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_EQUATIONS_JACOBIANS 
                            JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                              & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAPS( &
                              & equations_matrix_idx)%ptr
                            IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                              JACOBIAN_MATRIX=>JACOBIAN_TO_SOLVER_MAP%JACOBIAN_MATRIX
                              IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                                CALL SOLVER_MATRIX_JACOBIAN_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx,1.0_DP,JACOBIAN_MATRIX, &
                                  & ERR,ERROR,*999)
                              ELSE
                                CALL FlagError("Jacobian matrix is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="Jacobian to solver map is not associated for Jacobian number "// &
                                & TRIM(NumberToVString(equations_matrix_idx,"*",ERR,ERROR))//"."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ENDDO
                        ENDIF
                      ENDDO !equations_set_idx
                      !Loop over any interface conditions
                      DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                        !Loop over the interface matrices
                        DO interface_matrix_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%NUMBER_OF_INTERFACE_MATRICES
                          INTERFACE_TO_SOLVER_MAP=>SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                            & INTERFACE_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%INTERFACE_EQUATIONS_TO_SOLVER_MATRIX_MAPS( &
                            & interface_matrix_idx)%ptr
                          IF(ASSOCIATED(INTERFACE_TO_SOLVER_MAP)) THEN
                            INTERFACE_MATRIX=>INTERFACE_TO_SOLVER_MAP%INTERFACE_MATRIX
                            IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
                              CALL SOLVER_MATRIX_INTERFACE_MATRIX_ADD(SOLVER_MATRIX,interface_condition_idx,[1.0_DP,1.0_DP], &
                                & INTERFACE_MATRIX,ERR,ERROR,*999)
                            ELSE
                              CALL FlagError("The interface matrix is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("The interface matrix interface to solver map is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ENDDO !interface_matrix_idx
                      ENDDO !interface_condition_idx
                      !Update the solver matrix values
                      CALL DistributedMatrix_UpdateStart(SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                      IF(ASSOCIATED(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)) THEN
                        CALL DistributedMatrix_UpdateFinish(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
                      ENDIF
                      PREVIOUS_SOLVER_DISTRIBUTED_MATRIX=>SOLVER_DISTRIBUTED_MATRIX
                    ELSE
                      CALL FlagError("Solver matrix distributed matrix is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDIF !Update matrix
                ELSE
                  CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDDO !solver_matrix_idx
              IF(ASSOCIATED(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)) THEN
                CALL DistributedMatrix_UpdateFinish(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPUTimer(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                userElapsed=USER_TIME2(1)-USER_TIME1(1)
                systemElapsed=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
                  & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
                CALL Profiling_TimingsOutput(1,"Solver matrices assembly",userElapsed,systemElapsed,err,error,*999)
              ENDIF
            ENDIF
            !The solver matrices have only one residual vector
            NULLIFY(SOLVER_RESIDUAL_VECTOR)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
              !Assemble residual vector
              !We assemble residual vector before RHS vector, then when assembling the RHS vector we subtract
              !the RHS terms for fixed BCs from the residual vector as this residual evaluation uses a matrix
              !vector product of the full equations matrix rather than the reduced solver matrix
              IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPUTimer(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER_MATRICES%UPDATE_RESIDUAL) THEN
                SOLVER_RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                  !Initialise the residual to zero
                  CALL DistributedVector_AllValuesSet(SOLVER_RESIDUAL_VECTOR,0.0_DP,ERR,ERROR,*999)
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          NULLIFY(vectorEquations)
                          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                          vectorMatrices=>vectorEquations%vectorMatrices
                          IF(ASSOCIATED(vectorMatrices)) THEN
                            vectorMapping=>vectorEquations%vectorMapping
                            IF(ASSOCIATED(vectorMapping)) THEN
                              !Calculate the contributions from any linear matrices
                              linearMapping=>vectorMapping%linearMapping
                              IF(ASSOCIATED(linearMapping)) THEN
                                linearMatrices=>vectorMatrices%linearMatrices
                                IF(ASSOCIATED(linearMatrices)) THEN
                                  DO equations_matrix_idx=1,linearMatrices%numberOfLinearMatrices
                                    LINEAR_MATRIX=>linearMatrices%matrices(equations_matrix_idx)%ptr
                                    IF(ASSOCIATED(LINEAR_MATRIX)) THEN
                                      LINEAR_VARIABLE_TYPE=linearMapping%equationsMatrixToVarMaps(equations_matrix_idx)% &
                                        & variableType
                                      LINEAR_VARIABLE=>linearMapping%equationsMatrixToVarMaps(equations_matrix_idx)% &
                                        & variable
                                      IF(ASSOCIATED(LINEAR_VARIABLE)) THEN
                                        LINEAR_TEMP_VECTOR=>LINEAR_MATRIX%tempVector
                                        !Initialise the linear temporary vector to zero
                                        CALL DistributedVector_AllValuesSet(LINEAR_TEMP_VECTOR,0.0_DP,ERR,ERROR,*999)
                                        NULLIFY(DEPENDENT_VECTOR)
                                        CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,LINEAR_VARIABLE_TYPE, &
                                          & FIELD_VALUES_SET_TYPE,DEPENDENT_VECTOR,ERR,ERROR,*999)
                                        CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                                          & LINEAR_MATRIX%MATRIX,DEPENDENT_VECTOR,LINEAR_TEMP_VECTOR,ERR,ERROR,*999)
                                      ELSE
                                        CALL FlagError("Linear variable is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      LOCAL_ERROR="Linear matrix is not associated for linear matrix number "// &
                                        & TRIM(NumberToVString(equations_matrix_idx,"*",ERR,ERROR))//"."
                                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                    ENDIF
                                  ENDDO !equations_matrix_idx
                                ELSE
                                  CALL FlagError("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                              !Calculate the solver residual
                              nonlinearMapping=>vectorMapping%nonlinearMapping
                              IF(ASSOCIATED(nonlinearMapping)) THEN
                                nonlinearMatrices=>vectorMatrices%nonlinearMatrices
                                IF(ASSOCIATED(nonlinearMatrices)) THEN
                                  RESIDUAL_VECTOR=>nonlinearMatrices%RESIDUAL
                                  !Loop over the rows in the equations set
                                  DO equations_row_number=1,vectorMapping%totalNumberOfRows
                                    IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & NUMBER_OF_SOLVER_ROWS>0) THEN
                                      !Get the equations residual contribution
                                      CALL DistributedVector_ValuesGet(RESIDUAL_VECTOR,equations_row_number, &
                                        & RESIDUAL_VALUE,ERR,ERROR,*999)
                                      !Get the linear matrices contribution to the RHS values if there are any
                                      IF(ASSOCIATED(linearMapping)) THEN
                                        LINEAR_VALUE_SUM=0.0_DP
                                        DO equations_matrix_idx2=1,linearMatrices%numberOfLinearMatrices
                                          LINEAR_MATRIX=>linearMatrices%matrices(equations_matrix_idx2)%ptr
                                          LINEAR_TEMP_VECTOR=>LINEAR_MATRIX%tempVector
                                          CALL DistributedVector_ValuesGet(LINEAR_TEMP_VECTOR,equations_row_number, &
                                            & LINEAR_VALUE,ERR,ERROR,*999)
                                          LINEAR_VALUE_SUM=LINEAR_VALUE_SUM+LINEAR_VALUE
                                        ENDDO !equations_matrix_idx2
                                        RESIDUAL_VALUE=RESIDUAL_VALUE+LINEAR_VALUE_SUM
                                      ENDIF
                                      !Loop over the solver rows associated with this equations set residual row
                                      DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                        solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                          & solver_row_idx)
                                        row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                          & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                          & COUPLING_COEFFICIENTS(solver_row_idx)
                                        VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                        !Add in nonlinear residual values
                                        CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                          & ERR,ERROR,*999)
                                      ENDDO !solver_row_idx
                                    ENDIF
                                  ENDDO !equations_row_number
                                ELSE
                                  CALL FlagError("Equations matrices nonlinear matrices is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE IF(ASSOCIATED(linearMapping)) THEN
                                DO equations_row_number=1,vectorMapping%totalNumberOfRows
                                  IF(SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                      & NUMBER_OF_SOLVER_ROWS>0) THEN
                                    LINEAR_VALUE_SUM=0.0_DP
                                    DO equations_matrix_idx=1,linearMatrices%numberOfLinearMatrices
                                      LINEAR_MATRIX=>linearMatrices%matrices(equations_matrix_idx)%ptr
                                      LINEAR_TEMP_VECTOR=>LINEAR_MATRIX%tempVector
                                      CALL DistributedVector_ValuesGet(LINEAR_TEMP_VECTOR,equations_row_number, &
                                        & LINEAR_VALUE,ERR,ERROR,*999)
                                      LINEAR_VALUE_SUM=LINEAR_VALUE_SUM+LINEAR_VALUE
                                    ENDDO !equations_matrix_idx
                                    RESIDUAL_VALUE=LINEAR_VALUE_SUM
                                    !Loop over the solver rows associated with this equations set residual row
                                    DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                      solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                        & solver_row_idx)
                                      row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                        & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                        & COUPLING_COEFFICIENTS(solver_row_idx)
                                      VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                      !Add in nonlinear residual values
                                      CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                        & ERR,ERROR,*999)
                                    ENDDO !solver_row_idx
                                  ENDIF
                                ENDDO !equations_row_number
                              ENDIF
                            ELSE
                              CALL FlagError("Equations equations mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations equations matrices is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                  !Loop over the interface conditions
                  DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                    INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%ptr
                    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                      LAGRANGE_FIELD=>INTERFACE_CONDITION%LAGRANGE%LAGRANGE_FIELD
                      IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                        INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                          INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
                          IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
                            INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                            IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                              SELECT CASE(INTERFACE_CONDITION%METHOD)
                              CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
                                number_of_interface_matrices=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                                number_of_interface_matrices=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES-1
                              ENDSELECT
                              !Calculate the contributions from any interface matrices
                              DO interface_matrix_idx=1,number_of_interface_matrices
                                !Calculate the interface matrix-Lagrange vector product residual contribution
                                INTERFACE_MATRIX=>INTERFACE_MATRICES%matrices(interface_matrix_idx)%ptr
                                IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
                                  interface_variable_type=INTERFACE_MAPPING%LAGRANGE_VARIABLE_TYPE
                                  INTERFACE_VARIABLE=>INTERFACE_MAPPING%LAGRANGE_VARIABLE
                                  IF(ASSOCIATED(INTERFACE_VARIABLE)) THEN
                                    INTERFACE_TEMP_VECTOR=>INTERFACE_MATRIX%TEMP_VECTOR
                                    !Initialise the linear temporary vector to zero
                                    CALL DistributedVector_AllValuesSet(INTERFACE_TEMP_VECTOR,0.0_DP,ERR,ERROR,*999)
                                    NULLIFY(LAGRANGE_VECTOR)
                                    CALL FIELD_PARAMETER_SET_VECTOR_GET(LAGRANGE_FIELD,interface_variable_type, &
                                      & FIELD_VALUES_SET_TYPE,LAGRANGE_VECTOR,ERR,ERROR,*999)
                                    CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                                      & INTERFACE_MATRIX%MATRIX,LAGRANGE_VECTOR,INTERFACE_TEMP_VECTOR,ERR,ERROR,*999)
                                    !Add interface matrix residual contribution to the solver residual
                                    DO interface_row_number=1,INTERFACE_MATRIX%NUMBER_OF_ROWS
                                      IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                                        & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number)%NUMBER_OF_SOLVER_ROWS>0) THEN
                                        !Loop over the solver rows associated with this interface residual row
                                        !Currently earch interface matrix row has only one corresponding solver row number & coupling coefficient
                                        solver_row_number=SOLVER_MAPPING% & 
                                          & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                                          & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number)%SOLVER_ROW
                                        row_coupling_coefficient=SOLVER_MAPPING% &
                                          & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                                          & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number)%COUPLING_COEFFICIENT
                                        CALL DistributedVector_ValuesGet(INTERFACE_TEMP_VECTOR,interface_row_number, &
                                          & RESIDUAL_VALUE,ERR,ERROR,*999)
                                        VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                        !Add in nonlinear residual values
                                        CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                          & ERR,ERROR,*999)
                                      ENDIF
                                    ENDDO !interface_row_number
                                  ELSE
                                    CALL FlagError("Interface variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                  !Calculate the transposed interface matrix-dependent variable product residual contribution
                                  dependent_variable_type=INTERFACE_MAPPING% &
                                    & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE_TYPE
                                  DEPENDENT_VARIABLE=>INTERFACE_MAPPING% &
                                    & INTERFACE_MATRIX_ROWS_TO_VAR_MAPS(interface_matrix_idx)%VARIABLE
                                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                                    INTERFACE_TEMP_VECTOR=>INTERFACE_MATRIX%TEMP_TRANSPOSE_VECTOR
                                    !Initialise the linear temporary vector to zero
                                    CALL DistributedVector_AllValuesSet(INTERFACE_TEMP_VECTOR,0.0_DP,ERR,ERROR,*999)
                                    NULLIFY(DEPENDENT_VECTOR)
                                    DEPENDENT_FIELD=>DEPENDENT_VARIABLE%FIELD
                                    CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,dependent_variable_type, &
                                      & FIELD_VALUES_SET_TYPE,DEPENDENT_VECTOR,ERR,ERROR,*999)
                                    CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                                      & INTERFACE_MATRIX%MATRIX_TRANSPOSE,DEPENDENT_VECTOR,INTERFACE_TEMP_VECTOR,ERR,ERROR,*999)
                                    !Add interface matrix residual contribution to the solver residual.
                                    !The number of columns in the interface matrix is equivalent to the number of rows of the transposed interface matrices
                                    DO interface_row_number=1,INTERFACE_MATRICES%NUMBER_OF_COLUMNS
                                      IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                        & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(interface_row_number)%NUMBER_OF_SOLVER_ROWS>0) THEN
                                        !Loop over the solver rows associated with this interface residual row
                                        !Currently earch interface matrix row has only one corresponding solver row number & coupling coefficient
                                        solver_row_number=SOLVER_MAPPING% & 
                                          & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                          & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(interface_row_number)%SOLVER_ROW
                                        row_coupling_coefficient=SOLVER_MAPPING% & 
                                          & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                          & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(interface_row_number)%COUPLING_COEFFICIENT
                                        CALL DistributedVector_ValuesGet(INTERFACE_TEMP_VECTOR,interface_row_number, &
                                          & RESIDUAL_VALUE,ERR,ERROR,*999)
                                        VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                        !Add in nonlinear residual values
                                        CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                          & ERR,ERROR,*999)
                                      ENDIF
                                    ENDDO !interface_row_number
                                  ELSE
                                    CALL FlagError("Dependent variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  LOCAL_ERROR="Interface matrix is not associated for linear matrix number "// &
                                    & TRIM(NumberToVString(equations_matrix_idx,"*",ERR,ERROR))//"."
                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                ENDIF
                              ENDDO !interface_matrix_idx
                              SELECT CASE(INTERFACE_CONDITION%METHOD)
                              CASE(INTERFACE_CONDITION_PENALTY_METHOD)
                                interface_matrix_idx=INTERFACE_MAPPING%NUMBER_OF_INTERFACE_MATRICES
                                !Calculate the Lagrange-Lagrange vector product residual contribution from the penalty term
                                INTERFACE_MATRIX=>INTERFACE_MATRICES%matrices(interface_matrix_idx)%ptr
                                IF(ASSOCIATED(INTERFACE_MATRIX)) THEN
                                  interface_variable_type=INTERFACE_MAPPING%LAGRANGE_VARIABLE_TYPE
                                  INTERFACE_VARIABLE=>INTERFACE_MAPPING%LAGRANGE_VARIABLE
                                  IF(ASSOCIATED(INTERFACE_VARIABLE)) THEN
                                    INTERFACE_TEMP_VECTOR=>INTERFACE_MATRIX%TEMP_VECTOR
                                    !Initialise the linear temporary vector to zero
                                    CALL DistributedVector_AllValuesSet(INTERFACE_TEMP_VECTOR,0.0_DP,ERR,ERROR,*999)
                                    NULLIFY(LAGRANGE_VECTOR)
                                    CALL FIELD_PARAMETER_SET_VECTOR_GET(LAGRANGE_FIELD,interface_variable_type, &
                                      & FIELD_VALUES_SET_TYPE,LAGRANGE_VECTOR,ERR,ERROR,*999)
                                    CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                                      & INTERFACE_MATRIX%MATRIX,LAGRANGE_VECTOR,INTERFACE_TEMP_VECTOR,ERR,ERROR,*999)
                                    !Add interface matrix residual contribution to the solver residual
                                    DO interface_row_number=1,INTERFACE_MATRIX%NUMBER_OF_ROWS
                                      IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                        & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                                        & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number)%NUMBER_OF_SOLVER_ROWS>0) THEN
                                        !Loop over the solver rows associated with this interface residual row
                                        !Currently earch interface matrix row has only one corresponding solver row number & coupling coefficient
                                        solver_row_number=SOLVER_MAPPING% & 
                                          & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                                          & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number)%SOLVER_ROW
                                        row_coupling_coefficient=SOLVER_MAPPING% &
                                          & INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                          & INTERFACE_TO_SOLVER_MATRIX_MAPS_IM(interface_matrix_idx)% &
                                          & INTERFACE_ROW_TO_SOLVER_ROWS_MAP(interface_row_number)%COUPLING_COEFFICIENT
                                        CALL DistributedVector_ValuesGet(INTERFACE_TEMP_VECTOR,interface_row_number, &
                                          & RESIDUAL_VALUE,ERR,ERROR,*999)
                                        VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                        !Add in nonlinear residual values
                                        CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR,solver_row_number,VALUE, &
                                          & ERR,ERROR,*999)
                                      ENDIF
                                    ENDDO !interface_row_number
                                  ELSE
                                    CALL FlagError("Interface variable is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  LOCAL_ERROR="Interface matrix is not associated for linear matrix number "// &
                                    & TRIM(NumberToVString(equations_matrix_idx,"*",ERR,ERROR))//"."
                                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                ENDIF
                              ENDSELECT
                            ELSE
                              CALL FlagError("Interface mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Interface matrices is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Interface equations is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Interface Lagrange field is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !interface_condition_idx
                  !Start the update the solver residual vector values
                  CALL DistributedVector_UpdateStart(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
                ELSE
                  CALL FlagError("The solver residual vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                CALL DistributedVector_UpdateFinish(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPUTimer(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                userElapsed=USER_TIME2(1)-USER_TIME1(1)
                systemElapsed=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
                  & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
                CALL Profiling_TimingsOutput(1,"Solver residual assembly",userElapsed,systemElapsed,err,error,*999)
              ENDIF
            ENDIF
            NULLIFY(SOLVER_RHS_VECTOR)
            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_ONLY.OR. &
              & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
              !Assemble rhs vector
              IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPUTimer(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
              ENDIF
              IF(SOLVER_MATRICES%UPDATE_RHS_VECTOR) THEN
                SOLVER_RHS_VECTOR=>SOLVER_MATRICES%RHS_VECTOR
                IF(ASSOCIATED(SOLVER_RHS_VECTOR)) THEN
                  !Initialise the RHS to zero
                  CALL DistributedVector_AllValuesSet(SOLVER_RHS_VECTOR,0.0_DP,ERR,ERROR,*999)
                  NULLIFY(CHECK_DATA)
                  CALL DistributedVector_DataGet(SOLVER_RHS_VECTOR,CHECK_DATA,ERR,ERROR,*999)
                  SUBTRACT_FIXED_BCS_FROM_RESIDUAL=.FALSE.
                  IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                      & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                      & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
                    IF(SOLVER_MATRICES%UPDATE_RESIDUAL) THEN
                      IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                        SUBTRACT_FIXED_BCS_FROM_RESIDUAL=.TRUE.
                      ELSE
                        CALL FlagError("The solver residual vector is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ENDIF
                  !Loop over the equations sets
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          NULLIFY(vectorEquations)
                          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                          vectorMatrices=>vectorEquations%vectorMatrices
                          IF(ASSOCIATED(vectorMatrices)) THEN
                            vectorMapping=>vectorEquations%vectorMapping
                            IF(ASSOCIATED(vectorMapping)) THEN
                              sourceMapping=>vectorMapping%sourceMapping
                              IF(ASSOCIATED(sourceMapping)) THEN
                                sourceVector=>vectorMatrices%sourceVector
                                IF(ASSOCIATED(sourceVector)) THEN
                                  DISTRIBUTED_SOURCE_VECTOR=>sourceVector%VECTOR
                                ELSE
                                  CALL FlagError("Source vector vector is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ENDIF
                              rhsMapping=>vectorMapping%rhsMapping
                              IF(ASSOCIATED(rhsMapping)) THEN
                                NULLIFY(RHS_PARAMETERS)
                                RHS_VARIABLE_TYPE=rhsMapping%rhsVariableType
                                CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,RHS_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & RHS_PARAMETERS,ERR,ERROR,*999)
                                NULLIFY(CHECK_DATA)
                                CALL DistributedVector_DataGet(SOLVER_RHS_VECTOR,CHECK_DATA,ERR,ERROR,*999)    
                                rhsVector=>vectorMatrices%rhsVector
                                IF(ASSOCIATED(rhsVector)) THEN
                                  linearMapping=>vectorMapping%linearMapping
                                  nonlinearMapping=>vectorMapping%nonlinearMapping
                                  IF(ASSOCIATED(linearMapping)) THEN
                                    linearMatrices=>vectorMatrices%linearMatrices
                                    IF(ASSOCIATED(linearMatrices)) THEN
                                      ALLOCATE(DEPENDENT_PARAMETERS(linearMapping%numberOfLinearMatrixVariables),STAT=ERR)
                                      IF(ERR/=0) CALL FlagError("Could not allocate dependent_parameters.",ERR,ERROR,*999)
                                      DO variable_idx=1,linearMapping%numberOfLinearMatrixVariables
                                        variable_type=linearMapping%linearMatrixVariableTypes(variable_idx)
                                        NULLIFY(DEPENDENT_PARAMETERS(variable_idx)%ptr)
                                        CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                                          & DEPENDENT_PARAMETERS(variable_idx)%ptr,ERR,ERROR,*999)
                                      ENDDO !variable_idx
                                    ELSE
                                      CALL FlagError("Equations matrices linear matrices is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ENDIF
                                  BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                                  IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
!!TODO: what if the equations set doesn't have a RHS vector???
                                    RHS_VARIABLE=>rhsMapping%rhsVariable
                                    RHS_VARIABLE_TYPE=RHS_VARIABLE%VARIABLE_TYPE
                                    RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
                                    ! Check if there are any integrated values to add
                                    CALL FIELD_PARAMETER_SET_CREATED(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                      & FIELD_INTEGRATED_NEUMANN_SET_TYPE,HAS_INTEGRATED_VALUES,ERR,ERROR,*999)
                                    CALL FIELD_PARAMETER_SET_CREATED(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                      & FIELD_INTEGRATED_ROBIN_SET_TYPE,R_HAS_INTEGRATED_VALUES,ERR,ERROR,*999) !Elias
                                    EQUATIONS_RHS_VECTOR=>rhsVector%VECTOR
                                    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,RHS_VARIABLE, &
                                      & RHS_BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                                    IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN
                                      !Update RHS field by integrating any point Neumann conditions
                                      CALL BoundaryConditions_NeumannIntegrate(RHS_BOUNDARY_CONDITIONS, &
                                        & ERR,ERROR,*999)
                                      CALL BoundaryConditions_RobinIntegrate(RHS_BOUNDARY_CONDITIONS, &
                                        & ERR,ERROR,*999) !Elias
                                      !Loop over the rows in the equations set
                                      DO equations_row_number=1,vectorMapping%totalNumberOfRows
                                        !Get the source vector contribute to the RHS values if there are any
                                        IF(ASSOCIATED(sourceMapping)) THEN
                                          !Add in equations source values
                                          CALL DistributedVector_ValuesGet(DISTRIBUTED_SOURCE_VECTOR,equations_row_number, &
                                            & SOURCE_VALUE,ERR,ERROR,*999)
                                          !Loop over the solver rows associated with this equations set row
                                          DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS

                                            solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                              & solver_row_idx)

                                            row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                              & COUPLING_COEFFICIENTS(solver_row_idx)
                                            VALUE=1.0_DP*SOURCE_VALUE*row_coupling_coefficient
                                            !Calculates the contribution from each row of the equations matrix and adds to solver matrix
                                            CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                              & ERR,ERROR,*999)
                                          ENDDO !solver_row_idx
                                        ENDIF
                                        rhs_variable_dof=rhsMapping%equationsRowToRHSDofMap(equations_row_number)
                                        rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                                        rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%DOF_TYPES(rhs_global_dof)
                                        !Apply boundary conditions
                                        SELECT CASE(rhs_boundary_condition)
                                        CASE(BOUNDARY_CONDITION_DOF_FREE)
                                          !Add in equations RHS values
                                          CALL DistributedVector_ValuesGet(EQUATIONS_RHS_VECTOR,equations_row_number, &
                                            & RHS_VALUE,ERR,ERROR,*999)
                                          IF(HAS_INTEGRATED_VALUES) THEN
                                            !Add any Neumann integrated values, b = f + N q
                                            CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                              & FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                              & ERR,ERROR,*999)
                                            RHS_VALUE=RHS_VALUE+RHS_INTEGRATED_VALUE
                                          END IF

                                          IF(R_HAS_INTEGRATED_VALUES) THEN !Elias /*
                                            !Add any Robin integrated values, b = f + N q+R q_R
                                            CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                              & FIELD_INTEGRATED_ROBIN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                              & ERR,ERROR,*999)
                                            RHS_VALUE=RHS_VALUE+RHS_INTEGRATED_VALUE
                                          END IF ! Elias */

                                          !Loop over the solver rows associated with this equations set row
                                          DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                            & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                            solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                              & solver_row_idx)
                                            row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                              & COUPLING_COEFFICIENTS(solver_row_idx)
                                            VALUE=RHS_VALUE*row_coupling_coefficient
                                            CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                              & ERR,ERROR,*999)
                                          ENDDO !solver_row_idx
                                          !Set Dirichlet boundary conditions
                                          IF(ASSOCIATED(linearMapping).AND..NOT.ASSOCIATED(nonlinearMapping)) THEN
                                            !Loop over the dependent variables associated with this equations set row
                                            DO variable_idx=1,linearMapping%numberOfLinearMatrixVariables
                                              variable_type=linearMapping%linearMatrixVariableTypes(variable_idx)
                                              DEPENDENT_VARIABLE=>linearMapping%varToEquationsMatricesMaps( &
                                                & variable_type)%VARIABLE
                                              DEPENDENT_VARIABLE_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
                                              VARIABLE_DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                                              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE, &
                                                & DEPENDENT_BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                                              variable_dof=linearMapping%equationsRowToVariableDOFMaps( &
                                                & equations_row_number,variable_idx)
                                              variable_global_dof=VARIABLE_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(variable_dof)
                                              variable_boundary_condition=DEPENDENT_BOUNDARY_CONDITIONS%DOF_TYPES( &
                                                & variable_global_dof)
                                              IF(variable_boundary_condition==BOUNDARY_CONDITION_DOF_FIXED) THEN
                                                DEPENDENT_VALUE=DEPENDENT_PARAMETERS(variable_idx)%ptr(variable_dof)
                                                IF(ABS(DEPENDENT_VALUE)>=ZERO_TOLERANCE) THEN
                                                  DO equations_matrix_idx=1,linearMapping%varToEquationsMatricesMaps( &
                                                    & variable_type)%numberOfEquationsMatrices
                                                    equations_matrix_number=linearMapping%varToEquationsMatricesMaps( &
                                                      & variable_type)%equationsMatrixNumbers(equations_matrix_idx)
                                                    equationsMatrix=>linearMatrices%matrices(equations_matrix_number)%ptr
                                                    equations_column_number=linearMapping%varToEquationsMatricesMaps( &
                                                      & variable_type)%dofToColumnsMaps(equations_matrix_idx)%columnDOF( &
                                                      & variable_dof)
                                                    IF(ASSOCIATED(DEPENDENT_BOUNDARY_CONDITIONS%DIRICHLET_BOUNDARY_CONDITIONS)) THEN
                                                      IF(DEPENDENT_BOUNDARY_CONDITIONS%NUMBER_OF_DIRICHLET_CONDITIONS>0) THEN
                                                        DO dirichlet_idx=1,DEPENDENT_BOUNDARY_CONDITIONS% &
                                                          & NUMBER_OF_DIRICHLET_CONDITIONS
                                                          IF(DEPENDENT_BOUNDARY_CONDITIONS%DIRICHLET_BOUNDARY_CONDITIONS% &
                                                            & DIRICHLET_DOF_INDICES(dirichlet_idx)==equations_column_number) EXIT
                                                        ENDDO
                                                        SELECT CASE(equationsMatrix%storageType)
                                                        CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                                          DO dirichlet_row=1,vectorMatrices%totalNumberOfRows
                                                            CALL DistributedMatrix_ValuesGet(equationsMatrix%MATRIX, &
                                                              & dirichlet_row,equations_column_number,MATRIX_VALUE,ERR,ERROR,*999)
                                                            IF(ABS(MATRIX_VALUE)>=ZERO_TOLERANCE) THEN
                                                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                & dirichlet_row)%NUMBER_OF_SOLVER_ROWS
                                                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                  & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                  & dirichlet_row)%SOLVER_ROWS(solver_row_idx)
                                                                row_coupling_coefficient=SOLVER_MAPPING% &
                                                                  & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(dirichlet_row)% &
                                                                  & COUPLING_COEFFICIENTS(solver_row_idx)
                                                                VALUE=-1.0_DP*MATRIX_VALUE*DEPENDENT_VALUE*row_coupling_coefficient
                                                                CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR, &
                                                                  & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                IF(SUBTRACT_FIXED_BCS_FROM_RESIDUAL) THEN
                                                                  CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR, &
                                                                      & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                ENDIF
                                                              ENDDO !solver_row_idx
                                                            ENDIF
                                                          ENDDO !dirichlet_row
                                                        CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                                          dirichlet_row=equations_column_number
                                                          CALL DistributedMatrix_ValuesGet(equationsMatrix%MATRIX, &
                                                            & dirichlet_row,equations_column_number,MATRIX_VALUE,ERR,ERROR,*999)
                                                          IF(ABS(MATRIX_VALUE)>=ZERO_TOLERANCE) THEN
                                                            DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                              & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                              & dirichlet_row)%NUMBER_OF_SOLVER_ROWS
                                                              solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                & dirichlet_row)%SOLVER_ROWS(solver_row_idx)
                                                              row_coupling_coefficient=SOLVER_MAPPING% &
                                                                & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(dirichlet_row)% &
                                                                & COUPLING_COEFFICIENTS(solver_row_idx)
                                                              VALUE=-1.0_DP*MATRIX_VALUE*DEPENDENT_VALUE*row_coupling_coefficient
                                                              CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR, &
                                                                & solver_row_number,VALUE,ERR,ERROR,*999)
                                                              IF(SUBTRACT_FIXED_BCS_FROM_RESIDUAL) THEN
                                                                CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR, &
                                                                    & solver_row_number,VALUE,ERR,ERROR,*999)
                                                              ENDIF
                                                            ENDDO !solver_row_idx
                                                          ENDIF
                                                        CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                        CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                                          SPARSITY_INDICES=>DEPENDENT_BOUNDARY_CONDITIONS% &
                                                            & DIRICHLET_BOUNDARY_CONDITIONS%LINEAR_SPARSITY_INDICES( &
                                                            & equations_set_idx,equations_matrix_idx)%ptr
                                                          IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                                            DO equations_row_number2=SPARSITY_INDICES%SPARSE_COLUMN_INDICES( &
                                                              & dirichlet_idx),SPARSITY_INDICES%SPARSE_COLUMN_INDICES( &
                                                              & dirichlet_idx+1)-1
                                                              dirichlet_row=SPARSITY_INDICES%SPARSE_ROW_INDICES( &
                                                                & equations_row_number2)
                                                              CALL DistributedMatrix_ValuesGet(equationsMatrix%MATRIX, &
                                                                & dirichlet_row,equations_column_number,MATRIX_VALUE,ERR,ERROR,*999)
                                                              IF(ABS(MATRIX_VALUE)>=ZERO_TOLERANCE) THEN
                                                                DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                  & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                  & dirichlet_row)%NUMBER_OF_SOLVER_ROWS
                                                                  solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                    & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                    & dirichlet_row)%SOLVER_ROWS(solver_row_idx)
                                                                  row_coupling_coefficient=SOLVER_MAPPING% &
                                                                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(dirichlet_row)% &
                                                                    & COUPLING_COEFFICIENTS(solver_row_idx)
                                                                  VALUE=-1.0_DP*MATRIX_VALUE*DEPENDENT_VALUE* &
                                                                    & row_coupling_coefficient
                                                                  CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR, &
                                                                    & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                  IF(SUBTRACT_FIXED_BCS_FROM_RESIDUAL) THEN
                                                                    CALL DistributedVector_ValuesAdd(SOLVER_RESIDUAL_VECTOR, &
                                                                        & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                  ENDIF
                                                                ENDDO !solver_row_idx
                                                              ENDIF
                                                            ENDDO !equations_row_number2
                                                          ELSE
                                                            CALL FlagError("Sparsity indices are not associated.",ERR,ERROR,*999)
                                                          ENDIF
                                                        CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                        CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                        CASE DEFAULT
                                                          LOCAL_ERROR="The storage type of "// &
                                                            & TRIM(NumberToVString(equationsMatrix%storageType,"*", &
                                                            & ERR,ERROR))//" is invalid."
                                                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                                        END SELECT
                                                      ENDIF
                                                    ELSE
                                                      CALL FlagError("Dirichlet boundary conditions is not associated.",ERR, &
                                                        & ERROR,*999)
                                                    ENDIF
                                                  ENDDO !matrix_idx
                                                ENDIF
                                              ENDIF
                                            ENDDO !variable_idx
                                          ENDIF
                                        CASE(BOUNDARY_CONDITION_DOF_FIXED)
                                          RHS_VALUE=RHS_PARAMETERS(rhs_variable_dof)
                                          ! Add any integrated RHS values calculated from point Neumann conditions, b = f + N q
                                          IF(HAS_INTEGRATED_VALUES) THEN
                                            CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                              & FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                              & ERR,ERROR,*999)
                                            RHS_VALUE=RHS_VALUE+RHS_INTEGRATED_VALUE
                                          END IF

                                          IF(R_HAS_INTEGRATED_VALUES) THEN ! Elias */
                                            !Add any Robin integrated values, b = f + N q + R q_h
                                            CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                              & FIELD_INTEGRATED_ROBIN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                              & ERR,ERROR,*999)
                                            VALUE=VALUE+RHS_INTEGRATED_VALUE*row_coupling_coefficient
                                          END IF !Elias /*

                                          IF(ABS(RHS_VALUE)>=ZERO_TOLERANCE) THEN
                                            !Loop over the solver rows associated with this equations set row
                                            DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                              solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                                & solver_row_idx)
                                              row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                & COUPLING_COEFFICIENTS(solver_row_idx)
                                              !For nonlinear problems, f(x) - b = 0, and for linear, K x = b, so we always add the
                                              !right hand side field value to the solver right hand side vector
                                              VALUE=RHS_VALUE*row_coupling_coefficient
                                              CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                                & ERR,ERROR,*999)
                                            ENDDO !solver_row_idx
                                          ENDIF
                                        CASE(BOUNDARY_CONDITION_DOF_MIXED)
                                          !Set Robin or is it Cauchy??? boundary conditions
                                          CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The RHS boundary condition of "// &
                                            & TRIM(NumberToVString(rhs_boundary_condition,"*",ERR,ERROR))// &
                                            & " for RHS variable dof number "// &
                                            & TRIM(NumberToVString(rhs_variable_dof,"*",ERR,ERROR))//" is invalid."
                                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                        END SELECT
                                      ENDDO !equations_row_number
                                      IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                                        CALL DistributedVector_UpdateStart(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
                                        CALL DistributedVector_UpdateFinish(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)
                                      ENDIF
                                      NULLIFY(CHECK_DATA2)
                                      CALL DistributedVector_DataGet(EQUATIONS_RHS_VECTOR,CHECK_DATA2,ERR,ERROR,*999)    
                                      NULLIFY(CHECK_DATA3)
                                      CALL DistributedVector_DataGet(SOLVER_RHS_VECTOR,CHECK_DATA3,ERR,ERROR,*999)    
                                      NULLIFY(CHECK_DATA4)
                                      CALL DistributedVector_DataGet(SOLVER_RHS_VECTOR,CHECK_DATA4,ERR,ERROR,*999)    
                                    ELSE
                                      CALL FlagError("RHS boundary conditions variable is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Equations set boundary conditions is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                  IF(ASSOCIATED(linearMapping)) THEN
                                    DO variable_idx=1,linearMapping%numberOfLinearMatrixVariables
                                      variable_type=linearMapping%linearMatrixVariableTypes(variable_idx)
                                      CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,variable_type,FIELD_VALUES_SET_TYPE, &
                                        & DEPENDENT_PARAMETERS(variable_idx)%ptr,ERR,ERROR,*999)
                                    ENDDO !variable_idx
                                    IF(ALLOCATED(DEPENDENT_PARAMETERS)) DEALLOCATE(DEPENDENT_PARAMETERS)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Equations matrices RHS vector is not associated.",ERR,ERROR,*999)
                                ENDIF
                                CALL FIELD_PARAMETER_SET_DATA_RESTORE(DEPENDENT_FIELD,RHS_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & RHS_PARAMETERS,ERR,ERROR,*999)
                              ELSE
                                CALL FlagError("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Equations equations mapping is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Equations equations matrices is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                      ENDIF
                    ENDIF
                  ENDDO !equations_set_idx
                  !Add in any rows from any interface conditions
                  DO interface_condition_idx=1,SOLVER_MAPPING%NUMBER_OF_INTERFACE_CONDITIONS
                    INTERFACE_CONDITION=>SOLVER_MAPPING%INTERFACE_CONDITIONS(interface_condition_idx)%ptr
                    IF(ASSOCIATED(INTERFACE_CONDITION)) THEN
                      SELECT CASE(INTERFACE_CONDITION%METHOD)
                      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
                        INTERFACE_EQUATIONS=>INTERFACE_CONDITION%INTERFACE_EQUATIONS
                        IF(ASSOCIATED(INTERFACE_EQUATIONS)) THEN
                          INTERFACE_MAPPING=>INTERFACE_EQUATIONS%INTERFACE_MAPPING
                          IF(ASSOCIATED(INTERFACE_MAPPING)) THEN
                            INTERFACE_LAGRANGE=>INTERFACE_CONDITION%LAGRANGE
                            IF(ASSOCIATED(INTERFACE_LAGRANGE)) THEN
                              LAGRANGE_FIELD=>INTERFACE_LAGRANGE%LAGRANGE_FIELD
                              IF(ASSOCIATED(LAGRANGE_FIELD)) THEN
                                INTERFACE_RHS_MAPPING=>INTERFACE_MAPPING%RHS_MAPPING
                                IF(ASSOCIATED(INTERFACE_RHS_MAPPING)) THEN
                                  INTERFACE_MATRICES=>INTERFACE_EQUATIONS%INTERFACE_MATRICES
                                  IF(ASSOCIATED(INTERFACE_MATRICES)) THEN
                                    INTERFACE_RHS_VECTOR=>INTERFACE_MATRICES%RHS_VECTOR
                                    IF(ASSOCIATED(INTERFACE_RHS_VECTOR)) THEN
                                      !Worry about BCs on the Lagrange variables later.
                                      DO interface_column_number=1,INTERFACE_MAPPING%TOTAL_NUMBER_OF_COLUMNS
                                        CALL DistributedVector_ValuesGet(INTERFACE_RHS_VECTOR%RHS_VECTOR, &
                                          & interface_column_number,RHS_VALUE,ERR,ERROR,*999)
                                        !Loop over the solver rows this interface column is mapped to
                                        DO solver_row_idx=1,SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                                          & interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS( &
                                          & interface_column_number)%NUMBER_OF_SOLVER_ROWS
                                          solver_row_number=SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                                            & interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS( &
                                            & interface_column_number)%SOLVER_ROW
                                          row_coupling_coefficient=SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP( &
                                            & interface_condition_idx)%INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS( &
                                            & interface_column_number)%COUPLING_COEFFICIENT
                                          VALUE=RHS_VALUE*row_coupling_coefficient
                                          CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                            & ERR,ERROR,*999)
                                        ENDDO !solver_row_idx
                                      ENDDO !interface_column_idx
                                    ELSE
                                      CALL FlagError("Interface matrices RHS vector is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Interface equations interface matrices is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Interface mapping RHS mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Interface Lagrange field is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Interface Lagrange is not associated.",ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("Interface equations interface mapping is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Interface condition equations is not associated.",ERR,ERROR,*999)
                        ENDIF
                      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
                        CALL FlagError("Not implemented.",ERR,ERROR,*999)
                      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
                        CALL FlagError("Not implemented.",ERR,ERROR,*999)
                      CASE DEFAULT
                        LOCAL_ERROR="The interface condition method of "// &
                          & TRIM(NumberToVString(INTERFACE_CONDITION%METHOD,"*",ERR,ERROR))// &
                          & " is invalid."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Interface condition is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ENDDO !interface_condition_idx
                  !Start the update the solver RHS vector values
                  CALL DistributedVector_UpdateStart(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
                  NULLIFY(CHECK_DATA)
                  CALL DistributedVector_DataGet(SOLVER_RHS_VECTOR,CHECK_DATA,ERR,ERROR,*999)   
                ELSE
                  CALL FlagError("The solver RHS vector is not associated.",ERR,ERROR,*999)
                ENDIF
              ENDIF
              IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                CALL CPUTimer(USER_CPU,USER_TIME2,ERR,ERROR,*999)
                CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME2,ERR,ERROR,*999)
                userElapsed=USER_TIME2(1)-USER_TIME1(1)
                systemElapsed=SYSTEM_TIME2(1)-SYSTEM_TIME1(1)
                IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
                  & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
                CALL Profiling_TimingsOutput(1,"Solver RHS assembly",userElapsed,systemElapsed,err,error,*999)
              ENDIF
            ENDIF
            IF(ASSOCIATED(SOLVER_RHS_VECTOR)) THEN
              CALL DistributedVector_UpdateFinish(SOLVER_RHS_VECTOR,ERR,ERROR,*999)
            ENDIF
            !If required output the solver matrices
            IF(SOLVER%outputType>=SOLVER_MATRIX_OUTPUT) THEN
              CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SELECTION_TYPE,SOLVER_MATRICES,ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver solver matrices is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver matrices solution mapping is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solver equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("SOLVER_MATRICES_STATIC_ASSEMBLE")
    RETURN
999 IF(ALLOCATED(DEPENDENT_PARAMETERS)) DEALLOCATE(DEPENDENT_PARAMETERS)    
    ERRORSEXITS("SOLVER_MATRICES_STATIC_ASSEMBLE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_STATIC_ASSEMBLE
