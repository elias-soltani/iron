  !>Assembles the solver matrices and rhs from the dynamic equations.
  SUBROUTINE SOLVER_MATRICES_DYNAMIC_ASSEMBLE(SOLVER,SELECTION_TYPE,ERR,ERROR,*)
    
    !Argument variable
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(IN) :: SELECTION_TYPE !<The type of matrix selection \see SOLVER_MATRICES_ROUTINES_SelectMatricesTypes,SOLVER_MATRICES_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DYNAMIC_VARIABLE_TYPE,equations_matrix_idx,equations_row_number,equations_set_idx,LINEAR_VARIABLE_TYPE, &
      & rhs_boundary_condition,rhs_global_dof,rhs_variable_dof,rhs_variable_type,solver_row_idx,solver_row_number, &
      & solver_matrix_idx,residual_variable_dof,variable_boundary_condition,variable_type,equations_matrix_idx2, &
      & variable_idx,variable_global_dof,variable_dof,equations_row_number2,equations_matrix_number,DEPENDENT_VARIABLE_TYPE, &
      & equations_column_number,dirichlet_row,dirichlet_idx, &
      & interface_condition_idx,interface_matrix_idx,interface_column_number,interface_row_number, &
      & interface_variable_type,number_of_interface_matrices
    REAL(SP) :: systemElapsed,SYSTEM_TIME1(1),SYSTEM_TIME2(1),userElapsed,USER_TIME1(1),USER_TIME2(1)
    REAL(DP) :: DAMPING_MATRIX_COEFFICIENT,DELTA_T,DYNAMIC_VALUE,FIRST_UPDATE_FACTOR,RESIDUAL_VALUE, &
      & LINEAR_VALUE,LINEAR_VALUE_SUM,MASS_MATRIX_COEFFICIENT,RHS_VALUE,row_coupling_coefficient,PREVIOUS_RESIDUAL_VALUE, &
      & SECOND_UPDATE_FACTOR,SOURCE_VALUE,STIFFNESS_MATRIX_COEFFICIENT,VALUE,JACOBIAN_MATRIX_COEFFICIENT,ALPHA_VALUE, &
      & MATRIX_VALUE,DYNAMIC_DISPLACEMENT_FACTOR,DYNAMIC_VELOCITY_FACTOR,DYNAMIC_ACCELERATION_FACTOR,RHS_INTEGRATED_VALUE
    REAL(DP) :: MatrixCoefficients(2)=[0.0_DP,0.0_DP]
    REAL(DP), POINTER :: FIELD_VALUES_VECTOR(:),PREVIOUS_VALUES_VECTOR(:),PREVIOUS_VELOCITY_VECTOR(:), &
      & PREVIOUS_ACCELERATION_VECTOR(:),RHS_PARAMETERS(:)
    LOGICAL :: HAS_INTEGRATED_VALUES,R_HAS_INTEGRATED_VALUES !Elias
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: RHS_BOUNDARY_CONDITIONS,DEPENDENT_BOUNDARY_CONDITIONS
    TYPE(DistributedMatrixType), POINTER :: PREVIOUS_SOLVER_DISTRIBUTED_MATRIX,SOLVER_DISTRIBUTED_MATRIX
    TYPE(DistributedVectorType), POINTER :: DEPENDENT_VECTOR,DYNAMIC_TEMP_VECTOR,EQUATIONS_RHS_VECTOR,DISTRIBUTED_SOURCE_VECTOR, &
      & LINEAR_TEMP_VECTOR,PREDICTED_MEAN_ACCELERATION_VECTOR,PREDICTED_MEAN_DISPLACEMENT_VECTOR,PREDICTED_MEAN_VELOCITY_VECTOR, &
      & SOLVER_RHS_VECTOR, SOLVER_RESIDUAL_VECTOR,RESIDUAL_VECTOR,INCREMENTAL_VECTOR,INTERFACE_TEMP_VECTOR, &
      & LAGRANGE_VECTOR
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: RHS_DOMAIN_MAPPING,VARIABLE_DOMAIN_MAPPING
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: DAMPING_MATRIX,LINEAR_MATRIX,MASS_MATRIX,STIFFNESS_MATRIX,equationsMatrix
    TYPE(EquationsJacobianType), POINTER :: JACOBIAN_MATRIX
    TYPE(JACOBIAN_TO_SOLVER_MAP_TYPE), POINTER :: JACOBIAN_TO_SOLVER_MAP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,LAGRANGE_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DYNAMIC_VARIABLE,LINEAR_VARIABLE,RHS_VARIABLE,INTERFACE_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(BOUNDARY_CONDITIONS_SPARSITY_INDICES_TYPE), POINTER :: SPARSITY_INDICES
    
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: INTERFACE_CONDITION
    TYPE(INTERFACE_EQUATIONS_TYPE), POINTER :: INTERFACE_EQUATIONS
    TYPE(INTERFACE_LAGRANGE_TYPE), POINTER :: INTERFACE_LAGRANGE
    TYPE(INTERFACE_MAPPING_TYPE), POINTER :: INTERFACE_MAPPING
    TYPE(INTERFACE_MAPPING_RHS_TYPE), POINTER :: INTERFACE_RHS_MAPPING
    TYPE(INTERFACE_MATRICES_TYPE), POINTER :: INTERFACE_MATRICES
    TYPE(INTERFACE_MATRIX_TYPE), POINTER :: INTERFACE_MATRIX
    TYPE(INTERFACE_RHS_TYPE), POINTER :: INTERFACE_RHS_VECTOR
    TYPE(INTERFACE_TO_SOLVER_MAPS_TYPE), POINTER :: INTERFACE_TO_SOLVER_MAP

    REAL(DP), POINTER :: CHECK_DATA(:),PREVIOUS_RESIDUAL_PARAMETERS(:),CHECK_DATA2(:)
    !STABILITY_TEST under investigation
    LOGICAL :: STABILITY_TEST
    !.FALSE. guarantees weighting as described in OpenCMISS notes
    !.TRUE. weights mean predicted field rather than the whole NL contribution
    !-> to be removed later
    STABILITY_TEST=.FALSE.
   
    ENTERS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE",ERR,ERROR,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      NULLIFY(DYNAMIC_SOLVER)
      NULLIFY(SOLVER_EQUATIONS)
      NULLIFY(SOLVER_MAPPING)
      NULLIFY(SOLVER_MATRICES)
      !
      NULLIFY(BOUNDARY_CONDITIONS)
      NULLIFY(RHS_BOUNDARY_CONDITIONS)
      NULLIFY(DEPENDENT_BOUNDARY_CONDITIONS)
      NULLIFY(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)
      NULLIFY(SOLVER_DISTRIBUTED_MATRIX)
      NULLIFY(DEPENDENT_VECTOR)
      NULLIFY(DYNAMIC_TEMP_VECTOR)
      NULLIFY(EQUATIONS_RHS_VECTOR)
      NULLIFY(DISTRIBUTED_SOURCE_VECTOR)
      NULLIFY(LINEAR_TEMP_VECTOR)
      NULLIFY(PREDICTED_MEAN_ACCELERATION_VECTOR)
      NULLIFY(PREDICTED_MEAN_DISPLACEMENT_VECTOR)
      NULLIFY(PREDICTED_MEAN_VELOCITY_VECTOR)
      NULLIFY(SOLVER_RHS_VECTOR)
      NULLIFY(SOLVER_RESIDUAL_VECTOR)
      NULLIFY(RESIDUAL_VECTOR)
      NULLIFY(INCREMENTAL_VECTOR)
      NULLIFY(RHS_DOMAIN_MAPPING)
      NULLIFY(VARIABLE_DOMAIN_MAPPING)
      NULLIFY(EQUATIONS)
      NULLIFY(vectorMapping)
      NULLIFY(dynamicMapping)
      NULLIFY(nonlinearMapping)
      NULLIFY(linearMapping)
      NULLIFY(rhsMapping)
      NULLIFY(sourceMapping)
      NULLIFY(vectorMatrices)
      NULLIFY(dynamicMatrices)
      NULLIFY(nonlinearMatrices)
      NULLIFY(linearMatrices)
      NULLIFY(rhsVector)
      NULLIFY(sourceVector)
      NULLIFY(DAMPING_MATRIX)
      NULLIFY(LINEAR_MATRIX)
      NULLIFY(MASS_MATRIX)
      NULLIFY(STIFFNESS_MATRIX)
      NULLIFY(equationsMatrix)
      NULLIFY(JACOBIAN_MATRIX)
      NULLIFY(JACOBIAN_TO_SOLVER_MAP)
      NULLIFY(EQUATIONS_SET)
      NULLIFY(DEPENDENT_FIELD)
      NULLIFY(LAGRANGE_FIELD)
      NULLIFY(DYNAMIC_VARIABLE)
      NULLIFY(LINEAR_VARIABLE)
      NULLIFY(RHS_VARIABLE)
      NULLIFY(DEPENDENT_VARIABLE)
      NULLIFY(SOLVER_MATRIX)
      NULLIFY(INTERFACE_CONDITION)
      NULLIFY(INTERFACE_EQUATIONS)
      NULLIFY(INTERFACE_LAGRANGE)
      NULLIFY(INTERFACE_MAPPING)
      NULLIFY(INTERFACE_RHS_MAPPING)
      NULLIFY(INTERFACE_MATRICES)
      NULLIFY(INTERFACE_MATRIX)
      NULLIFY(INTERFACE_RHS_VECTOR)
      NULLIFY(INTERFACE_TO_SOLVER_MAP)
      NULLIFY(CHECK_DATA)
      NULLIFY(PREVIOUS_RESIDUAL_PARAMETERS)
      NULLIFY(CHECK_DATA2)
      
      !Determine which dynamic solver needs to be used
      IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
        DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
      ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN 
        DYNAMIC_SOLVER=>SOLVER%LINKING_SOLVER%DYNAMIC_SOLVER
      ELSE
        CALL FlagError("Dynamic solver solve type is not associated.",ERR,ERROR,*999)
      END IF
      IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
        IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) THEN
          DELTA_T=DYNAMIC_SOLVER%TIME_INCREMENT
          SELECT CASE(DYNAMIC_SOLVER%DEGREE)
          CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
            STIFFNESS_MATRIX_COEFFICIENT=1.0_DP*DYNAMIC_SOLVER%THETA(1)*DELTA_T
            DAMPING_MATRIX_COEFFICIENT=1.0_DP
            MASS_MATRIX_COEFFICIENT=0.0_DP
            JACOBIAN_MATRIX_COEFFICIENT=STIFFNESS_MATRIX_COEFFICIENT
            DYNAMIC_DISPLACEMENT_FACTOR=DELTA_T
          CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
            STIFFNESS_MATRIX_COEFFICIENT=1.0_DP*(DYNAMIC_SOLVER%THETA(2)*DELTA_T*DELTA_T)/2.0_DP
            DAMPING_MATRIX_COEFFICIENT=1.0_DP*DYNAMIC_SOLVER%THETA(1)*DELTA_T
            MASS_MATRIX_COEFFICIENT=1.0_DP
            JACOBIAN_MATRIX_COEFFICIENT=STIFFNESS_MATRIX_COEFFICIENT
            FIRST_UPDATE_FACTOR=DELTA_T
            DYNAMIC_DISPLACEMENT_FACTOR=DELTA_T*DELTA_T/2.0_DP
            DYNAMIC_VELOCITY_FACTOR=DELTA_T
          CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
            STIFFNESS_MATRIX_COEFFICIENT=1.0_DP*(DYNAMIC_SOLVER%THETA(3)*DELTA_T*DELTA_T*DELTA_T)/6.0_DP
            DAMPING_MATRIX_COEFFICIENT=1.0_DP*(DYNAMIC_SOLVER%THETA(2)*DELTA_T*DELTA_T)/2.0_DP
            MASS_MATRIX_COEFFICIENT=1.0_DP*DYNAMIC_SOLVER%THETA(1)*DELTA_T
            JACOBIAN_MATRIX_COEFFICIENT=STIFFNESS_MATRIX_COEFFICIENT
            FIRST_UPDATE_FACTOR=DELTA_T
            SECOND_UPDATE_FACTOR=DELTA_T*DELTA_T/2.0_DP
            DYNAMIC_DISPLACEMENT_FACTOR=DELTA_T*DELTA_T*DELTA_T/6.0_DP
            DYNAMIC_VELOCITY_FACTOR=DELTA_T*DELTA_T/2.0_DP
            DYNAMIC_ACCELERATION_FACTOR=DELTA_T
          CASE DEFAULT
            LOCAL_ERROR="The dynamic solver degree of "//TRIM(NumberToVString(DYNAMIC_SOLVER%DEGREE,"*",ERR,ERROR))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ENDIF
        SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
        IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
          SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
          IF(ASSOCIATED(SOLVER_MAPPING)) THEN
            SOLVER_MATRICES=>SOLVER_EQUATIONS%SOLVER_MATRICES
            IF(ASSOCIATED(SOLVER_MATRICES)) THEN
              !Assemble the solver matrices
              NULLIFY(PREVIOUS_SOLVER_DISTRIBUTED_MATRIX)
              NULLIFY(SOLVER_MATRIX)
              NULLIFY(SOLVER_DISTRIBUTED_MATRIX)
              NULLIFY(EQUATIONS)
              NULLIFY(vectorMatrices)
              NULLIFY(dynamicMatrices)
              NULLIFY(vectorMapping)
              NULLIFY(dynamicMapping)
              NULLIFY(STIFFNESS_MATRIX)
              NULLIFY(DAMPING_MATRIX)
              NULLIFY(MASS_MATRIX)

              IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
                IF(DYNAMIC_SOLVER%SOLVER_INITIALISED.OR.(.NOT.DYNAMIC_SOLVER%SOLVER_INITIALISED.AND. &
                  & ((DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_FIRST_ORDER.AND.DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
                  & (DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_SECOND_ORDER.AND.DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_SECOND_DEGREE)))) &
                  & THEN
                  !Assemble solver matrices
                  IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                    CALL CPUTimer(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                    CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
                  ENDIF
!               DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
!                 SOLVER_MATRIX=>SOLVER_MATRICES%matrices(solver_matrix_idx)%ptr
!                 END DO

                  solver_matrix_idx=1
                  IF(SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES==solver_matrix_idx) THEN
                    SOLVER_MATRIX=>SOLVER_MATRICES%matrices(1)%ptr
                    IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                      IF(SOLVER_MATRIX%UPDATE_MATRIX) THEN      
                        SOLVER_DISTRIBUTED_MATRIX=>SOLVER_MATRIX%MATRIX
                        IF(ASSOCIATED(SOLVER_DISTRIBUTED_MATRIX)) THEN
                          !Initialise matrix to zero
                          CALL DistributedMatrix_AllValuesSet(SOLVER_DISTRIBUTED_MATRIX,0.0_DP,ERR,ERROR,*999)
                          !Loop over the equations sets
                          DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                            EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)%EQUATIONS
                            IF(ASSOCIATED(EQUATIONS)) THEN
                              NULLIFY(vectorEquations)
                              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                              vectorMapping=>vectorEquations%vectorMapping
                              IF(ASSOCIATED(vectorMapping)) THEN
                                dynamicMapping=>vectorMapping%dynamicMapping
                                IF(ASSOCIATED(dynamicMapping)) THEN
                                  vectorMatrices=>vectorEquations%vectorMatrices
                                  IF(ASSOCIATED(vectorMatrices)) THEN
                                    dynamicMatrices=>vectorMatrices%dynamicMatrices
                                    IF(ASSOCIATED(dynamicMatrices)) THEN
                                      IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) THEN

                                        IF(dynamicMapping%stiffnessMatrixNumber/=0) THEN
                                          STIFFNESS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%stiffnessMatrixNumber)%ptr
                                          IF(ASSOCIATED(STIFFNESS_MATRIX)) THEN


                                            rhsMapping=>vectorMapping%rhsMapping !Elias /*
                                            IF(ASSOCIATED(rhsMapping)) THEN
                                              BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS 
                                              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN

!                                                rhs_variable_type=rhsMapping%rhsVariableType
                                                RHS_VARIABLE=>rhsMapping%rhsVariable
!                                                RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
!                                                CALL FIELD_PARAMETER_SET_CREATED(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
!                                                  & FIELD_INTEGRATED_ROBIN_SET_TYPE,R_HAS_INTEGRATED_VALUES,ERR,ERROR,*999)
                                                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,RHS_VARIABLE, &
                                                  & RHS_BOUNDARY_CONDITIONS,ERR,ERROR,*999) 
                                                IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN

                                                  !Add K_h to the stiffness matrix=K_e+K_h. 
                                                    !Later on this updated matrix is used in solver_rhs_vector as well.
                                                  CALL BoundaryConditions_RobinDynamicIntegrate(RHS_BOUNDARY_CONDITIONS, &
                                                    & STIFFNESS_MATRIX,ERR,ERROR,*999)   
                                                ELSE
                                                  CALL FlagError("RHS boundary conditions variable is not associated.",ERR,ERROR,*999)
                                                ENDIF
                                              ELSE
                                                CALL FlagError("Equations set boundary conditions is not associated.",ERR,ERROR,*999)
                                              ENDIF
                                            ELSE
                                              CALL FlagError("Equations mapping RHS mapping is not associated.",ERR,ERROR,*999)
                                            ENDIF !Elias */


                                            CALL SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx, &
                                              & STIFFNESS_MATRIX_COEFFICIENT,STIFFNESS_MATRIX,ERR,ERROR,*999)

                                          ELSE
                                            CALL FlagError("Dynamic stiffness matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF

                                        IF(dynamicMapping%dampingMatrixNumber/=0) THEN
                                          DAMPING_MATRIX=>dynamicMatrices%matrices(dynamicMapping%dampingMatrixNumber)%ptr
                                          IF(ASSOCIATED(DAMPING_MATRIX)) THEN
                                            CALL SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx, &
                                              & DAMPING_MATRIX_COEFFICIENT,DAMPING_MATRIX,ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic damping matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF

                                        IF(dynamicMapping%massMatrixNumber/=0) THEN
                                          MASS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%massMatrixNumber)%ptr
                                          IF(ASSOCIATED(MASS_MATRIX)) THEN
                                            CALL SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx, &
                                              & MASS_MATRIX_COEFFICIENT,MASS_MATRIX,ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic mass matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF

                                      ELSE
                                        IF(DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_SECOND_ORDER.AND. &
                                          & DYNAMIC_SOLVER%DEGREE==SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                                          IF(dynamicMapping%massMatrixNumber/=0) THEN
                                            MASS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%massMatrixNumber)%ptr
                                            IF(ASSOCIATED(MASS_MATRIX)) THEN
                                              CALL SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx, &
                                                & -1.0_DP,MASS_MATRIX,ERR,ERROR,*999)
                                            ELSE
                                              CALL FlagError("Dynamic stiffness matrix is not associated.",ERR,ERROR,*999)
                                            ENDIF
                                          ELSE
                                            CALL FlagError("Can not perform initial solve with no mass matrix.",ERR,ERROR,*999)
                                          ENDIF
                                        ELSE
                                          IF(dynamicMapping%dampingMatrixNumber/=0) THEN
                                            DAMPING_MATRIX=>dynamicMatrices%matrices(dynamicMapping%dampingMatrixNumber)%ptr
                                            IF(ASSOCIATED(DAMPING_MATRIX)) THEN
                                              CALL SOLVER_MATRIX_EQUATIONS_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx, &
                                                & -1.0_DP,DAMPING_MATRIX,ERR,ERROR,*999)
                                            ELSE
                                              CALL FlagError("Dynamic damping matrix is not associated.",ERR,ERROR,*999)
                                            ENDIF
                                          ELSE
                                            CALL FlagError("Can not perform initial solve with no damping matrix.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Equations equations matrices is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                ELSE
                                  nonlinearMapping=>vectorMapping%nonlinearMapping
                                  IF(ASSOCIATED(nonlinearMapping)) THEN
                                    vectorMatrices=>vectorEquations%vectorMatrices
                                    IF(.NOT.ASSOCIATED(vectorMatrices)) THEN
                                      CALL FlagError("Equations matrices not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                                  ENDIF
                                  !CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Equations equations mapping is not associated.",ERR,ERROR,*999)
                              ENDIF
                            ELSE
                              LOCAL_ERROR="Solver mapping equations is not associated for equations set number "// &
                                & TRIM(NumberToVString(equations_set_idx,"*",ERR,ERROR))//"."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                            NULLIFY(JACOBIAN_TO_SOLVER_MAP)
                            NULLIFY(JACOBIAN_MATRIX)
                            IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                              & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                              & SELECTION_TYPE==SOLVER_MATRICES_JACOBIAN_ONLY) THEN

                              !Now set the values from the equations Jacobian
                              nonlinearMatrices=>vectorMatrices%nonlinearMatrices
                              IF(ASSOCIATED(nonlinearMatrices)) THEN
                                DO equations_matrix_idx=1,nonlinearMatrices%numberOfJacobians
                                  JACOBIAN_TO_SOLVER_MAP=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                    & EQUATIONS_TO_SOLVER_MATRIX_MAPS_SM(solver_matrix_idx)%JACOBIAN_TO_SOLVER_MATRIX_MAPS( &
                                    & equations_matrix_idx)%ptr
                                  IF(ASSOCIATED(JACOBIAN_TO_SOLVER_MAP)) THEN
                                    JACOBIAN_MATRIX=>JACOBIAN_TO_SOLVER_MAP%JACOBIAN_MATRIX
                                    IF(ASSOCIATED(JACOBIAN_MATRIX)) THEN
                                      CALL SOLVER_MATRIX_JACOBIAN_MATRIX_ADD(SOLVER_MATRIX,equations_set_idx, & 
                                        & JACOBIAN_MATRIX_COEFFICIENT,JACOBIAN_MATRIX,ERR,ERROR,*999)
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
                                  SELECT CASE(INTERFACE_MATRIX%INTERFACE_MATRIX_TIME_DEPENDENCE_TYPE)
                                  CASE(INTERFACE_MATRIX_STATIC)
                                    MatrixCoefficients(1)=STIFFNESS_MATRIX_COEFFICIENT
                                  CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                                    MatrixCoefficients(1)=DAMPING_MATRIX_COEFFICIENT
                                  CASE DEFAULT
                                    CALL FlagError("Not implemented.",Err,Error,*999)
                                  END SELECT
                                  IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                    SELECT CASE(INTERFACE_MATRIX%INTERFACE_MATRIX_TRANSPOSE_TIME_DEPENDENCE_TYPE)
                                    CASE(INTERFACE_MATRIX_STATIC)
                                      MatrixCoefficients(2)=STIFFNESS_MATRIX_COEFFICIENT
                                    CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                                      MatrixCoefficients(2)=DAMPING_MATRIX_COEFFICIENT
                                    CASE DEFAULT
                                      CALL FlagError("Not implemented.",Err,Error,*999)
                                    END SELECT
                                  ELSE
                                    MatrixCoefficients(2)=0.0_DP
                                  ENDIF
                                  CALL SOLVER_MATRIX_INTERFACE_MATRIX_ADD(SOLVER_MATRIX,interface_condition_idx, &
                                    & MatrixCoefficients,INTERFACE_MATRIX,Err,Error,*999)
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

                        IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                          IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                        ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN 
                          IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                        ELSE
                          CALL FlagError("Dynamic solver solve type is not associated.",ERR,ERROR,*999)
                        END IF

                      ENDIF !Update matrix
                    ELSE
                      CALL FlagError("Solver matrix is not associated.",ERR,ERROR,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Invalid number of solver matrices.",ERR,ERROR,*999)
                  ENDIF
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
              ENDIF

              NULLIFY(SOLVER_RHS_VECTOR)
              IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_LINEAR_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_RHS_ONLY) THEN
                IF(DYNAMIC_SOLVER%SOLVER_INITIALISED.OR.(.NOT.DYNAMIC_SOLVER%SOLVER_INITIALISED.AND. &
                  & ((DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_FIRST_ORDER.AND.DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
                  & (DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_SECOND_ORDER.AND.DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_SECOND_DEGREE)))) &
                  & THEN
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
                      !Get the solver variables data                  
                      NULLIFY(CHECK_DATA)
                      CALL DistributedVector_DataGet(SOLVER_RHS_VECTOR,CHECK_DATA,ERR,ERROR,*999)             
                      !Loop over the equations sets
                      DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                        EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                        IF(ASSOCIATED(EQUATIONS_SET)) THEN
                          NULLIFY(DEPENDENT_FIELD) 
                          DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          EQUATIONS=>EQUATIONS_SET%EQUATIONS
                          IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                            IF(ASSOCIATED(EQUATIONS)) THEN
                              NULLIFY(vectorEquations)
                              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                              vectorMatrices=>vectorEquations%vectorMatrices
                              IF(ASSOCIATED(vectorMatrices)) THEN
                                vectorMapping=>vectorEquations%vectorMapping
                                IF(ASSOCIATED(vectorMapping)) THEN

                                  dynamicMapping=>vectorMapping%dynamicMapping
                                  IF(ASSOCIATED(dynamicMapping)) THEN
                                    DYNAMIC_VARIABLE_TYPE=dynamicMapping%dynamicVariableType
                                    !Calculate the dynamic contributions
                                    DYNAMIC_VARIABLE=>dynamicMapping%dynamicVariable
                                    IF(ASSOCIATED(DYNAMIC_VARIABLE)) THEN
                                      dynamicMatrices=>vectorMatrices%dynamicMatrices
                                      IF(ASSOCIATED(dynamicMatrices)) THEN
                                        DYNAMIC_TEMP_VECTOR=>dynamicMatrices%tempVector
                                        !Initialise the dynamic temporary vector to zero
                                        CALL DistributedVector_AllValuesSet(DYNAMIC_TEMP_VECTOR,0.0_DP,ERR,ERROR,*999)
                                        IF(dynamicMapping%stiffnessMatrixNumber/=0) THEN
                                          STIFFNESS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%stiffnessMatrixNumber)%ptr
                                          IF(ASSOCIATED(STIFFNESS_MATRIX)) THEN
                                            NULLIFY(PREDICTED_MEAN_DISPLACEMENT_VECTOR)
                                            CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,PREDICTED_MEAN_DISPLACEMENT_VECTOR, &
                                              & ERR,ERROR,*999)
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                              & -1.0_DP,STIFFNESS_MATRIX%MATRIX, &
!                                              & -DYNAMIC_SOLVER%THETA(1),STIFFNESS_MATRIX%MATRIX, &
                                              & PREDICTED_MEAN_DISPLACEMENT_VECTOR,DYNAMIC_TEMP_VECTOR,ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic stiffness matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                        IF(dynamicMapping%dampingMatrixNumber/=0.AND. &
                                          & DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
                                          DAMPING_MATRIX=>dynamicMatrices%matrices(dynamicMapping%dampingMatrixNumber)%ptr
                                          IF(ASSOCIATED(DAMPING_MATRIX)) THEN
                                            NULLIFY(PREDICTED_MEAN_VELOCITY_VECTOR)
                                            CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,PREDICTED_MEAN_VELOCITY_VECTOR, &
                                              & ERR,ERROR,*999)
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                              & -1.0_DP,DAMPING_MATRIX%MATRIX,PREDICTED_MEAN_VELOCITY_VECTOR,DYNAMIC_TEMP_VECTOR, &
                                              & ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic damping matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                        IF(dynamicMapping%massMatrixNumber/=0.AND. &
                                          & DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                                          MASS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%massMatrixNumber)%ptr
                                          IF(ASSOCIATED(MASS_MATRIX)) THEN
                                            NULLIFY(PREDICTED_MEAN_ACCELERATION_VECTOR)
                                            CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              & FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE,PREDICTED_MEAN_ACCELERATION_VECTOR, &
                                              & ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic mass matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                      ELSE
                                        CALL FlagError("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Dynamic variable is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ELSE
                                    nonlinearMapping=>vectorMapping%nonlinearMapping
                                    IF(ASSOCIATED(nonlinearMapping)) THEN
                                      !Default to FIELD_U_VARIABLE_TYPE
                                      DYNAMIC_VARIABLE_TYPE=FIELD_U_VARIABLE_TYPE
                                      IF(ASSOCIATED(DYNAMIC_TEMP_VECTOR)) NULLIFY(DYNAMIC_TEMP_VECTOR)
                                    ELSE
                                      CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                    !CALL FlagError("Equations mapping dynamic mapping is not associated.",ERR,ERROR,*999)
                                  ENDIF
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
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                              & 1.0_DP,LINEAR_MATRIX%MATRIX,DEPENDENT_VECTOR,LINEAR_TEMP_VECTOR,ERR,ERROR,*999)
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
                                    CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,RHS_VARIABLE_TYPE, &
                                      & FIELD_VALUES_SET_TYPE,RHS_PARAMETERS,ERR,ERROR,*999)
                                    rhsVector=>vectorMatrices%rhsVector
                                    IF(ASSOCIATED(rhsVector)) THEN
                                      BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
  !!TODO: what if the equations set doesn't have a RHS vector???
                                        rhs_variable_type=rhsMapping%rhsVariableType
                                        RHS_VARIABLE=>rhsMapping%rhsVariable
                                        RHS_DOMAIN_MAPPING=>RHS_VARIABLE%DOMAIN_MAPPING
                                        CALL FIELD_PARAMETER_SET_CREATED(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                          & FIELD_INTEGRATED_NEUMANN_SET_TYPE,HAS_INTEGRATED_VALUES,ERR,ERROR,*999)
                                        CALL FIELD_PARAMETER_SET_CREATED(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                          & FIELD_INTEGRATED_ROBIN_SET_TYPE,R_HAS_INTEGRATED_VALUES,ERR,ERROR,*999) !Elias
                                        EQUATIONS_RHS_VECTOR=>rhsVector%VECTOR
                                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,RHS_VARIABLE, &
                                          & RHS_BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                                        IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN
                                          !Update RHS field by integrating any point Neumann or Robin conditions
                                          CALL BoundaryConditions_NeumannIntegrate(RHS_BOUNDARY_CONDITIONS, &
                                            & ERR,ERROR,*999)

                                          CALL BoundaryConditions_RobinIntegrate(RHS_BOUNDARY_CONDITIONS, &
                                            & ERR,ERROR,*999) !Elias

                                          !Loop over the rows in the equations set
                                          DO equations_row_number=1,vectorMapping%totalNumberOfRows
                                            !Get the dynamic contribution to the RHS values
                                          !
                                            IF(ASSOCIATED(DYNAMIC_TEMP_VECTOR)) THEN
                                              CALL DistributedVector_ValuesGet(DYNAMIC_TEMP_VECTOR,equations_row_number, &
                                                & DYNAMIC_VALUE,ERR,ERROR,*999)
                                            ELSE
                                              DYNAMIC_VALUE=0.0_DP
                                            ENDIF
                                          !
                                            !Get the linear matrices contribution to the RHS values if there are any
                                            IF(ASSOCIATED(linearMapping)) THEN
                                              LINEAR_VALUE_SUM=0.0_DP
                                              DO equations_matrix_idx=1,linearMatrices%numberOfLinearMatrices
                                                LINEAR_MATRIX=>linearMatrices%matrices(equations_matrix_idx)%ptr
                                                LINEAR_TEMP_VECTOR=>LINEAR_MATRIX%tempVector
                                                CALL DistributedVector_ValuesGet(LINEAR_TEMP_VECTOR,equations_row_number, &
                                                  & LINEAR_VALUE,ERR,ERROR,*999)
                                                LINEAR_VALUE_SUM=LINEAR_VALUE_SUM+LINEAR_VALUE
                                              ENDDO !equations_matrix_idx
                                              DYNAMIC_VALUE=DYNAMIC_VALUE+LINEAR_VALUE_SUM
                                            ENDIF
                                            !Get the source vector contribute to the RHS values if there are any
                                            IF(ASSOCIATED(sourceMapping)) THEN
                                              !Add in equations source values
                                              CALL DistributedVector_ValuesGet(DISTRIBUTED_SOURCE_VECTOR,equations_row_number, &
                                                & SOURCE_VALUE,ERR,ERROR,*999)
                                              DYNAMIC_VALUE=DYNAMIC_VALUE+SOURCE_VALUE
                                            ENDIF
                                            !Get the nonlinear vector contribute to the RHS values if nonlinear solve
                                            IF(.NOT.STABILITY_TEST) THEN 
                                              IF(SOLVER%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN 
                                                nonlinearMapping=>vectorMapping%nonlinearMapping
                                                  IF(ASSOCIATED(nonlinearMapping)) THEN
                                                   NULLIFY(PREVIOUS_RESIDUAL_PARAMETERS)
                                                   CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                                     & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,PREVIOUS_RESIDUAL_PARAMETERS,ERR,ERROR, &
                                                     & *999)  
                                                   residual_variable_dof=nonlinearMapping% & 
                                                     & equationsRowToResidualDOFMap(equations_row_number)
                                                   PREVIOUS_RESIDUAL_VALUE=-1.0_DP*PREVIOUS_RESIDUAL_PARAMETERS & 
                                                     & (residual_variable_dof)
                                                   DYNAMIC_VALUE=DYNAMIC_VALUE+PREVIOUS_RESIDUAL_VALUE*(1.0_DP-DYNAMIC_SOLVER% & 
                                                     & THETA(1))
                                                  ENDIF
                                              END IF
                                            END IF
                                            !Loop over the solver rows associated with this equations set row
                                            DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                              & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                              solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                                & solver_row_idx)
                                              row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                & COUPLING_COEFFICIENTS(solver_row_idx)
                                               VALUE=DYNAMIC_VALUE*row_coupling_coefficient
                                               CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                                & ERR,ERROR,*999)
                                            ENDDO !solver_row_idx
                                          ENDDO !equations_row_number
                                          
                                          SELECT CASE(DYNAMIC_SOLVER%DEGREE)
                                          CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                                            NULLIFY(FIELD_VALUES_VECTOR)
                                            NULLIFY(PREVIOUS_VALUES_VECTOR)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_VALUES_SET_TYPE,FIELD_VALUES_VECTOR,ERR,ERROR,*999)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_PREVIOUS_VALUES_SET_TYPE,PREVIOUS_VALUES_VECTOR,ERR,ERROR,*999)
                                          CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                                            NULLIFY(FIELD_VALUES_VECTOR)
                                            NULLIFY(PREVIOUS_VALUES_VECTOR)
                                            NULLIFY(PREVIOUS_VELOCITY_VECTOR)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_VALUES_SET_TYPE,FIELD_VALUES_VECTOR,ERR,ERROR,*999)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_PREVIOUS_VALUES_SET_TYPE,PREVIOUS_VALUES_VECTOR,ERR,ERROR,*999)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_PREVIOUS_VELOCITY_SET_TYPE,PREVIOUS_VELOCITY_VECTOR,ERR,ERROR,*999)
                                          CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                                            NULLIFY(FIELD_VALUES_VECTOR)
                                            NULLIFY(PREVIOUS_VALUES_VECTOR)
                                            NULLIFY(PREVIOUS_VELOCITY_VECTOR)
                                            NULLIFY(PREVIOUS_ACCELERATION_VECTOR)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_VALUES_SET_TYPE,FIELD_VALUES_VECTOR,ERR,ERROR,*999)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_PREVIOUS_VALUES_SET_TYPE,PREVIOUS_VALUES_VECTOR,ERR,ERROR,*999)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_PREVIOUS_VELOCITY_SET_TYPE,PREVIOUS_VELOCITY_VECTOR,ERR,ERROR,*999)
                                            CALL FIELD_PARAMETER_SET_DATA_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                              FIELD_PREVIOUS_ACCELERATION_SET_TYPE,PREVIOUS_ACCELERATION_VECTOR,ERR,ERROR,*999)
                                          CASE DEFAULT
                                            LOCAL_ERROR="The dynamic solver degree of "// &
                                              & TRIM(NumberToVString(DYNAMIC_SOLVER%DEGREE,"*",ERR,ERROR))// &
                                              & " is invalid."
                                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                          END SELECT

                                          DO equations_row_number=1,vectorMapping%totalNumberOfRows
                                            !Get the dynamic contribution to the the RHS values
                                            rhs_variable_dof=rhsMapping%equationsRowToRHSDofMap(equations_row_number)
                                            rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                                            rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%DOF_TYPES(rhs_global_dof)
                                            !Apply boundary conditions
                                            SELECT CASE(rhs_boundary_condition)
                                            CASE(BOUNDARY_CONDITION_DOF_FREE)
                                              !Get the equations RHS values
                                              CALL DistributedVector_ValuesGet(EQUATIONS_RHS_VECTOR,equations_row_number, &
                                                & RHS_VALUE,ERR,ERROR,*999)
                                              IF(HAS_INTEGRATED_VALUES) THEN
                                                !Add any Neumann integrated values, b = f + N q+R q_h
                                                CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                                  & FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                                  & ERR,ERROR,*999)
                                                RHS_VALUE=RHS_VALUE+RHS_INTEGRATED_VALUE
                                              END IF

                                              IF(R_HAS_INTEGRATED_VALUES) THEN !Elias /*
                                                !Add any Robin integrated values, b = f + N q+R q_h
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
                                              !Note: the Dirichlet boundary conditions are implicitly included by doing a matrix
                                              !vector product above with the dynamic stiffness matrix and the mean predicited
                                              !displacement vector
                                              !
                                              !This is only true for nonlinear cases and linear cases with fixed values at the boundaries
                                              !
                                              !For changing linear boundary conditions the following needs to be added 
                                              !             
                                              IF(DYNAMIC_SOLVER%UPDATE_BC)THEN
                                                !Set Dirichlet boundary conditions
                                                IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                                                !for linear case only |
!                                                 IF(ASSOCIATED(linearMapping).AND..NOT.ASSOCIATED(nonlinearMapping)) THEN
                                                  !Loop over the dependent variables associated with this equations set row
!                                                   DO variable_idx=1,dynamicMapping%numberOfLinearMatrixVariables
                                                    variable_idx=1
!                                                     variable_type=dynamicMapping%dynamicVariableTypeS(variable_idx)
                                                    variable_type=dynamicMapping%dynamicVariableType
                                                    DEPENDENT_VARIABLE=>dynamicMapping%varToEquationsMatricesMaps( &
                                                      & variable_type)%VARIABLE
                                                    DEPENDENT_VARIABLE_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
                                                    VARIABLE_DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                                                    CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE, &
                                                      & DEPENDENT_BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                                                    variable_dof=dynamicMapping%equationsRowToVariableDOFMaps( &
                                                      & equations_row_number)
                                                    variable_global_dof=VARIABLE_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(variable_dof)
                                                    variable_boundary_condition=DEPENDENT_BOUNDARY_CONDITIONS%DOF_TYPES( &
                                                      & variable_global_dof)

                                                    IF(variable_boundary_condition==BOUNDARY_CONDITION_DOF_FIXED) THEN
                                                      SELECT CASE(DYNAMIC_SOLVER%DEGREE)
                                                      CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                                                        ALPHA_VALUE=(FIELD_VALUES_VECTOR(variable_dof)- &
                                                          & PREVIOUS_VALUES_VECTOR(variable_dof))/ &
                                                          & DYNAMIC_DISPLACEMENT_FACTOR
                                                      CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                                                        ALPHA_VALUE=(FIELD_VALUES_VECTOR(variable_dof)- &
                                                          & PREVIOUS_VALUES_VECTOR(variable_dof)- &
                                                          & DYNAMIC_DISPLACEMENT_FACTOR*PREVIOUS_VELOCITY_VECTOR(variable_dof))/ &
                                                          & DYNAMIC_VELOCITY_FACTOR
                                                      CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                                                        ALPHA_VALUE=(FIELD_VALUES_VECTOR(variable_dof)- &
                                                          & PREVIOUS_VALUES_VECTOR(variable_dof)- &
                                                          & DYNAMIC_DISPLACEMENT_FACTOR*PREVIOUS_VELOCITY_VECTOR(variable_dof) - &
                                                          & DYNAMIC_VELOCITY_FACTOR*PREVIOUS_ACCELERATION_VECTOR(variable_dof))/ &
                                                          & DYNAMIC_ACCELERATION_FACTOR
                                                       CASE DEFAULT
                                                        LOCAL_ERROR="The dynamic solver degree of "// &
                                                          & TRIM(NumberToVString(DYNAMIC_SOLVER%DEGREE,"*",ERR,ERROR))// &
                                                          & " is invalid."
                                                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                                      END SELECT
                                                         

                                                      IF(ABS(ALPHA_VALUE)>=ZERO_TOLERANCE) THEN
                                                        DO equations_matrix_idx=1,dynamicMapping%varToEquationsMatricesMaps( &
                                                          & variable_type)%numberOfEquationsMatrices
                                                          equations_matrix_number=dynamicMapping%varToEquationsMatricesMaps( &
                                                            & variable_type)%equationsMatrixNumbers(equations_matrix_idx)
                                                          IF(equations_matrix_number==dynamicMapping%stiffnessMatrixNumber) & 
                                                            & THEN 
                                                             ALPHA_VALUE=ALPHA_VALUE*STIFFNESS_MATRIX_COEFFICIENT
                                                          ENDIF
                                                          IF(equations_matrix_number==dynamicMapping%dampingMatrixNumber) &
                                                            & THEN 
                                                             ALPHA_VALUE=ALPHA_VALUE*DAMPING_MATRIX_COEFFICIENT
                                                          ENDIF
                                                          IF(equations_matrix_number==dynamicMapping%massMatrixNumber) &
                                                            & THEN 
                                                             ALPHA_VALUE=ALPHA_VALUE*MASS_MATRIX_COEFFICIENT
                                                          ENDIF
                                                          equationsMatrix=>dynamicMatrices% &
                                                            & MATRICES(equations_matrix_number)%ptr
                                                          equations_column_number=dynamicMapping% &
                                                            & varToEquationsMatricesMaps(variable_type)% &
                                                            & dofToColumnsMaps(equations_matrix_idx)% &
                                                            & columnDOF(variable_dof)
                                                          IF(ASSOCIATED(DEPENDENT_BOUNDARY_CONDITIONS% &
                                                            & DIRICHLET_BOUNDARY_CONDITIONS)) THEN
                                                            IF(DEPENDENT_BOUNDARY_CONDITIONS% &
                                                              & NUMBER_OF_DIRICHLET_CONDITIONS>0) THEN
                                                              DO dirichlet_idx=1,DEPENDENT_BOUNDARY_CONDITIONS% &
                                                                & NUMBER_OF_DIRICHLET_CONDITIONS
                                                                IF(DEPENDENT_BOUNDARY_CONDITIONS% &
                                                                  & DIRICHLET_BOUNDARY_CONDITIONS% &
                                                                  & DIRICHLET_DOF_INDICES(dirichlet_idx)== &
                                                                  & equations_column_number) EXIT
                                                              ENDDO
                                                              SELECT CASE(equationsMatrix%storageType)
                                                              CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                                                DO dirichlet_row=1,vectorMatrices%totalNumberOfRows
                                                                  CALL DistributedMatrix_ValuesGet(equationsMatrix% &
                                                                    & MATRIX,dirichlet_row,equations_column_number, &
                                                                    & MATRIX_VALUE,ERR,ERROR,*999)
                                                                  IF(ABS(MATRIX_VALUE)>=ZERO_TOLERANCE) THEN
                                                                    DO solver_row_idx=1,SOLVER_MAPPING% &
                                                                      & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% & 
                                                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                      & dirichlet_row)%NUMBER_OF_SOLVER_ROWS
                                                                      solver_row_number=SOLVER_MAPPING% &
                                                                        & EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                        & equations_set_idx)% &
                                                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                        & dirichlet_row)%SOLVER_ROWS(solver_row_idx)
                                                                      row_coupling_coefficient=SOLVER_MAPPING% &
                                                                        & EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                        & equations_set_idx)% &
                                                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                        & dirichlet_row)%COUPLING_COEFFICIENTS( &
                                                                        & solver_row_idx)
                                                                      VALUE=-1.0_DP*MATRIX_VALUE*ALPHA_VALUE* &
                                                                        & row_coupling_coefficient
                                                                      CALL DistributedVector_ValuesAdd( &
                                                                        & SOLVER_RHS_VECTOR, &
                                                                        & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                    ENDDO !solver_row_idx
                                                                  ENDIF
                                                                ENDDO !dirichlet_row
                                                              CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                                                dirichlet_row=equations_column_number
                                                                CALL DistributedMatrix_ValuesGet(equationsMatrix% &
                                                                  & MATRIX,dirichlet_row,equations_column_number, &
                                                                  & MATRIX_VALUE,ERR,ERROR,*999)
                                                                IF(ABS(MATRIX_VALUE)>=ZERO_TOLERANCE) THEN
                                                                  DO solver_row_idx=1,SOLVER_MAPPING% &
                                                                    & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% & 
                                                                    & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                    & dirichlet_row)%NUMBER_OF_SOLVER_ROWS
                                                                    solver_row_number=SOLVER_MAPPING% &
                                                                      & EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                      & equations_set_idx)% &
                                                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                      & dirichlet_row)%SOLVER_ROWS(solver_row_idx)
                                                                    row_coupling_coefficient=SOLVER_MAPPING% &
                                                                      & EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                      & equations_set_idx)% &
                                                                      & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                      & dirichlet_row)%COUPLING_COEFFICIENTS( &
                                                                      & solver_row_idx)
                                                                    VALUE=-1.0_DP*MATRIX_VALUE*ALPHA_VALUE* &
                                                                      & row_coupling_coefficient
                                                                    CALL DistributedVector_ValuesAdd( &
                                                                      & SOLVER_RHS_VECTOR, &
                                                                      & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                  ENDDO !solver_row_idx
                                                                ENDIF

                                                              CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                                                CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                              CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                                                CALL FlagError("Not implemented.",ERR,ERROR,*999)
                                                              CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                                                SPARSITY_INDICES=>DEPENDENT_BOUNDARY_CONDITIONS% &
                                                                  & DIRICHLET_BOUNDARY_CONDITIONS%DYNAMIC_SPARSITY_INDICES( &
                                                                  & equations_set_idx,equations_matrix_idx)%ptr
                                                                IF(ASSOCIATED(SPARSITY_INDICES)) THEN
                                                                  DO equations_row_number2=SPARSITY_INDICES% & 
                                                                    & SPARSE_COLUMN_INDICES(dirichlet_idx), &
                                                                    & SPARSITY_INDICES%SPARSE_COLUMN_INDICES( &
                                                                    & dirichlet_idx+1)-1
                                                                    dirichlet_row=SPARSITY_INDICES%SPARSE_ROW_INDICES( &
                                                                      & equations_row_number2)
                                                                    CALL DistributedMatrix_ValuesGet(equationsMatrix% &
                                                                      & MATRIX,dirichlet_row,equations_column_number, &
                                                                      & MATRIX_VALUE,ERR,ERROR,*999)
                                                                    IF(ABS(MATRIX_VALUE)>=ZERO_TOLERANCE) THEN
                                                                      DO solver_row_idx=1,SOLVER_MAPPING% &
                                                                        & EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% & 
                                                                        & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                        & dirichlet_row)%NUMBER_OF_SOLVER_ROWS
                                                                        solver_row_number=SOLVER_MAPPING% &
                                                                          & EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                          & equations_set_idx)% &
                                                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                          & dirichlet_row)%SOLVER_ROWS(solver_row_idx)
                                                                        row_coupling_coefficient=SOLVER_MAPPING% &
                                                                          & EQUATIONS_SET_TO_SOLVER_MAP( &
                                                                          & equations_set_idx)% &
                                                                          & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS( &
                                                                          & dirichlet_row)%COUPLING_COEFFICIENTS( &
                                                                          & solver_row_idx)
                                                                        VALUE=-1.0_DP*MATRIX_VALUE*ALPHA_VALUE* &
                                                                          & row_coupling_coefficient
                                                                        CALL DistributedVector_ValuesAdd( &
                                                                          & SOLVER_RHS_VECTOR, &
                                                                          & solver_row_number,VALUE,ERR,ERROR,*999)
                                                                      ENDDO !solver_row_idx
                                                                    ENDIF
                                                                  ENDDO !equations_row_number2
                                                                ELSE
                                                                  CALL FlagError("Sparsity indices are not associated.", &
                                                                    & ERR,ERROR,*999)
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
                                                            CALL FlagError("Dirichlet boundary conditions is &
                                                              & not associated.",ERR,ERROR,*999)
                                                          ENDIF
                                                        ENDDO !matrix_idx
                                                      ENDIF
                                                    ENDIF
!                                                   ENDDO !variable_idx
                                                ENDIF
                                              ENDIF

                                            CASE(BOUNDARY_CONDITION_DOF_FIXED)
                                              !Set Neumann and Robin boundary conditions
                                              !Loop over the solver rows associated with this equations set row
                                              DO solver_row_idx=1,SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%NUMBER_OF_SOLVER_ROWS
                                                solver_row_number=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(equations_set_idx)% &
                                                  & EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)%SOLVER_ROWS( &
                                                  & solver_row_idx)
                                                row_coupling_coefficient=SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP( &
                                                  & equations_set_idx)%EQUATIONS_ROW_TO_SOLVER_ROWS_MAPS(equations_row_number)% &
                                                  & COUPLING_COEFFICIENTS(solver_row_idx)
                                                VALUE=RHS_PARAMETERS(rhs_variable_dof)*row_coupling_coefficient
                                                IF(HAS_INTEGRATED_VALUES) THEN
                                                  !Add any Neumann integrated values, b = f + N q + R q_h
                                                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                                    & FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                                    & ERR,ERROR,*999)
                                                  VALUE=VALUE+RHS_INTEGRATED_VALUE*row_coupling_coefficient
                                                END IF
                                                IF(R_HAS_INTEGRATED_VALUES) THEN
                                                  !Add any Robin integrated values, b = f + N q + R q_h
                                                  CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(RHS_VARIABLE%FIELD,RHS_VARIABLE_TYPE, &
                                                    & FIELD_INTEGRATED_ROBIN_SET_TYPE,rhs_variable_dof,RHS_INTEGRATED_VALUE, &
                                                    & ERR,ERROR,*999)
                                                  VALUE=VALUE+RHS_INTEGRATED_VALUE*row_coupling_coefficient
                                                END IF
                                                CALL DistributedVector_ValuesAdd(SOLVER_RHS_VECTOR,solver_row_number,VALUE, &
                                                  & ERR,ERROR,*999)
                                              ENDDO !solver_row_idx
                                            CASE(BOUNDARY_CONDITION_DOF_MIXED)
                                              !Set Robin or is it Cauchy??? boundary conditions.
                                              CALL FlagError("Mixed Boundary Conditions Not implemented.",ERR,ERROR,*999)
                                            CASE DEFAULT
                                              LOCAL_ERROR="The RHS boundary condition of "// &
                                                & TRIM(NumberToVString(rhs_boundary_condition,"*",ERR,ERROR))// &
                                                & " for RHS variable dof number "// &
                                                & TRIM(NumberToVString(rhs_variable_dof,"*",ERR,ERROR))//" is invalid."
                                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                                            END SELECT
                                          ENDDO !equations_row_number
                                        ELSE
                                          CALL FlagError("RHS boundary conditions variable is not associated.",ERR,ERROR,*999)
                                        ENDIF
                                      ELSE
                                        CALL FlagError("Equations set boundary conditions is not associated.",ERR,ERROR,*999)
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
                            CALL FlagError("Equations set dependent field is not associated.",ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
                        ENDIF
                      ENDDO !equations_set_idx
              !!!! TODO TODO !!!! ???
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
              !        
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
              END IF

              NULLIFY(SOLVER_RESIDUAL_VECTOR)
              IF(SELECTION_TYPE==SOLVER_MATRICES_ALL.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
                & SELECTION_TYPE==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
                IF(DYNAMIC_SOLVER%SOLVER_INITIALISED.OR.(.NOT.DYNAMIC_SOLVER%SOLVER_INITIALISED.AND. &
                  & ((DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_FIRST_ORDER.AND.DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_FIRST_DEGREE).OR. & 
                  & (DYNAMIC_SOLVER%ORDER==SOLVER_DYNAMIC_SECOND_ORDER.AND.DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_SECOND_DEGREE)))) &
                  & THEN
                  !Assemble residual vector
                  IF(SOLVER%outputType>=SOLVER_TIMING_OUTPUT) THEN
                    CALL CPUTimer(USER_CPU,USER_TIME1,ERR,ERROR,*999)
                    CALL CPUTimer(SYSTEM_CPU,SYSTEM_TIME1,ERR,ERROR,*999)
                  ENDIF
                  IF(SOLVER_MATRICES%UPDATE_RESIDUAL) THEN
                    SOLVER_RESIDUAL_VECTOR=>SOLVER_MATRICES%RESIDUAL
                    IF(ASSOCIATED(SOLVER_RESIDUAL_VECTOR)) THEN
                      !Initialise the residual to zero
                      CALL DistributedVector_AllValuesSet(SOLVER_RESIDUAL_VECTOR,0.0_DP,ERR,ERROR,*999)
                      !Get the solver variables data
                      NULLIFY(CHECK_DATA)
                      CALL DistributedVector_DataGet(SOLVER_RESIDUAL_VECTOR,CHECK_DATA,ERR,ERROR,*999)
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
                                  dynamicMapping=>vectorMapping%dynamicMapping
                                  IF(ASSOCIATED(dynamicMapping)) THEN
                                    DYNAMIC_VARIABLE_TYPE=dynamicMapping%dynamicVariableType
                                    !Calculate the dynamic contributions
                                    DYNAMIC_VARIABLE=>dynamicMapping%dynamicVariable
                                    IF(ASSOCIATED(DYNAMIC_VARIABLE)) THEN
                                      dynamicMatrices=>vectorMatrices%dynamicMatrices
                                      IF(ASSOCIATED(dynamicMatrices)) THEN
                                        DYNAMIC_TEMP_VECTOR=>dynamicMatrices%tempVector
                                        !Initialise the dynamic temporary vector to zero
                                        CALL DistributedVector_AllValuesSet(DYNAMIC_TEMP_VECTOR,0.0_DP,ERR,ERROR,*999)
                                        NULLIFY(INCREMENTAL_VECTOR)
                                        !Define the pointer to the INCREMENTAL_VECTOR
                                        CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE, &
                                          & FIELD_INCREMENTAL_VALUES_SET_TYPE,INCREMENTAL_VECTOR,ERR,ERROR,*999)
                                        IF(dynamicMapping%stiffnessMatrixNumber/=0) THEN
                                          STIFFNESS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%stiffnessMatrixNumber)%ptr
                                          IF(ASSOCIATED(STIFFNESS_MATRIX)) THEN
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, & 
                                              & STIFFNESS_MATRIX_COEFFICIENT,STIFFNESS_MATRIX%MATRIX,INCREMENTAL_VECTOR, & 
                                              & DYNAMIC_TEMP_VECTOR,ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic stiffness matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                        IF(dynamicMapping%dampingMatrixNumber/=0.AND. &
                                          & DYNAMIC_SOLVER%DEGREE>=SOLVER_DYNAMIC_FIRST_DEGREE) THEN
                                          DAMPING_MATRIX=>dynamicMatrices%matrices(dynamicMapping%dampingMatrixNumber)%ptr
                                          IF(ASSOCIATED(DAMPING_MATRIX)) THEN
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                              & DAMPING_MATRIX_COEFFICIENT,DAMPING_MATRIX%MATRIX,INCREMENTAL_VECTOR, & 
                                              & DYNAMIC_TEMP_VECTOR,ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic damping matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                        IF(dynamicMapping%massMatrixNumber/=0.AND. &
                                          & DYNAMIC_SOLVER%DEGREE>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                                          MASS_MATRIX=>dynamicMatrices%matrices(dynamicMapping%massMatrixNumber)%ptr
                                          IF(ASSOCIATED(MASS_MATRIX)) THEN
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                              & MASS_MATRIX_COEFFICIENT,MASS_MATRIX%MATRIX,INCREMENTAL_VECTOR, & 
                                              & DYNAMIC_TEMP_VECTOR,ERR,ERROR,*999)
                                          ELSE
                                            CALL FlagError("Dynamic mass matrix is not associated.",ERR,ERROR,*999)
                                          ENDIF
                                        ENDIF
                                      ELSE
                                        CALL FlagError("Dynamic variable is not associated.",ERR,ERROR,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Equations matrices dynamic matrices is not associated.",ERR,ERROR,*999)
                                    ENDIF
                                  ENDIF
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
                                            CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                              & 1.0_DP,LINEAR_MATRIX%MATRIX,DEPENDENT_VECTOR,LINEAR_TEMP_VECTOR,ERR,ERROR,*999)
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
                                          IF(STABILITY_TEST) THEN
                                            RESIDUAL_VALUE=RESIDUAL_VALUE
                                          ELSE
                                            RESIDUAL_VALUE=RESIDUAL_VALUE*DYNAMIC_SOLVER%THETA(1)
                                          ENDIF
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
                                          IF(ASSOCIATED(dynamicMapping)) THEN
                                            !Get the dynamic contribution to the residual values
                                            CALL DistributedVector_ValuesGet(DYNAMIC_TEMP_VECTOR,equations_row_number, &
                                              & DYNAMIC_VALUE,ERR,ERROR,*999)
                                               RESIDUAL_VALUE=RESIDUAL_VALUE+DYNAMIC_VALUE
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
                                  ELSE
                                    CALL FlagError("Equations mapping nonlinear mapping is not associated.",ERR,ERROR,*999)
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
                                        
                                    !
                                        SELECT CASE(INTERFACE_MATRIX%INTERFACE_MATRIX_TIME_DEPENDENCE_TYPE)
                                        CASE(INTERFACE_MATRIX_STATIC)
                                          MatrixCoefficients(1)=STIFFNESS_MATRIX_COEFFICIENT
                                        CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                                          MatrixCoefficients(1)=DAMPING_MATRIX_COEFFICIENT
                                        CASE DEFAULT
                                          CALL FlagError("Not implemented.",Err,Error,*999)
                                        END SELECT
                                        IF(INTERFACE_MATRIX%HAS_TRANSPOSE) THEN
                                          SELECT CASE(INTERFACE_MATRIX%INTERFACE_MATRIX_TRANSPOSE_TIME_DEPENDENCE_TYPE)
                                          CASE(INTERFACE_MATRIX_STATIC)
                                            MatrixCoefficients(2)=STIFFNESS_MATRIX_COEFFICIENT
                                          CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                                            MatrixCoefficients(2)=DAMPING_MATRIX_COEFFICIENT
                                          CASE DEFAULT

                                            CALL FlagError("Not implemented.",Err,Error,*999)
                                          END SELECT
                                        ELSE
                                          MatrixCoefficients(2)=0.0_DP
                                        ENDIF
                                    !
                                        
                                        
                                       ! CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                                       !   & INTERFACE_MATRIX%MATRIX,LAGRANGE_VECTOR,INTERFACE_TEMP_VECTOR,ERR,ERROR,*999)
                                        CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                          & MatrixCoefficients(1),INTERFACE_MATRIX%MATRIX,LAGRANGE_VECTOR,INTERFACE_TEMP_VECTOR, &
                                          & ERR,ERROR,*999)
                                        
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
                                        !hard-coded for now TODO under the assumption that the first equations set is the solid
                                        !equations set and the second equations set is the fluid equations set
                                        !FSI only - needs to be extended/generalized for other coupled problems TODO
                                        IF(interface_matrix_idx==1) THEN
                                          CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,dependent_variable_type, &
                                            & FIELD_INCREMENTAL_VALUES_SET_TYPE,DEPENDENT_VECTOR,ERR,ERROR,*999)
                                        ELSE
                                          CALL FIELD_PARAMETER_SET_VECTOR_GET(DEPENDENT_FIELD,dependent_variable_type, &
                                            & FIELD_VALUES_SET_TYPE,DEPENDENT_VECTOR,ERR,ERROR,*999)
                                        ENDIF
                                       ! CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                                       !   & INTERFACE_MATRIX%MATRIX_TRANSPOSE,DEPENDENT_VECTOR,INTERFACE_TEMP_VECTOR,ERR,ERROR,*999)
                                        CALL DistributedMatrix_MatrixByVectorAdd(DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                                          & MatrixCoefficients(2),INTERFACE_MATRIX%MATRIX_TRANSPOSE,DEPENDENT_VECTOR, &
                                          & INTERFACE_TEMP_VECTOR,ERR,ERROR,*999)
                                        
                                        !Add interface matrix residual contribution to the solver residual.
                                        !The number of columns in the interface matrix is equivalent to the number of rows of the transposed interface matrices
                                        DO interface_row_number=1,INTERFACE_MATRICES%NUMBER_OF_COLUMNS
                                          IF(SOLVER_MAPPING%INTERFACE_CONDITION_TO_SOLVER_MAP(interface_condition_idx)% &
                                            & INTERFACE_COLUMN_TO_SOLVER_ROWS_MAPS(interface_row_number)% &
                                            & NUMBER_OF_SOLVER_ROWS>0) THEN
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
                                         !   IF(interface_matrix_idx==1) THEN
                                         !     VALUE=RESIDUAL_VALUE*row_coupling_coefficient/DELTA_T
                                         !   ELSE
                                              VALUE=RESIDUAL_VALUE*row_coupling_coefficient
                                         !   ENDIF
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
                  !
                      !Start the update the solver residual vector values
                      CALL DistributedVector_UpdateStart(SOLVER_RESIDUAL_VECTOR,ERR,ERROR,*999)

                      NULLIFY(CHECK_DATA2)
                      CALL DistributedVector_DataGet(SOLVER_RESIDUAL_VECTOR,CHECK_DATA2,ERR,ERROR,*999)

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
              ENDIF

              IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) THEN
                !Set the first part of the next time step. Note that we do not have to add in the previous time value as it is
                !already there from when we copied the values to the previous time step.
                !Loop over the equations sets
                IF(DYNAMIC_SOLVER%DEGREE>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                        EQUATIONS=>EQUATIONS_SET%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          NULLIFY(vectorEquations)
                          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                          vectorMapping=>vectorEquations%vectorMapping
                          IF(ASSOCIATED(vectorMapping)) THEN
                            dynamicMapping=>vectorMapping%dynamicMapping
                            IF(ASSOCIATED(dynamicMapping)) THEN
                              DYNAMIC_VARIABLE_TYPE=dynamicMapping%dynamicVariableType
                              SELECT CASE(DYNAMIC_SOLVER%DEGREE)
                              CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                                !Do nothing. Increment will be added after the solve.
                              CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                                CALL FIELD_PARAMETER_SETS_ADD(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE,FIRST_UPDATE_FACTOR, &
                                  & FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                              CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                                CALL FIELD_PARAMETER_SETS_ADD(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE,[FIRST_UPDATE_FACTOR, &
                                  & SECOND_UPDATE_FACTOR],[FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE], &
                                  & FIELD_VALUES_SET_TYPE,ERR,ERROR,*999)
                                CALL FIELD_PARAMETER_SETS_ADD(DEPENDENT_FIELD,DYNAMIC_VARIABLE_TYPE,FIRST_UPDATE_FACTOR, &
                                  & FIELD_PREVIOUS_ACCELERATION_SET_TYPE,FIELD_VELOCITY_VALUES_SET_TYPE,ERR,ERROR,*999)
                              CASE DEFAULT
                                LOCAL_ERROR="The dynamic solver degree of "// &
                                  & TRIM(NumberToVString(DYNAMIC_SOLVER%DEGREE,"*",ERR,ERROR))//" is invalid."
                                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                              END SELECT
                            ELSE
                              LOCAL_ERROR="Equations mapping dynamic mapping is not associated for equations set index number "// &
                                & TRIM(NumberToVString(equations_set_idx,"*",ERR,ERROR))//"."
                              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                            ENDIF
                          ELSE
                            LOCAL_ERROR="Equations equations mapping is not associated for equations set index number "// &
                              & TRIM(NumberToVString(equations_set_idx,"*",ERR,ERROR))//"."
                            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                          ENDIF
                        ELSE
                          LOCAL_ERROR="Equations set equations is not associated for equations set index number "// &
                            & TRIM(NumberToVString(equations_set_idx,"*",ERR,ERROR))//"."
                          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                        ENDIF
                      ELSE
                        LOCAL_ERROR="Equations set dependent field is not associated for equations set index number "// &
                          & TRIM(NumberToVString(equations_set_idx,"*",ERR,ERROR))//"."
                        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index number "// &
                        & TRIM(NumberToVString(equations_set_idx,"*",ERR,ERROR))//"."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ENDIF
              ENDIF
              !If required output the solver matrices
              IF(SOLVER%outputType>=SOLVER_MATRIX_OUTPUT) THEN
                CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SELECTION_TYPE,SOLVER_MATRICES,ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Solver solver matrices is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Solver solver equations is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver dynamic solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE")
    RETURN
999 ERRORSEXITS("SOLVER_MATRICES_DYNAMIC_ASSEMBLE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE SOLVER_MATRICES_DYNAMIC_ASSEMBLE
