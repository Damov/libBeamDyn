!..................................................................................................................................
! LICENSING                                                                                                                         
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    Glue is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!  
!    ADD DESCRIPTION
!	
!    References:
!
!
!**********************************************************************************************************************************
MODULE BeamDyn_Library
    use iso_c_binding

   USE BeamDyn_SP
   USE BeamDyn_Types

   USE NWTC_Library

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                     :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                    :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   REAL(DbKi)                         :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                         :: t_initial        ! time at initialization
!   REAL(DbKi)                         :: t_final          ! time at simulation end 
   REAL(DbKi)                         :: t_global         ! global-loop time marker

!   INTEGER(IntKi)                     :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                     :: n_t_global       ! global-loop time counter

   INTEGER(IntKi)                     :: pc_max           ! 1:explicit loose; 2:pc loose
   INTEGER(IntKi)                     :: pc               ! counter for pc iterations

   INTEGER(IntKi)                     :: BD_interp_order     ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InitInputType)           :: BD_InitInput
   TYPE(BD_ParameterType)           :: BD_Parameter
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousStateDeriv
   TYPE(BD_InitOutputType)          :: BD_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   TYPE(BD_OtherStateType)          :: BD_OtherState

   TYPE(BD_InputType),ALLOCATABLE  :: BD_Input(:)
   REAL(DbKi),        ALLOCATABLE  :: BD_InputTimes(:)

   TYPE(BD_OutputType),ALLOCATABLE  :: BD_Output(:)
   REAL(DbKi),ALLOCATABLE           :: BD_OutputTimes(:)

   TYPE(BD_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(BD_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState_pred
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState_pred
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState_pred

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   INTEGER(IntKi),PARAMETER:: QiDisUnit = 20
   INTEGER(IntKi),PARAMETER:: QiRootUnit = 30
   INTEGER(IntKi),PARAMETER:: QiReacUnit = 40

   REAL(ReKi):: temp_H(3,3)
   REAL(ReKi):: temp_cc(3)
   REAL(ReKi):: temp_R(3,3)
   REAL(ReKi):: temp_vec(3)
   REAL(ReKi):: temp1,temp2
   REAL(ReKi):: temp_a,temp_b,temp_force
   REAL(DbKi):: start, finish
   Integer(IntKi)                     :: temp_count
!   REAL(ReKi),PARAMETER    :: PI = 3.1415926D0

   !
   ! OTHER VARS for interfacing BeamDyn and OpenFOAM
   !
   logical, parameter :: writeOutput = .false.
   logical, parameter :: timing = .false.
   logical, parameter :: debug = .false.

   Integer(IntKi) :: order_elem

   real(DbKi), dimension(:,:), allocatable :: F_foam, M_foam
   real(DbKi), dimension(:,:), allocatable :: F_foam_t0, M_foam_t0
   real(DbKi), dimension(:,:), allocatable :: F_foam_t00, M_foam_t00
   integer, dimension(:), allocatable :: foamUpdatedNode

CONTAINS

!   7000 FORMAT (ES12.5,9ES21.12)

!********************************************************************************************!
!
SUBROUTINE init(t0_foam,dt_foam) bind(c, name="beamDynStart")
!
!********************************************************************************************!
    !real(DbKi), intent(in) :: dt_foam
    real(R8Ki), intent(in) :: t0_foam, dt_foam

   write(*,*) 'BEAMDYN: Input t0= ',t0_foam
   write(*,*) 'BEAMDYN: Input dt= ',dt_foam

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   !vvvvvvvv Edit variables below this line vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

!   t_initial = 0.0D+00
    t_initial = t0_foam
!   t_final   = 1.0D+02

!   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
!   dt_global = 5.0D-03

! Let user specify a variable time step at every step() call
   dt_global = dt_foam

!   BD_InitInput%InputFile = 'BeamDyn_Input_5MW_New.inp'
   BD_InitInput%InputFile = 'beam.input'

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BD_interp_order = 2

   !^^^^^^^^ Edit variables above this line ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1)) 
   ALLOCATE(BD_InputTimes(BD_interp_order + 1)) 

   ALLOCATE(BD_Output(BD_interp_order + 1)) 
   ALLOCATE(BD_OutputTimes(BD_interp_order + 1)) 

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------
   if( writeOutput ) then
        OPEN(unit = QiDisUnit, file = 'QiDisp_AM2.out', status = 'REPLACE',ACTION = 'WRITE')
        OPEN(unit = QiRootUnit,file = 'QiRoot_AM2.out', status = 'REPLACE',ACTION = 'WRITE')
        OPEN(unit = QiReacUnit,file = 'QiReac_AM2.out', status = 'REPLACE',ACTION = 'WRITE')
   end if


   BD_InitInput%RootName  = TRIM(BD_Initinput%InputFile)
   ALLOCATE(BD_InitInput%gravity(3)) 
   BD_InitInput%gravity(1) = 0.0D0 !-9.80665
   BD_InitInput%gravity(2) = 0.0D0 
   BD_InitInput%gravity(3) = 0.0D0 

   ! position vector for the blade coordinate system origin
   ALLOCATE(BD_InitInput%GlbPos(3)) 
   BD_InitInput%GlbPos(1) = 0.0D+00
   BD_InitInput%GlbPos(2) = 0.0D+01
   BD_InitInput%GlbPos(3) = 0.0D0

   ! transformation for blade coordinate system (blade is along x1)
   ALLOCATE(BD_InitInput%GlbRot(3,3)) 
   BD_InitInput%GlbRot(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 0.0
   temp_vec(3) = 0.0 !4.0D0*TAN((3.1415926D0/2.0D0)/4.0D0)
   CALL CrvMatrixR(temp_vec,temp_R)
   BD_InitInput%GlbRot(1:3,1:3) = temp_R(1:3,1:3)

   CALL BeamDyn_Init(BD_InitInput        &
                   , BD_Input(1)         &
                   , BD_Parameter        &
                   , BD_ContinuousState  &
                   , BD_DiscreteState    &
                   , BD_ConstraintState  &
                   , BD_OtherState       &
                   , BD_Output(1)        &
                   , dt_global             &
                   , BD_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )


   CALL BD_CopyInput(  BD_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BD_CopyOutput( BD_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(Mod1_Input)
   DO i = 1, BD_interp_order + 1  
      BD_InputTimes(i) = t_initial - (i - 1) * dt_global
      BD_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

   DO i = 1, BD_interp_order
      CALL BD_CopyInput (BD_Input(i),  BD_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      CALL BD_CopyOutput (BD_Output(i),  BD_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   ENDDO

   write(*,*) '********************'
   write(*,*) 'Constant inputs:'
   write(*,*) '  blade length =', BD_Parameter%blade_length
   write(*,*) '  nodes/elem   =', BD_Parameter%node_elem
   write(*,*) '  dof/node     =', BD_Parameter%dof_node
   write(*,*) '  tot # elems  =', BD_Parameter%elem_total
   write(*,*) '  tot # nodes  =', BD_Parameter%node_total
   write(*,*) '  tot # dofs   =', BD_Parameter%dof_total
   write(*,*) '  # gauss pts  =', BD_Parameter%ngp
   write(*,*) '  interp order =', BD_interp_order !max is 2 for integration in TIME domain
   write(*,*) '********************'

   order_elem = BD_Parameter%node_elem-1

   write(*,*) ''
   write(*,*) 'BEAMDYN: Allocating distributed loads arrays'
   allocate(       F_foam(3, order_elem) )
   allocate(       M_foam(3, order_elem) )
   allocate(    F_foam_t0(3, order_elem) )
   allocate(    M_foam_t0(3, order_elem) )
!   allocate(   F_foam_t00(3, order_elem) )
!   allocate(   M_foam_t00(3, order_elem) )
   allocate( foamUpdatedNode(order_elem) )
   F_foam(:,:) = 0.0
   M_foam(:,:) = 0.0
   F_foam_t0(:,:) = 0.0
   M_foam_t0(:,:) = 0.0
!   F_foam_t00(:,:) = 0.0
!   M_foam_t00(:,:) = 0.0
   foamUpdatedNode(:) = 0

!   write(*,*) ''
!   write(*,*) 'BEAMDYN: Calling step() to initialize solution --is this necessary?'
!   call step(dt_foam)

END SUBROUTINE init


!********************************************************************************************!
!
!SUBROUTINE step() bind(c, name="beamDynStep")
SUBROUTINE step(dt_foam) bind(c, name="beamDynStep")
!
!********************************************************************************************!
    real(R8Ki), intent(in) :: dt_foam

    ! SANITY CHECK
    if(n_t_global > 0 .and. .not. all(foamUpdatedNode==1)) then
        write(*,*) 'BEAMDYN: Warning, loads data probably out of sync!', foamUpdatedNode
    end if
    foamUpdatedNode(:) = 0

    if(timing) CALL CPU_TIME(start)

!   DO n_t_global = 0, n_t_final-1

      !WRITE(*,*) "BEAMDYN: Step ", n_t_global, &
      !           " Time-step size= ", dt_foam, &
      !           "  t= ", t_global
      WRITE(*,'(" BEAMDYN:  Step ",i7,"  Time-step size= ",es16.8,"  t0= ",f16.8)') &
            n_t_global, dt_foam, t_global
      if(.not. dt_foam==dt_global) then
          write(*,*) "BEAMDYN: --Note-- Time step has changed!", dt_global, dt_foam
      end if

! original:
!      CALL BD_InputSolve( t_global,                   BD_Input(1), BD_InputTimes(1), BD_Parameter, ErrStat, ErrMsg)
!      CALL BD_InputSolve( t_global + dt_global,       BD_Input(2), BD_InputTimes(2), BD_Parameter, ErrStat, ErrMsg)
!      CALL BD_InputSolve( t_global + 2.0D0*dt_global, BD_Input(3), BD_InputTimes(3), BD_Parameter, ErrStat, ErrMsg)

!                         time                      root motion   "don't worry"      fixed params
!                                                 + applied force 
!                         input                     in/output     in/output          input         input       input       output   output
!     CALL BD_InputSolve( t_global,                 BD_Input(1),  BD_InputTimes(1),  BD_Parameter, F_foam_t00, M_foam_t00, ErrStat, ErrMsg)
!     CALL BD_InputSolve( t_global + dt_foam,       BD_Input(2),  BD_InputTimes(2),  BD_Parameter, F_foam_t0,  M_foam_t0,  ErrStat, ErrMsg)
!     CALL BD_InputSolve( t_global + 2.0D0*dt_foam, BD_Input(3),  BD_InputTimes(3),  BD_Parameter, F_foam,     M_foam,     ErrStat, ErrMsg)

      CALL BD_InputSolve( t_global + dt_foam,       BD_Input(2),  BD_InputTimes(2),  BD_Parameter, F_foam, M_foam, ErrStat, ErrMsg)

!--------------------------------------------
! Compute initial condition given root motion
!--------------------------------------------
      IF(n_t_global .EQ. 0) THEN
          CALL BD_InputSolve( t_global, BD_Input(1), BD_InputTimes(1), BD_Parameter, F_foam_t0, M_foam_t0, ErrStat, ErrMsg)
          CALL BD_InitialCondition(BD_Input(1),BD_Parameter,BD_ContinuousState,ErrStat,ErrMsg)
      ENDIF
!------------------------------
! END Compute initial condition
!------------------------------

     ! updates the BD_Output(1) structure AT TIME t
     !                        input     input        input         in/output           input
     !                        input     in/output    in/output     output              output
!    CALL BeamDyn_CalcOutput( t_global, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
!                            BD_ConstraintState, BD_OtherState,  BD_Output(1), ErrStat, ErrMsg)


     ! time marching from t to t + dt occurs here!
     ! updates BD_ContinuousState, i.e. 'x', which contains the displacements (including orientations) 
     !   and velocities (including angular velocities) AT TIME t + dt
     CALL BeamDyn_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, &
                               BD_ContinuousState, &
                               BD_DiscreteState, BD_ConstraintState, &
                               BD_OtherState, ErrStat, ErrMsg )

     ! copies solution from UPDATED BD_ContinuousState to the BD_Output(2) structure AT TIME t + dt
     !                        input     input        input         in/output           input
     !                        input     in/output    in/output     output              output
     CALL BeamDyn_CalcOutput( t_global, BD_Input(2), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                              BD_ConstraintState, BD_OtherState,  BD_Output(2), ErrStat, ErrMsg)


      if( writeOutput ) then
          6000 FORMAT (ES12.5,6ES21.12)

          CALL CrvExtractCrv(BD_OutPut(2)%BldMotion%Orientation(1:3,1:3,BD_Parameter%node_total),temp_cc)
          WRITE(QiDisUnit,6000) t_global,&
                               &BD_OutPut(2)%BldMotion%TranslationDisp(1:3,BD_Parameter%node_total),&
                               &temp_cc(1:3)

          WRITE(QiRootUnit,6000) t_global,&
                               &BD_OutPut(2)%BldForce%Force(1:3,1),&
                               &BD_OutPut(2)%BldForce%Moment(1:3,1)

          WRITE(QiReacUnit,6000) t_global,&
                               &BD_OutPut(2)%ReactionForce%Force(1:3,1),&
                               &BD_OutPut(2)%ReactionForce%Moment(1:3,1)
      end if

      ! update the global time

!      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial
      t_global = t_global+dt_foam + t_initial

!   ENDDO

    if( timing ) then
        CALL CPU_TIME(finish)
        !print '("Time = ",f6.3," seconds.")',finish-start
        WRITE(*,*) 'BEAMDYN: Start ', start
        WRITE(*,*) 'BEAMDYN: Finish ', finish
        WRITE(*,*) 'BEAMDYN: Time ', finish-start
    end if

   if( debug ) then
   ! calculate final time normalized rms error
    !   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(Mod1_Parameter%method)))
       CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

       !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
       CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

    !  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
    !  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )
    end if

    n_t_global = n_t_global + 1

END SUBROUTINE step


!********************************************************************************************!
!********************************************************************************************!
!
!SUBROUTINE foamUpdateNode(inode0,F,M) bind(c, name="beamDynSetDistributedLoadAtNode")
SUBROUTINE foamUpdateLoad(ig,F,M) bind(c, name="beamDynSetDistributedLoad")
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(in) :: ig !inode0
    real(R8Ki), dimension(3), intent(in) :: F, M
    integer :: idx

!    if(inode0 < 0 .or. inode0 >= BD_Parameter%node_total) then
!        write(*,*) 'Invalid node number',inode0,BD_Parameter%node_total
!        return
!    end if
!    idx = inode0 + 1

    if(ig < 0 .or. ig >= order_elem) then
        write(*,*) 'Invalid gauss pt',ig,order_elem
        return
    end if
    idx = ig + 1

    foamUpdatedNode(idx) = foamUpdatedNode(idx) + 1

!     F_foam_t00(1:3,idx) = F_foam_t0(1:3,idx)
!     M_foam_t00(1:3,idx) = M_foam_t0(1:3,idx)
 
    F_foam_t0(1:3,idx) = F_foam(1:3,idx)
    M_foam_t0(1:3,idx) = M_foam(1:3,idx)
 
    F_foam(1:3,idx) = F
    M_foam(1:3,idx) = M

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE get_nnodes(nnodes) bind(c, name="beamDynGetNnodes")
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(out) :: nnodes
    nnodes = BD_Parameter%node_total
    !write(*,*) 'get_nnodes : ', nnodes

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE get_disp_at_node(inode,pos,crv) bind(c, name="beamDynGetDispAtNode")
!
!**** NOTE ****
! 'pos' IS ACTUALLY THE DISPLACEMENT
! 'crv' IS THE ROTATED (FINAL) ORIENTATION = R*R0, R=transformation, R0=initial orientation
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(in) :: inode
    real(R8Ki), dimension(3), intent(out) :: pos, crv

    real(R8Ki), parameter :: four=4
    integer :: i

    if( inode > BD_Parameter%node_total ) then
        write(*,*) 'BEAMDYN WARNING: Attempted to retrieve displacement at station', &
            inode, '>', BD_Parameter%node_total
        pos(:) = 0.0
        crv(:) = 0.0
        return
    end if

    ! note: BD_Output(1) is at previous time step (i.e. time t)
    !       BD_Output(2) coincides with time t + dt and is updated through BeamDyn_UpdateStates
    pos(:) = BD_Output(2)%BldMotion%TranslationDisp(1:3,inode)
    CALL CrvExtractCrv(BD_Output(2)%BldMotion%Orientation(1:3,1:3,inode),crv)

    ! convert from rotation parameter to rotation angle
    do i=1,3
        ! P = 4 Tan(\phi/4.0); P = rotation param, phi = rotation angle
        crv(i) = four * atan( crv(i)/four )
    end do

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE get_node0_position(inode0,pos,crv) bind(c, name="beamDynGetNode0Position")
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(in) :: inode0
    real(R8Ki), dimension(3), intent(out) :: pos, crv
    integer :: j

    call get_disp_at_node(inode0+1,pos,crv)

    do j=1,3
        pos(j) = pos(j) + BD_Parameter%uuN0(inode0*BD_Parameter%dof_node+j,1) ! assume only 1 element for now
    end do

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE get_ngp(ngp) bind(c, name="beamDynGetNgp")
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(out) :: ngp
    ngp = BD_Parameter%ngp

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
!SUBROUTINE get_all_GLL_position(gp, gw) bind(c, name="beamDynGetGLLPts")
SUBROUTINE get_all_gauss_pts(gp, gw) bind(c, name="beamDynGetGaussPts")
    real(R8Ki), dimension(BD_Parameter%ngp), intent(out) :: gp, gw
    call BldGaussPointWeight(BD_Parameter%ngp, gp, gw)

END SUBROUTINE

SUBROUTINE get_all_GLL_position(gllp, gllw) bind(c, name="beamDynGetGLLPts")
    real(R8Ki), dimension(BD_Parameter%ngp+1), intent(out) :: gllp, gllw
    call BD_gen_gll_LSGL(BD_Parameter%ngp, gllp, gllw)
END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
! based on SUBROUTINE BldComputeJacobianLSGL(rr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,jacobian)
SUBROUTINE get_shape_functions(s,Ps) bind(c, name="beamDynGetShapeFunctions")
!
!********************************************************************************************!
!********************************************************************************************!
!   REAL(ReKi),    INTENT(IN   )::  rr            ! rrth Gauss point location 
!   REAL(ReKi),    INTENT(IN   )::  Nuu0(:)       ! Element nodal initial position
!   REAL(ReKi),    INTENT(IN   )::  gp(:)         ! Gauss point location
   real(ReKi),    intent(in   )::  s             ! probe location
!   REAL(ReKi),    INTENT(IN   )::  GLL_temp(:)   ! Gauss-Lobatto-Legendre point location (i.e. nodes)
!   INTEGER(IntKi),INTENT(IN   )::  node_elem     ! Number of node per element
!   INTEGER(IntKi),INTENT(IN   )::  dof_node      ! Number of DoF per node
!   INTEGER(IntKi),INTENT(IN   )::  ngp           ! Total number of Gauss point
!   INTEGER(IntKi)              ::  igp           ! ith Gauss point
!   REAL(ReKi),    INTENT(  OUT):: hhx(:)         ! Shape function
!   REAL(ReKi),    INTENT(  OUT):: hpx(:)         ! derivative of shape function
   REAL(ReKi), intent(out) :: Ps(order_elem+1) ! shape function
   REAL(ReKi)              :: dPhis(order_elem+1)
   REAL(ReKi)              :: jac            ! Jacobian of element

   ! Local variables
   REAL(ReKi)            :: Gup0(3)
   INTEGER(IntKi)        :: inode
   INTEGER(IntKi)        :: temp_id
   INTEGER(IntKi)        :: i,j,k,l
   REAL(ReKi)            :: dnum,den
   REAL(ReKi), PARAMETER :: eps = 1.0D-08

!   REAL(ReKi)            :: dPhis(order_elem+1,ngp), Ps(order_elem+1,ngp)

   REAL(ReKi), dimension(BD_Parameter%node_elem) :: GLL_temp, w_temp   ! Gauss-Lobatto-Legendre point location (i.e. nodes)

   call BD_gen_gll_LSGL(order_elem,GLL_temp,w_temp)

! Calculate the shape function at each Gauss pt
!   CALL diffmtc(order_elem,ngp,gp,GLL_temp,igp,hhx,hpx)
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!   Calculate Lagrangian interpolant tensor at ngp points where basis functions
!   are assumed to be associated with (order_elem+1) GLL points on [-1,1]
!
   do l = 1,order_elem+1
       den = 1.
!     do j = 1,ngp
!       dPhis(l,j) = 0.
!       if ((abs(gp(j)-1.).LE.eps).AND.(l.EQ.order_elem+1)) then
!         dPhis(l,j) = float((order_elem+1)*order_elem)/4.
!       elseif ((abs(gp(j)+1.).LE.eps).AND.(l.EQ.1)) then
!         dPhis(l,j) = -float((order_elem+1)*order_elem)/4.
!       elseif (abs(gp(j)-GLL_temp(l)).LE.eps) then
!         dPhis(l,j) = 0.
       dPhis(l) = 0.
       if ((abs(s-1.).LE.eps).AND.(l.EQ.order_elem+1)) then
         dPhis(l) = float((order_elem+1)*order_elem)/4.
       elseif ((abs(s-1.).LE.eps).AND.(l.EQ.1)) then
         dPhis(l) = -float((order_elem+1)*order_elem)/4.
       elseif (abs(s-GLL_temp(l)).LE.eps) then
         dPhis(l) = 0.
       else
         do i = 1,order_elem+1
           if (i.NE.l) then
             den = den*(GLL_temp(l)-GLL_temp(i))
           endif
           dnum = 1.
           do k = 1,order_elem+1
             if ((k.NE.l).AND.(k.NE.i).AND.(i.NE.l)) then
!               dnum = dnum*(gp(j)-GLL_temp(k))
               dnum = dnum*(s-GLL_temp(k))
             elseif (i.EQ.l) then
               dnum = 0.
             endif
           enddo
!           dPhis(l,j) = dPhis(l,j) + dnum
           dPhis(l) = dPhis(l) + dnum
         enddo
!         dPhis(l,j) = dPhis(l,j)/den
         dPhis(l) = dPhis(l)/den
       endif
!     enddo
   enddo

   do l = 1,order_elem+1
!     do j = 1,ngp
!       Ps(l,j) = 0.
       Ps(l) = 0.
       dnum = 1.
       den = 1.
!       if(abs(gp(j)-GLL_temp(l)).LE.eps) then
!         Ps(l,j) = 1.
       if(abs(s-GLL_temp(l)).LE.eps) then
         Ps(l) = 1.
       else
         do k = 1,order_elem+1
           if (k.NE.l) then
             den = den*(GLL_temp(l) - GLL_temp(k))
!             dnum = dnum*(gp(j) - GLL_temp(k))
             dnum = dnum*(s - GLL_temp(k))
           endif
         enddo
!         Ps(l,j) = dnum/den
         Ps(l) = dnum/den
       endif
!     enddo
   enddo
   
!   DO i=1,order_elem+1
!      hhx(i) = Ps(i,igp)
!      hpx(i) = dPhis(i,igp)
!   ENDDO

! end of modified diffmtc code
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!  Gup0 = 0.0D0
!  DO inode=1,BD_Parameter%node_elem
!      temp_id = (inode-1)*BD_Parameter%dof_node
!      DO i=1,3
!          !Gup0(i) = Gup0(i) + hpx(inode)*Nuu0(temp_id+i)
!          Gup0(i) = Gup0(i) + dPhis(inode)*BD_Parameter%uuN0(temp_id+i,1) ! assume 1 element for now
!      ENDDO
!  ENDDO
!
!  jac = SQRT(DOT_PRODUCT(Gup0,Gup0))
!  write(*,*) 'DEBUG h',s,Ps,jac
!  do inode=1,BD_Parameter%node_elem
!    Ps(inode) = jac*Ps(inode)
!  end do

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE calc_deriv( N, x, y, yp )
! used by interpolate_disp, deprecated as of 5/18/15
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(in) :: N
    real(8), dimension(N), intent(in) :: x, y
    real(8), dimension(N), intent(out) :: yp

    yp(1) = (y(2)-y(1)) / (x(2)-x(1))
    do i=2,N-1
        yp(i) = (y(i+1)-y(i-1)) / (x(i+1)-x(i-1))
    end do
    yp(N) = (y(N)-y(N-1)) / (x(N)-x(N-1))

END SUBROUTINE


!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE interpolate_disp(N,s,u,v,w,kx,ky,kz) bind(c, name="beamDynGetDisp")
! deprecated as of 5/18/15
!
!********************************************************************************************!
!********************************************************************************************!
    integer, intent(in)      :: N !number of output values
    real(R8Ki), dimension(N), intent(in) :: s !parameterized position
    real(R8Ki), dimension(N), intent(inout) :: u, v, w !displacements
    real(R8Ki), dimension(N), intent(inout) :: kx,ky,kz !curvatures

    ! variables passed to the spline library
    real(8), dimension(BD_Parameter%node_total) :: x, y, yp
    real(8) :: dummy
    real(8), dimension(4,BD_Parameter%node_total) :: c_u1, c_u2, c_u3
    real(8), dimension(4,BD_Parameter%node_total) :: c_k1, c_k2, c_k3

    !real(R8Ki), dimension(BD_Parameter%node_total) :: sn
    real(R8Ki), dimension(BD_Parameter%node_total,3) :: posn, crvn
    integer :: i, nnodes

    write(*,*) 'BEAMDYN: Warning, beamDynGetDisp is deprecated!'

    nnodes = BD_Parameter%node_total

    ! update nodal values
    do i=1,nnodes
        x(i) = real(i-1,8)/real(nnodes-1,8) ! assume equally spaced nodes
        call get_disp_at_node(i,posn(i,:),crvn(i,:)) ! update pos, crv
    end do

    ! interpolate displacement u1
    call calc_deriv( nnodes, x, posn(:,1), yp ) !update yp
    call spline_hermite_set( nnodes, x, posn(:,1), yp, c_u1) !sets up interpolant c
    do i=1,N
        call spline_hermite_val( nnodes, x, c_u1, s(i), u(i), dummy ) ! evaluate interpolants
    end do

    ! interpolate displacement u2
    call calc_deriv( nnodes, x, posn(:,2), yp ) !update yp
    call spline_hermite_set( nnodes, x, posn(:,2), yp, c_u2) !sets up interpolant c
    do i=1,N
        call spline_hermite_val( nnodes, x, c_u2, s(i), v(i), dummy ) ! evaluate interpolants
    end do

    ! interpolate displacement u3
    call calc_deriv( nnodes, x, posn(:,3), yp ) !update yp
    call spline_hermite_set( nnodes, x, posn(:,3), yp, c_u3) !sets up interpolant c
    do i=1,N
        call spline_hermite_val( nnodes, x, c_u3, s(i), w(i), dummy ) ! evaluate interpolants
    end do

    ! interpolate curvature k1
    call calc_deriv( nnodes, x, crvn(:,1), yp ) !update yp
    call spline_hermite_set( nnodes, x, crvn(:,1), yp, c_k1) !sets up interpolant c
    do i=1,N
        call spline_hermite_val( nnodes, x, c_k1, s(i), kx(i), dummy ) ! evaluate interpolants
    end do

    ! interpolate curvature k2
    call calc_deriv( nnodes, x, crvn(:,2), yp ) !update yp
    call spline_hermite_set( nnodes, x, crvn(:,2), yp, c_k2) !sets up interpolant c
    do i=1,N
        call spline_hermite_val( nnodes, x, c_k2, s(i), ky(i), dummy ) ! evaluate interpolants
    end do

    ! interpolate curvature k3
    call calc_deriv( nnodes, x, crvn(:,3), yp ) !update yp
    call spline_hermite_set( nnodes, x, crvn(:,3), yp, c_k3) !sets up interpolant c
    do i=1,N
        call spline_hermite_val( nnodes, x, c_k3, s(i), kz(i), dummy ) ! evaluate interpolants
    end do

END SUBROUTINE

   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

!********************************************************************************************!
!********************************************************************************************!
!
SUBROUTINE wrap_up() bind(c, name="beamDynEnd")
!
!********************************************************************************************!
!********************************************************************************************!
   write(*,*) 'BEAMDYN: Cleaning up'

   CALL BeamDyn_End( BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                    BD_ConstraintState, BD_OtherState, BD_Output(1), ErrStat, ErrMsg )

   DO i = 2, BD_interp_order+1
      CALL BD_DestroyInput(BD_Input(i), ErrStat, ErrMsg )
      CALL BD_DestroyOutput(BD_Output(i), ErrStat, ErrMsg )
   ENDDO

   DEALLOCATE(BD_Input)
   DEALLOCATE(BD_InputTimes)
   DEALLOCATE(BD_Output)
   DEALLOCATE(BD_OutputTimes)
   DEALLOCATE(BD_InitInput%gravity)

   deallocate(F_foam)
   deallocate(M_foam)
   deallocate(F_foam_t0)
   deallocate(M_foam_t0)
!   deallocate(F_foam_t00)
!   deallocate(M_foam_t00)
   deallocate(foamUpdatedNode)

   if( writeOutput ) then
       CLOSE (QiDisUnit)
       CLOSE (QiRootUnit)
       CLOSE (QiReacUnit)
   end if

!CLOSE (QiHUnit)

END SUBROUTINE wrap_up


END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!  PROBLEM SPECIFIC MODIFICATIONS BELOW HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    
!    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    
!    |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |    
!   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /   \ /
!    V     V     V     V     V     V     V     V     V     V     V     V     V     V     V     V


!  TODO: specifies root motion!
!   SUBROUTINE BD_InputSolve( t, u, ut, p, ErrStat, ErrMsg)
   SUBROUTINE BD_InputSolve( t, u, ut, p, F, M, ErrStat, ErrMsg)
 
   USE BeamDyn_SP
   USE BeamDyn_Types

   REAL(DbKi),                     INTENT(IN   ):: t            ! time
   TYPE(BD_InputType),             INTENT(INOUT):: u            ! root motion
   REAL(DbKi),                     INTENT(INOUT):: ut           ! (don't worry about it)
   TYPE(BD_ParameterType),         INTENT(IN   ):: p            ! fixed parameters for the simulation
   REAL(DbKi), dimension(3,p%node_total), intent(in) :: F, M
   INTEGER(IntKi),                 INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)          :: i                ! do-loop counter
!   REAL(ReKi)              :: temp_vec(3)
!   REAL(ReKi)              :: temp_vec2(3)
!   REAL(ReKi)              :: temp_rr(3)
!   REAL(ReKi)              :: temp_pp(3)
!   REAL(ReKi)              :: temp_qq(3)
!   REAL(ReKi)              :: temp_R(3,3)
   REAL(ReKi)              :: pai
   REAL(ReKi)              :: omega

   write(*,*) 'BEAMDYN: BD_InputSolve called to specify problem-specific root motion'

   ErrStat = ErrID_None
   ErrMsg  = ''

   pai = ACOS(-1.0D0)
   omega = pai/5.0D0
!   temp_rr(:)     = 0.0D0
!   temp_pp(:)     = 0.0D0
!   temp_qq(:)     = 0.0D0
!   temp_R(:,:)    = 0.0D0
   ! gather point forces and line forces

!------------------
!  Rotating beam
!------------------
   ! Point mesh: RootMotion 
   ! Calculate root displacements and rotations
   u%RootMotion%TranslationDisp(:,:)  = 0.0D0
   u%RootMotion%Orientation(:,:,:) = 0.0D0
   DO i=1,3
       u%RootMotion%Orientation(i,i,1) = 1.0D0
   ENDDO
   ! END Calculate root displacements and rotations

   ! Calculate root translational and angular velocities
   u%RootMotion%TranslationVel(:,:) = 0.0D0
   u%RootMotion%RotationVel(:,:) = 0.0D0
   ! END Calculate root translational and angular velocities


   ! Calculate root translational and angular accelerations
   u%RootMotion%TranslationAcc(:,:) = 0.0D0
   u%RootMotion%RotationAcc(:,:) = 0.0D0
   ! END Calculate root translational and angular accelerations
!------------------
! End rotating beam
!------------------

   u%PointLoad%Force(:,:)  = 0.0D0
   u%PointLoad%Moment(:,:) = 0.0D0
   
! FOR DRIVER CODE
!   IF(t .LT. 0.05) THEN
!       u%PointLoad%Force(2,p%node_total) = 8.0D+04*t
!   ELSEIF(t .GE. 0.05 .AND. t .LT. 0.1) THEN
!       u%PointLoad%Force(2,p%node_total) = -8.0D+04*t+8D+03
!   ELSE
!       u%PointLoad%Force(2,p%node_total) = 0.0D0
!   ENDIF

   ! LINE2 mesh: DistrLoad
   u%DistrLoad%Force(:,:)  = 0.0D0
   u%DistrLoad%Moment(:,:) = 0.0D0

! FOR COUPLING WITH OPENFOAM
! NOTE:
! distributed loads have size (1:3,1:order_elem+2), by the FAST convention
!   entries (1:3,2:order_elem+1) correspond to order_elem number of gauss pts
!   do j=1,p%node_total
   do j=1,p%node_elem-1
       do i=1,3
            u%DistrLoad%Force(i,j+1) = F(i,j)
           u%DistrLoad%Moment(i,j+1) = M(i,j)
       end do
   end do

   ut = t


   END SUBROUTINE BD_InputSolve


   ! TODO: problem specific
   ! called once for ICs
   SUBROUTINE BD_InitialCondition(u,p,x,ErrStat,ErrMsg)

   USE BeamDyn_SP
   USE BeamDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InputType),                     INTENT(IN   ):: u            ! root motion
   TYPE(BD_ParameterType),                 INTENT(IN   ):: p            ! fixed parameters for the simulation
   TYPE(BD_ContinuousStateType),           INTENT(INOUT):: x            ! contains q: displacement, dqdt: velocity
   INTEGER(IntKi),                         INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                           INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                       :: i
   INTEGER(IntKi)                                       :: j
   INTEGER(IntKi)                                       :: k
   INTEGER(IntKi)                                       :: temp_id
   INTEGER(IntKi)                                       :: temp_id2
   REAL(ReKi)                                           :: temp66(6,6)

   write(*,*) 'BEAMDYN: BD_InitialCondition called to specify problem-specific initial condition'

   write(*,*) 'BEAMDYN: --TODO-- restarts not handled yet!'

   ErrStat = ErrID_None
   ErrMsg  = ''
   ! Initial displacements and rotations
   x%q(:) = 0.0D0
   ! Initial velocities and angular velocities
!   DO i=1,p%elem_total
!       DO j=1,p%node_elem
!           temp_id = (i-1)*p%dof_node*p%node_elem+(j-1)*p%dof_node
!           temp_id2= (j-1)*p%dof_node
!           x%dqdt(temp_id+1:temp_id+3) = &
!           MATMUL(Tilde(u%RootMotion%RotationVel(1:3,1)),p%GlbPos(1:3)+MATMUL(p%GlbRot,p%uuN0(temp_id2+1:temp_id2+3,i)))
!           x%dqdt(temp_id+4:temp_id+6) = u%RootMotion%RotationVel(1:3,1)
!           CALL MotionTensor(p%GlbRot,p%GlbPos,temp66,1)
!           x%dqdt(temp_id+1:temp_id+6) = MATMUL(temp66,x%dqdt(temp_id+1:temp_id+6))
!       ENDDO
!   ENDDO
    x%dqdt(:) = 0.0D0
   END SUBROUTINE BD_InitialCondition
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
