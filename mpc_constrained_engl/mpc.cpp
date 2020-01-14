/**************************************************************************************************
 * Class for MPC without constraint
 *
 *  The plant to be controlled is a Linear Time-Invariant System:
 *          x(k+1)  = A*x(k) + B*u(k)   ; x = Nx1, u = Mx1
 *          z(k)    = C*x(k)            ; z = Zx1
 *
 *
 ** Calculate prediction of x(k+1..k+Hp) constants ************************************************
 *
 *      Prediction of state variable of the system:
 *        x(k+1..k+Hp) = PSI*x(k) + OMEGA*u(k-1) + THETA*dU(k..k+Hu-1)                  ...{MPC_1}
 *
 *        Constants:
 *          PSI   = [A A^2 ... A^Hp]'                                         : (Hp*N)xN
 *          OMEGA = [B B+A*B ... Sigma(i=0->Hp-1)A^i*B]'                      : (Hp*N)xM
 *          THETA = [         B               0  ....           0            ]
 *                  [       B+A*B             B   .             0            ]
 *                  [         .               .    .            B            ]: (Hp*N)x(Hu*M)
 *                  [         .               .     .           .            ]
 *                  [Sigma(i=0->Hp-1)(A^i*B)  .  ....  Sigma(i=0->Hp-Hu)A^i*B]
 *
 *
 ** Calculate prediction of z(k+1..k+Hp) constants ************************************************
 *
 *      Prediction of output of the system:
 *        z(k+1..k+Hp) = (Cz*PSI)*x(k) + (Cz*OMEGA)*u(k-1) + (Cz*THETA)*dU(k..k+Hu-1)   ...{MPC_2}
 *
 *        Constants:
 *          Cz      : [C 0 0 .. 0]    ; 0 = zero(ZxN)       : (Hp*Z)x(Hp*N)
 *                    [0 C 0 .. 0]
 *                    [0 0 C .. 0]
 *                    [. . . .. 0]
 *                    [0 0 0 .. C]
 *
 *          CPSI   = Cz*PSI                                 : (Hp*Z)xN
 *          COMEGA = Cz*OMEGA                               : (Hp*Z)xM
 *          CTHETA = Cz*THETA                               : (Hp*Z)x(Hu*M)
 *
 *
 ** Calculate offline optimization variable H *****************************************************
 * 
 *      H = CTHETA'*Q*CTHETA + R                                                        ...{MPC_3}
 * 
 ** MPC update algorithm **************************************************************************
 *
 *      Formulation of plant error prediction
 *          E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                    ...{MPC_4}
 * 
 *      Calculate online optimization variable G:
 *          G = 2*CTHETA'*Q*E(k)                                                        ...{MPC_5}
 * 
 *      Create the inequality constraints                                               ...{MPC_6}
 *          dUmin <= dU(k) <= dUmax                                                     ...{MPC_6a}
 *          Umin <= U(kj <= Umax                                                        ...{MPC_6b}
 *          Zmin <= Z(kj <= Zmax                                                        ...{MPC_6c}
 * 
 *      MPC solution:
 *          (a) For unconstrained MPC:
 *              d [dU(k)'*H*dU(k) - G'*dU(k)]
 *              ----------------------------- = 0   -->   2*H*dU(k)-G = 0
 *                      d[dU(k)]
 *
 *              --> dU(k)_optimal = 1/2 * H^-1 * G
 *              --> https://github.com/pronenewbits/Arduino_Unconstrained_MPC_Library
 *
 *          (b) For constrained MPC (quadrating programming):
 *                  min   dU(k)'*H*dU(k) - G'*dU(k)     ; subject to inequality constraints
 *                 dU(k)
 * 
 *              Reconditioning H & G matrices into a standard QP form:
 *                x_opt = arg Min  1/2*x'*Q*x + c'*x    ; subject to ineqLHS*x <= ineqRHS
 *                           x
 *              
 *              Then solve by Active Set:
 *                  dU_opt(k) = ActiveSet(2H, -G, ineqLHS, ineqRHS)                     ...{MPC_7}
 *
 *      Integrate the du(k) to get u(k):
 *              u(k) = u(k-1) + du(k)                                                   ...{MPC_8}
 *
 * 
 * 
 * See https://github.com/pronenewbits for more!
 *************************************************************************************************/
#include "mpc.h"


MPC::MPC(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR)
{
    vReInit(A, B, C, _bobotQ, _bobotR);
}

void MPC::vReInit(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR)
{
    this->A = A;
    this->B = B;
    this->C = C;
    Q.vIsiDiagonal(_bobotQ);
    R.vIsiDiagonal(_bobotR);

#if (0)
    /*  Calculate prediction of x(k+1..k+Hp) constants
     *
     *      Prediction of state variable of the system:
     *        x(k+1..k+Hp) = PSI*x(k) + OMEGA*u(k-1) + THETA*dU(k..k+Hu-1)                  ...{MPC_1}
     *
     *        Constants:
     *          PSI   = [A A^2 ... A^Hp]'                                         : (Hp*N)xN
     *          OMEGA = [B B+A*B ... Sigma(i=0->Hp-1)A^i*B]'                      : (Hp*N)xM
     *          THETA = [         B               0  ....           0            ]
     *                  [       B+A*B             B   .             0            ]
     *                  [         .               .    .            B            ]: (Hp*N)x(Hu*M)
     *                  [         .               .     .           .            ]
     *                  [Sigma(i=0->Hp-1)(A^i*B)  .  ....  Sigma(i=0->Hp-Hu)A^i*B]
     *
     */
    Matrix _PSI     ((MPC_HP_LEN*SS_X_LEN), SS_X_LEN);
    Matrix _OMEGA   ((MPC_HP_LEN*SS_X_LEN), SS_U_LEN);
    Matrix _THETA   ((MPC_HP_LEN*SS_X_LEN), (MPC_HU_LEN*SS_U_LEN));

    Matrix _Apow(SS_X_LEN, SS_X_LEN);
    /* PSI      : [  A  ]
     *            [ A^2 ]
     *            [  .  ]                                                   : (Hp*N) x N
     *            [  .  ]
     *            [A^Hp ]
     */
    _Apow = A;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        _PSI = _PSI.InsertSubMatrix(_Apow, _i*SS_X_LEN, 0);
        _Apow = _Apow * A;
    }
    
    
    /* OMEGA    : [          B          ]
     *            [        B+A*B        ]
     *            [          .          ]                                   : (Hp*N) x M
     *            [          .          ]
     *            [Sigma(i=0->Hp-1)A^i*B]
     */
    Matrix _tempSigma(SS_X_LEN, SS_U_LEN);
    _Apow.vSetIdentitas();
    _tempSigma = B;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        _OMEGA = _OMEGA.InsertSubMatrix(_tempSigma, _i*SS_X_LEN, 0);
        _Apow = _Apow * A;
        _tempSigma = _tempSigma + (_Apow*B);
    }
    
    
    /* THETA    : [         B               0  ....           0            ]
     *            [       B+A*B             B   .             0            ]
     *            [         .               .    .            B            ]: (Hp*N)x(Hu*M)
     *            [         .               .     .           .            ]
     *            [Sigma(i=0->Hp-1)A^i*B    .  ....  Sigma(i=0->Hp-Hu)A^i*B]
     *
     *          : [OMEGA   [0 OMEGA(0:(len(OMEGA)-len(B)),:)]'  ....  [0..0 OMEGA(0:(len(OMEGA)-((Hp-Hu)*len(B))),:)]']
     */
    for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
        _THETA = _THETA.InsertSubMatrix(_OMEGA, _i*SS_X_LEN, _i*SS_U_LEN, (MPC_HP_LEN*SS_X_LEN)-(_i*SS_X_LEN), SS_U_LEN);
    }


    /* Calculate prediction of z(k+1..k+Hp) constants
     *
     *      Prediction of output of the system:
     *        z(k+1..k+Hp) = (Cz*PSI)*x(k) + (Cz*OMEGA)*u(k-1) + (Cz*THETA)*dU(k..k+Hu-1)   ...{MPC_2}
     *
     *        Constants:
     *          Cz      : [C 0 0 .. 0]    ; 0 = zero(ZxN)       : (Hp*Z)x(Hp*N)
     *                    [0 C 0 .. 0]
     *                    [0 0 C .. 0]
     *                    [. . . .. 0]
     *                    [0 0 0 .. C]
     *
     *          CPSI   = Cz*PSI                                 : (Hp*Z)xN
     *          COMEGA = Cz*OMEGA                               : (Hp*Z)xM
     *          CTHETA = Cz*THETA                               : (Hp*Z)x(Hu*M)
     */
    Matrix Cz(MPC_HP_LEN*SS_Z_LEN, MPC_HP_LEN*SS_X_LEN);
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        Cz = Cz.InsertSubMatrix(C, _i*SS_Z_LEN, _i*SS_X_LEN);
    }
    CPSI    = Cz * _PSI;
    COMEGA  = Cz * _OMEGA;
    CTHETA  = Cz * _THETA;
#else
    Matrix _Apow(SS_X_LEN, SS_X_LEN);
    /* CPSI     : [ C *   A  ]
     *            [ C *  A^2 ]
     *            [     .    ]                                                   : (Hp*N) x N
     *            [     .    ]
     *            [ C * A^Hp ]
     */
    _Apow = A;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        CPSI = CPSI.InsertSubMatrix((C*_Apow), _i*SS_Z_LEN, 0);
        _Apow = _Apow * A;
    }
    
    
    /* OMEGA    : [           C * B          ]
     *            [         C * B+A*B        ]
     *            [             .            ]                                   : (Hp*N) x M
     *            [             .            ]
     *            [ C * Sigma(i=0->Hp-1)A^i*B]
     */
    Matrix _tempSigma(SS_X_LEN, SS_U_LEN);
    _Apow.vSetIdentitas();
    _tempSigma = B;
    for (int32_t _i = 0; _i < MPC_HP_LEN; _i++) {
        COMEGA = COMEGA.InsertSubMatrix((C*_tempSigma), _i*SS_Z_LEN, 0);
        _Apow = _Apow * A;
        _tempSigma = _tempSigma + (_Apow*B);
    }
    
    /* THETA    : [         B               0  ....           0            ]
     *            [       B+A*B             B   .             0            ]
     *            [         .               .    .            B            ]: (Hp*N)x(Hu*M)
     *            [         .               .     .           .            ]
     *            [Sigma(i=0->Hp-1)A^i*B    .  ....  Sigma(i=0->Hp-Hu)A^i*B]
     *
     *          : [OMEGA   [0 OMEGA(0:(len(OMEGA)-len(B)),:)]'  ....  [0..0 OMEGA(0:(len(OMEGA)-((Hp-Hu)*len(B))),:)]']
     */
    for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
        CTHETA = CTHETA.InsertSubMatrix(COMEGA, _i*SS_Z_LEN, _i*SS_U_LEN, (MPC_HP_LEN*SS_Z_LEN)-(_i*SS_Z_LEN), SS_U_LEN);
    }
#endif
    
    
    /* Calculate offline optimization variable H:
     *  H = CTHETA'*Q*CTHETA + R                                                            ...{MPC_3}
     * 
     * Note: we need to reconditioning into standard QP form (Q_qp = 2H_mpc, c_qp = -G_mpc)
     */
    H = 2*((CTHETA.Transpose()) * Q * CTHETA) + R;
}


bool MPC::bUpdate(Matrix &SP, Matrix &x, Matrix &u)
{
    Matrix Err((MPC_HP_LEN*SS_Z_LEN), 1);
    Matrix G((MPC_HU_LEN*SS_U_LEN), 1);
    
    
    /* Formulation of plant error prediction:
     *  E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                            ...{MPC_4} 
     */
    Err = SP - CPSI*x - COMEGA*u;
    
    
    /* Calculate online optimization variable G:
     *  G = 2*CTHETA'*Q*E(k)                                                                ...{MPC_5}
     * 
     * Note: we need to reconditioning into standard QP form (Q_qp = 2H_mpc, c_qp = -G_mpc)
     */
    G = -(CTHETA.Transpose()) * Q * Err * 2.0;
    
    
    /*  Create the inequality constraints                                                   ...{MPC_6} */
    Matrix Constraint_LHS   {(INEQ_LEN_DU + INEQ_LEN_U + INEQ_LEN_Z), (MPC_HU_LEN*SS_U_LEN)};
    Matrix Constraint_RHS   {(INEQ_LEN_DU + INEQ_LEN_U + INEQ_LEN_Z), 1};
    #if (defined(MPC_USE_CONSTRAINT_DU) || defined(MPC_USE_CONSTRAINT_U) || defined(MPC_USE_CONSTRAINT_Z))
        int32_t _i32offsetCnstr = 0;
    #endif
    #if defined(MPC_USE_CONSTRAINT_DU)
        /* for dUmin <= dU(kj <= dUmax                                                      ...{MPC_6a} */
        {
            Matrix CnstDu_LHS   {(2*(MPC_HU_LEN*SS_U_LEN)), (MPC_HU_LEN*SS_U_LEN)};
            Matrix CnstDu_RHS   {(2*(MPC_HU_LEN*SS_U_LEN)), 1};
            Matrix CnstDu_Identity((MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN));
            CnstDu_Identity.vSetIdentitas();
            CnstDu_LHS = CnstDu_LHS.InsertSubMatrix(-CnstDu_Identity, 0, 0);
            CnstDu_LHS = CnstDu_LHS.InsertSubMatrix(CnstDu_Identity, (MPC_HU_LEN*SS_U_LEN), 0);
            
            Matrix CnstDu_MaxMin((MPC_HU_LEN*SS_U_LEN), 1);
            CnstDu_MaxMin.vIsiHomogen(MPC_MIN_DU);
            CnstDu_RHS = CnstDu_RHS.InsertSubMatrix(-CnstDu_MaxMin, 0, 0);
            CnstDu_MaxMin.vIsiHomogen(MPC_MAX_DU);
            CnstDu_RHS = CnstDu_RHS.InsertSubMatrix(CnstDu_MaxMin, (MPC_HU_LEN*SS_U_LEN), 0);
        
            /* Append to Constraint_LHS & Constraint_RHS */
            Constraint_LHS = Constraint_LHS.InsertSubMatrix(CnstDu_LHS, _i32offsetCnstr, 0);
            Constraint_RHS = Constraint_RHS.InsertSubMatrix(CnstDu_RHS, _i32offsetCnstr, 0);
            _i32offsetCnstr += CnstDu_LHS.i32getBaris();
        }
    #endif
    #if defined(MPC_USE_CONSTRAINT_U)
        /* for Umin <= U(kj <= Umax                                                         ...{MPC_6b} */
        {
            Matrix CnstU_LHS   {(2*(MPC_HU_LEN*SS_U_LEN)), (MPC_HU_LEN*SS_U_LEN)};
            Matrix CnstU_RHS   {(2*(MPC_HU_LEN*SS_U_LEN)), 1};
            Matrix U_Identity(SS_U_LEN, SS_U_LEN);
            U_Identity.vSetIdentitas();
            for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
                for (int32_t _j = 0; _j < MPC_HU_LEN; _j++) {
                    if (_j <= _i) {
                        /* for Umin <= U(k) */
                        CnstU_LHS = CnstU_LHS.InsertSubMatrix(-U_Identity, (_i*SS_U_LEN), (_j*SS_U_LEN));
                        /* for U(k) <= Umax */
                        CnstU_LHS = CnstU_LHS.InsertSubMatrix(U_Identity,  (MPC_HU_LEN*SS_U_LEN)+(_i*SS_U_LEN), (_j*SS_U_LEN));
                    }
                }
            }
            Matrix CnstU_MaxMin((MPC_HU_LEN*SS_U_LEN), 1);
            for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
                for (int32_t _j = 0; _j < SS_U_LEN; _j++) {
                    CnstU_MaxMin[(_i*SS_U_LEN)+_j][0] = -MPC_MIN_U + u[_j][0];
                }
            }
            CnstU_RHS = CnstU_RHS.InsertSubMatrix(CnstU_MaxMin, 0, 0);
            for (int32_t _i = 0; _i < MPC_HU_LEN; _i++) {
                for (int32_t _j = 0; _j < SS_U_LEN; _j++) {
                    CnstU_MaxMin[(_i*SS_U_LEN)+_j][0] = MPC_MAX_U - u[_j][0];
                }
            }
            CnstU_RHS = CnstU_RHS.InsertSubMatrix(CnstU_MaxMin, (MPC_HU_LEN*SS_U_LEN), 0);
        
            /* Append to Constraint_LHS & Constraint_RHS */
            Constraint_LHS = Constraint_LHS.InsertSubMatrix(CnstU_LHS, _i32offsetCnstr, 0);
            Constraint_RHS = Constraint_RHS.InsertSubMatrix(CnstU_RHS, _i32offsetCnstr, 0);
            _i32offsetCnstr += CnstU_LHS.i32getBaris();
        }
    #endif
    #if defined(MPC_USE_CONSTRAINT_Z)
        /* for Zmin <= Z(kj <= Zmax                                                         ...{MPC_6c} */
        {
            Matrix CnstZ_LHS   {(2*(MPC_HP_LEN*SS_Z_LEN)), (MPC_HU_LEN*SS_U_LEN)};
            CnstZ_LHS = CnstZ_LHS.InsertSubMatrix(-CTHETA, 0, 0);
            CnstZ_LHS = CnstZ_LHS.InsertSubMatrix(CTHETA, (MPC_HP_LEN*SS_Z_LEN), 0);
            
            Matrix CnstZ_RHS   {(2*(MPC_HP_LEN*SS_Z_LEN)), 1};
            
            Matrix CnstZ_Min((MPC_HP_LEN*SS_Z_LEN), 1);
            CnstZ_Min = CPSI*x + COMEGA*u - MPC_MIN_Z;
            CnstZ_RHS = CnstZ_RHS.InsertSubMatrix(CnstZ_Min, 0, 0);
            
            Matrix CnstZ_Max((MPC_HP_LEN*SS_Z_LEN), 1);
            CnstZ_Max = -CPSI*x + -COMEGA*u + MPC_MAX_Z;
            CnstZ_RHS = CnstZ_RHS.InsertSubMatrix(CnstZ_Max, (MPC_HP_LEN*SS_Z_LEN), 0);
        
            /* Append to Constraint_LHS & Constraint_RHS */
            Constraint_LHS = Constraint_LHS.InsertSubMatrix(CnstZ_LHS, _i32offsetCnstr, 0);
            Constraint_RHS = Constraint_RHS.InsertSubMatrix(CnstZ_RHS, _i32offsetCnstr, 0);
            _i32offsetCnstr += CnstZ_LHS.i32getBaris();
        }
    #endif
    
    
    #if (0)
        /* Experiment: disable constraints, the MPC will behave like unconstrained MPC */
        #warning("Contstraints bypassed (no constraints)");
        Constraint_LHS = Matrix(0,0);
        Constraint_RHS = Matrix(0,0);
    #endif
    
    
    /*  Formulation of the optimal control problem:
     *
     *      For constrained MPC (quadrating programming):
     *          min   dU(k)'*H*dU(k) - G'*dU(k)     ; subject to inequality constraints
     *             dU(k)
     * 
     *      Reconditioning H & G matrices into a standard QP form:
     *          x_opt = arg Min  1/2*x'*Q*x + c'*x  ; subject to ineqLHS*x <= ineqRHS
     *                     x
     * 
     *      Then solve by Active Set:
     *          dU_opt(k) = ActiveSet(2H, -G, ineqLHS, ineqRHS)                         ...{MPC_7}
     * 
     * 
     * 1/2*Q = H    --> Q = 2*H     (the reconditioning has been done when we create H & G)
     * c' = -G'     --> c = -G
     */
    if (!bActiveSet(DU, H, G, Constraint_LHS, Constraint_RHS, cntIterActiveSet)) {
        /* return false; */
        DU.vIsiNol();
        return false;
    }
    
    /* Integrate the du(k) to get u(k):
     *  u(k) = u(k-1) + du(k)                                                           ...{MPC_8} 
     */
    Matrix DU_Out(SS_U_LEN, 1);
    for (int32_t _i = 0; _i < SS_U_LEN; _i++) {
        DU_Out[_i][0] = DU[_i][0];
    }
    u = u + DU_Out;
    
    return true;
}


/* Active Set solver for Quadratic Programming problem in the form:
 * 
 *      x_opt = arg Min  1/2*x'*Q*x + c'*x    ; subject to ineqLHS*x <= ineqRHS
 *                 x
 * 
 * 
 * The Active Set: search x, by solving this minimization problem:
 * 
 *  min. 1/2*dx'*Q*dx + (Q*x+c)'*dx   , subject to ineqLHS_dx*x = 0
 *   dx
 * 
 *  Integrate x(iter+1) = x(iter) + dx(iter)
 * 
 *  Until KKT conditions is satisfied:
 *      1. Q*dx + ineqLHS_dx*(dLambda) = -(Q*x+c)   (dLambda = Lagrange multiplier of above minimization solution)
 *      2. -ineqLHS_dx*dx = 0
 * 
 */
bool MPC::bActiveSet(Matrix &x, Matrix &Q, Matrix &c, Matrix &ineqLHS, Matrix &ineqRHS, int16_t &_i16iterActiveSet)
{
    bool _flagConstActive[ineqRHS.i32getBaris()];    /* Contains information about which inequality is active (if true, then that inequality is active). */
    bool _dlambdaPos, _dxNol;
    Matrix dx(x.i32getBaris(), 1);
    
    /* TODO: Make sure x initial value is inside feasible region (e.g. using Linear Programming).
     *       For now we set it to zero --> no mathematical guarantee!
     */
    x.vIsiNol();
    
    /* In the beginning, every constraints is non-active */
    for (int32_t _i = 0; _i < ineqRHS.i32getBaris(); _i++) {
        _flagConstActive[_i] = false;
    }
    
    _i16iterActiveSet = 0;  /* Reset counter iteration Active Set */
    
    do {
        uint16_t _u16cntConstActive = 0;
        for (int32_t _i = 0; _i < ineqLHS.i32getBaris(); _i++) {
            if (_flagConstActive[_i] == true) {
                _u16cntConstActive++;
            }
        }
        
        /* Construct Active Set solution matrix (the KKT matrix):
         * 
         *  _KKT_LHS(k) * [  dx   ] = _KKT_RHS(k)
         *                [dLambda]
         * 
         *  _KKT_LHS(k) = [        Q           _ineqActiveLHS(k)']
         *                [_ineqActiveLHS(k)           0         ]
         * 
         *  _KKT_RHS(k) = [-(Q*x+c)]
         *                [    0   ]
         */
        int32_t _indexConstActive[_u16cntConstActive];                      /* Contains information about index of the active inequality constraints in ineqLHS */
        Matrix _ineqActiveLHS(_u16cntConstActive, ineqLHS.i32getKolom());   /* Aggregation of the row vectors of ineqLHS that is 'active' */
        
        for (int32_t _i = 0, _iterConstActive = 0; _i < ineqLHS.i32getBaris(); _i++) {
            if (_flagConstActive[_i] == true) {
                _ineqActiveLHS = _ineqActiveLHS.InsertSubMatrix(ineqLHS, _iterConstActive, 0, _i, 0, 1, ineqLHS.i32getKolom());
                _indexConstActive[_iterConstActive] = _i;

                _iterConstActive++;
                ASSERT((_iterConstActive <= _u16cntConstActive), "Bug on the active set: Create _ineqActiveLHS");
            }
        }
        Matrix _KKT_LHS((Q.i32getBaris()+_u16cntConstActive), (Q.i32getKolom()+_u16cntConstActive));
        Matrix _KKT_RHS((Q.i32getBaris()+_u16cntConstActive), 1);
        
        _KKT_LHS = _KKT_LHS.InsertSubMatrix(Q, 0, 0);
        _KKT_LHS = _KKT_LHS.InsertSubMatrix(_ineqActiveLHS, Q.i32getBaris(), 0);
        _KKT_LHS = _KKT_LHS.InsertSubMatrix(_ineqActiveLHS.Transpose(), 0, Q.i32getKolom());
        _KKT_RHS = _KKT_RHS.InsertSubMatrix(-((Q*x)+c), 0, 0);
        
        
        /* [   dx  ] = [        Q           _ineqActiveLHS(k)']^-1 * [-(Q*x+c)]
         * [dLambda]   [_ineqActiveLHS(k)           0         ]      [    0   ]
         * 
         */
        Matrix _KKTvector(_KKT_LHS.i32getBaris(), 1);
        _KKTvector = _KKT_LHS.Invers()*_KKT_RHS;
        if (!_KKTvector.bCekMatrixValid()) {
            return false;
        }
        dx = dx.InsertSubMatrix(_KKTvector, 0, 0, x.i32getBaris(), 1);
        
        Matrix dLambda(_u16cntConstActive, 1);
        dLambda = dLambda.InsertSubMatrix(_KKTvector, 0, 0, x.i32getBaris(), 0, _u16cntConstActive, 1);
        
        /* Karush–Kuhn–Tucker conditions:
         *  1. dx = 0
         *  2. dLambda >= 0
         */
        _dxNol = true;
        for (int32_t _i = 0; _i < Q.i32getBaris(); _i++) {
            if (fabs(dx[_i][0]) > float_prec(float_prec_ZERO_ECO)) {
                _dxNol = false;
            }
        }
        
        float_prec _lowestLambda = 0.0;
        int32_t _idxIneqLowestLambda = -1;
        _dlambdaPos = true;
        for (int32_t _i = 0; _i < _u16cntConstActive; _i++) {
            if (dLambda[_i][0] < -float_prec(float_prec_ZERO_ECO)) {
                _dlambdaPos = false;
                
                /* Some constraints become inactive, search constraint with lowest Lagrange multiplier so we can remove it */
                if (dLambda[_i][0] < _lowestLambda) {
                    _lowestLambda = dLambda[_i][0];
                    _idxIneqLowestLambda = _indexConstActive[_i];
                }
            }
        }    
        if (!_dlambdaPos) {
            /* remove the constraint that have lowest Lagrange multiplier */
            ASSERT((_idxIneqLowestLambda != -1), "Bug on the active set: remove ineq active");
            _flagConstActive[_idxIneqLowestLambda] = false;
        }
        
        if (_dxNol && _dlambdaPos) {
            break;
        }
        
        /* Search for violated constraints. If any, then get the lowest scalar value that make the searching direction feasible again */
        float_prec _alphaMin = 1.0;
        int32_t _idxAlphaMin = -1;
        for (int32_t _i = 0; _i < ineqLHS.i32getBaris(); _i++) {
            if (_flagConstActive[_i] == false) {
                /* For all constraints that is not active, we check if that constraint is violated */
                Matrix _ineqNonActiveLHS(1, ineqLHS.i32getKolom());
                _ineqNonActiveLHS = _ineqNonActiveLHS.InsertSubMatrix(ineqLHS, 0, 0, _i, 0, 1, ineqLHS.i32getKolom());
                
                if (float_prec((_ineqNonActiveLHS * dx)[0][0]) > float_prec_ZERO_ECO) {
                    /* We have a violated constraint. Make that constraint a candidate of active constraint */
                    float_prec _alphaCandidate = - ((_ineqNonActiveLHS * x)[0][0] - ineqRHS[_i][0]) / (_ineqNonActiveLHS * dx)[0][0];
                    
                    if (_alphaCandidate < _alphaMin) {
                        _alphaMin = _alphaCandidate;
                        _idxAlphaMin = _i;
                    }
                }
            }
        }
        x = x + (dx * _alphaMin);
        
        if (_idxAlphaMin != -1) {
            _flagConstActive[_idxAlphaMin] = true;
        }
        
        _i16iterActiveSet++;
        if (_i16iterActiveSet > MPC_MAXIMUM_ACTIVE_SET_ITERATION) {
            /* Maximum number of iteration reached, terminate and use last solution */
            break;
        }
    } while ((!_dxNol) || (!_dlambdaPos));
    
    return true;
}


