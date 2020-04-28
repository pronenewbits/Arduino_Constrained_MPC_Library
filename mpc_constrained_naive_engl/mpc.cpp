/**************************************************************************************************
 * Class for MPC with constraint
 *
 *  The plant to be controlled is a Linear Time-Invariant System:
 *          x(k+1)  = A*x(k) + B*u(k)   ; x = Nx1, u = Mx1
 *          z(k)    = C*x(k)            ; z = Zx1
 *
 *
 ** Calculate prediction of z(k+1..k+Hp) constants ************************************************
 *
 *      Prediction of state variable of the system:
 *        z(k+1..k+Hp) = (CPSI)*x(k) + (COMEGA)*u(k-1) + (CTHETA)*dU(k..k+Hu-1)         ...{MPC_1}
 *
 *        Constants:
 *          CPSI   = [CA C(A^2) ... C(A^Hp)]'                                       : (Hp*N)xN
 *          COMEGA = [CB C(B+A*B) ... C*Sigma(i=0->Hp-1)A^i*B]'                     : (Hp*N)xM
 *          CTHETA = [         CB                0  ....           0              ]
 *                   [       C(B+A*B)           CB   .             0              ]
 *                   [           .               .    .           CB              ] : (Hp*N)x(Hu*M)
 *                   [           .               .     .           .              ]
 *                   [C*Sigma(i=0->Hp-1)(A^i*B)  .  ....  C*Sigma(i=0->Hp-Hu)A^i*B]
 *
 *
 ** Calculate offline optimization variable H *****************************************************
 * 
 *      H = CTHETA'*Q*CTHETA + R                                                        ...{MPC_2}
 * 
 ** MPC update algorithm **************************************************************************
 *
 *      Formulation of plant error prediction
 *          E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                    ...{MPC_3}
 * 
 *      Calculate online optimization variable G:
 *          G = 2*CTHETA'*Q*E(k)                                                        ...{MPC_4}
 * 
 *      Create the inequality constraints                                               ...{MPC_5}
 *          dUmin <= dU(k) <= dUmax                                                     ...{MPC_5a}
 *          Umin <= U(kj <= Umax                                                        ...{MPC_5b}
 *          Zmin <= Z(kj <= Zmax                                                        ...{MPC_5c}
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
 *                  dU_opt(k) = ActiveSet(2H, -G, ineqLHS, ineqRHS)                     ...{MPC_6}
 *
 *      Integrate the du(k) to get u(k):
 *              u(k) = u(k-1) + du(k)                                                   ...{MPC_7}
 *
 * 
 * 
 * See https://github.com/pronenewbits for more!
 *************************************************************************************************/
#include "mpc.h"


MPC::MPC(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR, void (*vCreateConstraintsLHS)(Matrix &, Matrix &), void (*vCreateConstraintsRHS)(Matrix &, Matrix &, Matrix &, Matrix &, Matrix &))
{
    vReInit(A, B, C, _bobotQ, _bobotR, vCreateConstraintsLHS, vCreateConstraintsRHS);
}

void MPC::vReInit(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR, void (*vCreateConstraintsLHS)(Matrix &, Matrix &), void (*vCreateConstraintsRHS)(Matrix &, Matrix &, Matrix &, Matrix &, Matrix &))
{
    this->A = A;
    this->B = B;
    this->C = C;
    this->vCreateConstraintsLHS = vCreateConstraintsLHS;
    this->vCreateConstraintsRHS = vCreateConstraintsRHS;
    CnstLHS = Matrix(0,0);
    CnstRHS = Matrix(0,0);
    Q.vSetDiag(_bobotQ);
    R.vSetDiag(_bobotR);
    
    /*  Calculate prediction of z(k+1..k+Hp) constants
     *
     *      Prediction of state variable of the system:
     *        z(k+1..k+Hp) = (CPSI)*x(k) + (COMEGA)*u(k-1) + (CTHETA)*dU(k..k+Hu-1)         ...{MPC_1}
     *
     *        Constants:
     *          CPSI   = [CA C(A^2) ... C(A^Hp)]'                                       : (Hp*N)xN
     *          COMEGA = [CB C(B+A*B) ... C*Sigma(i=0->Hp-1)A^i*B]'                     : (Hp*N)xM
     *          CTHETA = [         CB                0  ....           0              ]
     *                   [       C(B+A*B)           CB   .             0              ]
     *                   [           .               .    .           CB              ] : (Hp*N)x(Hu*M)
     *                   [           .               .     .           .              ]
     *                   [C*Sigma(i=0->Hp-1)(A^i*B)  .  ....  C*Sigma(i=0->Hp-Hu)A^i*B]
     *
     */
    Matrix _Apow(SS_X_LEN, SS_X_LEN);
    /* CPSI     : [ C *   A  ]
     *            [ C *  A^2 ]
     *            [     .    ]                                                   : (Hp*N) x N
     *            [     .    ]
     *            [ C * A^Hp ]
     */
    _Apow = A;
    for (int16_t _i = 0; _i < MPC_HP_LEN; _i++) {
        CPSI = CPSI.InsertSubMatrix((C*_Apow), _i*SS_Z_LEN, 0);
        _Apow = _Apow * A;
    }
    
    /* COMEGA   : [          C * (B)         ]
     *            [        C * (B+A*B)       ]
     *            [             .            ]                                   : (Hp*N) x M
     *            [             .            ]
     *            [ C * Sigma(i=0->Hp-1)A^i*B]
     */
    Matrix _tempSigma(SS_X_LEN, SS_U_LEN);
    _Apow.vSetIdentity();
    _tempSigma = B;
    for (int16_t _i = 0; _i < MPC_HP_LEN; _i++) {
        COMEGA = COMEGA.InsertSubMatrix((C*_tempSigma), _i*SS_Z_LEN, 0);
        _Apow = _Apow * A;
        _tempSigma = _tempSigma + (_Apow*B);
    }
    
    /* CTHETA   : [          C * (B)              0         ....              0             ]
     *            [       C * (B+A*B)           C * (B)      .                0             ]
     *            [            .                  .           .             C * (B)         ]: (Hp*N)x(Hu*M)
     *            [            .                  .            .              .             ]
     *            [C * Sigma(i=0->Hp-1)A^i*B      .         ....  C * Sigma(i=0->Hp-Hu)A^i*B]
     *
     *          : [COMEGA   [0 COMEGA(0:(len(COMEGA)-len(B)),:)]'  ....  [0..0 COMEGA(0:(len(COMEGA)-((Hp-Hu)*len(B))),:)]']
     */
    for (int16_t _i = 0; _i < MPC_HU_LEN; _i++) {
        CTHETA = CTHETA.InsertSubMatrix(COMEGA, _i*SS_Z_LEN, _i*SS_U_LEN, (MPC_HP_LEN*SS_Z_LEN)-(_i*SS_Z_LEN), SS_U_LEN);
    }
    
    
    /* Calculate offline optimization variable H:
     *  H = CTHETA'*Q*CTHETA + R                                                            ...{MPC_2}
     */
    H = ((CTHETA.Transpose()) * Q * CTHETA) + R;
}


bool MPC::bUpdate(Matrix &SP, Matrix &x, Matrix &u)
{
    Matrix Err((MPC_HP_LEN*SS_Z_LEN), 1, Matrix::NoInitMatZero);
    Matrix G((MPC_HU_LEN*SS_U_LEN), 1, Matrix::NoInitMatZero);
    
    /* Formulation of plant error prediction:
     *  E(k) = SP(k) - CPSI*x(k) - COMEGA*u(k-1)                                            ...{MPC_3} 
     */
    Err = SP - CPSI*x - COMEGA*u;
    
    
    /* Calculate online optimization variable G:
     *  G = 2*CTHETA'*Q*E(k)                                                                ...{MPC_4}
     */
    G = 2.0 * (CTHETA.Transpose()) * Q * Err;
    
    
    /*  Create the inequality constraints                                                   ...{MPC_5}
     *          dUmin <= dU(k) <= dUmax                                                     ...{MPC_5a}
     *          Umin <= U(kj <= Umax                                                        ...{MPC_5b}
     *          Zmin <= Z(kj <= Zmax                                                        ...{MPC_5c}
     */
    vCreateConstraintsLHS(CnstLHS, CTHETA);
    vCreateConstraintsRHS(CnstRHS, COMEGA, CPSI, u, x);
    
    
    #if (0)
        /* Experiment: disable constraints, the MPC will behave like unconstrained MPC */
        #warning("Contstraints bypassed (no constraints)");
        CnstLHS = Matrix(0,0);
        CnstRHS = Matrix(0,0);
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
     *          dU_opt(k) = ActiveSet(2H, -G, ineqLHS, ineqRHS)                             ...{MPC_6}
     * 
     * 
     * 1/2*Q = H    --> Q = 2*H
     * c' = -G'     --> c = -G
     */
    if (!bActiveSet(DU, 2*H, -G, CnstLHS, CnstRHS, cntIterActiveSet)) {
        /* return false; */
        DU.vSetToZero();
        return false;
    }
    
    
    /* Integrate the du(k) to get u(k):
     *  u(k) = u(k-1) + du(k)                                                               ...{MPC_7} 
     */
    Matrix DU_Out(SS_U_LEN, 1);
    for (int16_t _i = 0; _i < SS_U_LEN; _i++) {
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
bool MPC::bActiveSet(Matrix &x, const Matrix &Q, const Matrix &c, const Matrix &ineqLHS, const Matrix &ineqRHS, int16_t &_i16iterActiveSet)
{
    bool _flagConstActive[ineqRHS.i16getRow()];    /* Contains information about which inequality is active (if true, then that inequality is active). */
    bool _dlambdaPos, _dxNol;
    Matrix dx(x.i16getRow(), 1, Matrix::NoInitMatZero);
    
    
    /* TODO: Make sure x initial value is inside feasible region (e.g. using Linear Programming).
     *       For now we set it to zero --> no mathematical guarantee!
     */
    x.vSetToZero();
    /* In the beginning, every constraints is non-active */
    for (int16_t _i = 0; _i < ineqRHS.i16getRow(); _i++) {
        _flagConstActive[_i] = false;
    }
    
    
    _i16iterActiveSet = 0;  /* Reset counter iteration Active Set */
    do {
        /* Construct active set matrix Aw (_ineqActiveLHS) */
        uint16_t _u16cntConstActive = 0;
        for (int16_t _i = 0; _i < ineqLHS.i16getRow(); _i++) {
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
        int16_t _indexConstActive[_u16cntConstActive];                                          /* Contains information about index of the active inequality constraints in ineqLHS */
        Matrix _ineqActiveLHS(_u16cntConstActive, ineqLHS.i16getCol(), Matrix::NoInitMatZero);  /* Aggregation of the row vectors of ineqLHS that is 'active', i.e. the Aw matrix */
        
        for (int16_t _i = 0, _iterConstActive = 0; _i < ineqLHS.i16getRow(); _i++) {
            if (_flagConstActive[_i] == true) {
                _ineqActiveLHS = _ineqActiveLHS.InsertSubMatrix(ineqLHS, _iterConstActive, 0, _i, 0, 1, ineqLHS.i16getCol());
                _indexConstActive[_iterConstActive] = _i;

                _iterConstActive++;
                ASSERT((_iterConstActive <= _u16cntConstActive), "Bug on the active set: Create _ineqActiveLHS");
            }
        }
        Matrix _KKT_LHS((Q.i16getRow()+_u16cntConstActive), (Q.i16getCol()+_u16cntConstActive));
        Matrix _KKT_RHS((Q.i16getRow()+_u16cntConstActive), 1);
        
        _KKT_LHS = _KKT_LHS.InsertSubMatrix(Q, 0, 0);
        _KKT_LHS = _KKT_LHS.InsertSubMatrix(_ineqActiveLHS, Q.i16getRow(), 0);
        _KKT_LHS = _KKT_LHS.InsertSubMatrix(_ineqActiveLHS.Transpose(), 0, Q.i16getCol());
        _KKT_RHS = _KKT_RHS.InsertSubMatrix(-((Q*x)+c), 0, 0);
        
        /* [   dx  ] = [        Q           _ineqActiveLHS(k)']^-1 * [-(Q*x+c)]
         * [dLambda]   [_ineqActiveLHS(k)           0         ]      [    0   ]
         * 
         */
        Matrix _KKTvector(_KKT_LHS.Invers()*_KKT_RHS);
        if (!_KKTvector.bMatrixIsValid()) {
            return false;
        }
        dx = dx.InsertSubMatrix(_KKTvector, 0, 0, x.i16getRow(), 1);
        
        Matrix dLambda(_u16cntConstActive, 1, Matrix::NoInitMatZero);
        dLambda = dLambda.InsertSubMatrix(_KKTvector, 0, 0, x.i16getRow(), 0, _u16cntConstActive, 1);
        
        
        /* Check for Karush–Kuhn–Tucker conditions ------------------------------------------------------------------------------------------------ */
        /* Karush–Kuhn–Tucker conditions:
         *  1. dx == 0
         *  2. dLambda >= 0
         */
        /* Search for dx == 0 ----------------- */
        _dxNol = true;
        for (int16_t _i = 0; _i < Q.i16getRow(); _i++) {
            if (fabs(dx[_i][0]) > float_prec(float_prec_ZERO_ECO)) {
                _dxNol = false;
                break;
            }
        }
        /* Search for dLambda >= 0 ------------ */
        float_prec _lowestLambda = 0.0;
        int16_t _idxIneqLowestLambda = -1;
        _dlambdaPos = true;
        for (int16_t _i = 0; _i < _u16cntConstActive; _i++) {
            if (dLambda[_i][0] < _lowestLambda) {
                /* Some constraints become inactive, search constraint with lowest Lagrange multiplier so we can remove it */
                _lowestLambda = dLambda[_i][0];
                _idxIneqLowestLambda = _indexConstActive[_i];
                _dlambdaPos = false;
            }
        }
        /* Check it! -------------------------- */
        if (_dxNol && _dlambdaPos) {
            break;
        } 
        
        
        if (!_dlambdaPos) {
            /* There is a constraint that become inactive ----------------------------------------------------------------------------------------- */

            /* remove the constraint that have lowest Lagrange multiplier */
            ASSERT((_idxIneqLowestLambda != -1), "Bug on the active set: remove ineq active");
            _flagConstActive[_idxIneqLowestLambda] = false;
        }
        
        
        if (!_dxNol) {
            /* It means the x(k+1) will move ------------------------------------------------------------------------------------------------------ */
            
            /* Search for violated constraints. If any, then get the lowest scalar value _alpha that makes 
             *  the searching direction feasible again
             */
            float_prec _alphaMin = 1.0;
            int16_t _idxAlphaMin = -1;
            
            Matrix _ineqNonActiveLHS_dx(ineqLHS * dx);
            Matrix _ineqNonActiveLHS_x(ineqLHS * x);
            
            for (int16_t _i = 0; _i < ineqLHS.i16getRow(); _i++) {
                if (_flagConstActive[_i] == false) {
                    /* For all constraints that is not active, we check if that constraint is violated */
                    if (_ineqNonActiveLHS_dx[_i][0] > float_prec_ZERO_ECO) {
                        /* We have a violated constraint. Make that constraint a candidate of active constraint */
                        float_prec _alphaCandidate = - (_ineqNonActiveLHS_x[_i][0] - ineqRHS[_i][0]) / _ineqNonActiveLHS_dx[_i][0];
                        
                        if (_alphaCandidate < _alphaMin) {
                            _alphaMin = _alphaCandidate;
                            _idxAlphaMin = _i;
                        }
                    }
                }
            }
            if (_idxAlphaMin != -1) {
                _flagConstActive[_idxAlphaMin] = true;
            }
            
            /* Should we do ASSERT(_alphaMin > 0.0) here? */
            if (_alphaMin > 0) {
                x = x + (dx * _alphaMin);
            }
        }
        
        
        /* Check for maximum iteration conditions ------------------------------------------------------------------------------------------------- */
        _i16iterActiveSet++;
        if (_i16iterActiveSet > MPC_MAXIMUM_ACTIVE_SET_ITERATION) {
            /* Maximum number of iteration reached, terminate and use last solution */
            break;
        }
    } while (1);
    
    return true;
}


