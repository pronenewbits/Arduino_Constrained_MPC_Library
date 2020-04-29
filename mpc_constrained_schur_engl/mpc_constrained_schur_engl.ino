#include <Wire.h>
#include <elapsedMillis.h>
#include "konfig.h"
#include "matrix.h"
#include "mpc.h"


void vCreateConstraintsLHS(Matrix &CnstLHS, const Matrix &CTHETA);
void vCreateConstraintsRHS(Matrix &CnstRHS, const Matrix &COMEGA, const Matrix &CPSI, const Matrix &u_prev, const Matrix &x);


elapsedMillis timerMPC;     /* Timer for sampling time */
uint64_t u64compuTime;      /* For benchmark */
char bufferTxSer[100];      /* For serial printing */


/* Plant system */
Matrix A(SS_X_LEN, SS_X_LEN);
Matrix B(SS_X_LEN, SS_U_LEN);
Matrix C(SS_Z_LEN, SS_X_LEN);
Matrix SP((MPC_HP_LEN*SS_Z_LEN), 1);
Matrix x(SS_X_LEN, 1);
Matrix u(SS_U_LEN, 1);
Matrix z(SS_Z_LEN, 1);

int16_t i16iterSP = 0;
MPC MPC_HIL(A, B, C, 1, 0.01, vCreateConstraintsLHS, vCreateConstraintsRHS);
uint64_t _maxu64compuTime = 0;
void setup() {
    /* serial to display data */
    Serial.begin(115200);
    while(!Serial) {}
    pinMode(LED_BUILTIN, OUTPUT);
    
    
    /* Ref: https://www.mathworks.com/help/control/ug/mimo-state-space-models.html#buv3tp8-1
     * 
     * State-Space Model of Jet Transport Aircraft
     *  This example shows how to build a MIMO model of a jet transport. Because the development of a physical model 
     *  for a jet aircraft is lengthy, only the state-space equations are presented here. See any standard text in 
     *  aviation for a more complete discussion of the physics behind aircraft flight.
     * The jet model during cruise flight at MACH = 0.8 and H = 40,000 ft. is
     * 
     * (The model has two inputs and two outputs. The units are radians for beta (sideslip angle) and phi (bank angle) and 
     * radians/sec for yaw (yaw rate) and roll (roll rate). The rudder and aileron deflections are in degrees.)
     */ 
    A[0][0] = -0.0558;      A[0][1] = -0.9968;      A[0][2] =  0.0802;      A[0][3] = 0.0415;
    A[1][0] =  0.5980;      A[1][1] = -0.1150;      A[1][2] = -0.0318;      A[1][3] = 0.0000;
    A[2][0] = -3.0500;      A[2][1] =  0.3880;      A[2][2] = -0.4650;      A[2][3] = 0.0000;
    A[3][0] =  0.0000;      A[3][1] =  0.0805;      A[3][2] =  1.0000;      A[3][3] = 0.0000;
    
    B[0][0] =  0.0073;      B[0][1] =  0.0000;
    B[1][0] = -0.4750;      B[1][1] =  0.0077;
    B[2][0] =  0.1530;      B[2][1] =  0.1430;
    B[3][0] =  0.0000;      B[3][1] =  0.0000;
    
    C[0][0] =  0.0000;      C[0][1] =  1.0000;      C[0][2] =  0.0000;      C[0][3] =  0.0000;
    C[1][0] =  0.0000;      C[1][1] =  0.0000;      C[1][2] =  0.0000;      C[1][3] =  1.0000;
    
    MPC_HIL.vReInit(A, B, C, 10, 0.05, vCreateConstraintsLHS, vCreateConstraintsRHS);
}



void loop() {
    if (timerMPC > SS_DT_MILIS) {
        /* ================================ Updating Set Point ================================= */
        Matrix SP_NEXT(SS_Z_LEN, 1);
        if (i16iterSP < 100-MPC_HP_LEN+1) {
            SP_NEXT[0][0] = 3.14/2.;
            SP_NEXT[1][0] = 1;
        } else if (i16iterSP < 200-MPC_HP_LEN+1) {
            SP_NEXT[0][0] = 3.14/2.;
            SP_NEXT[1][0] = -3;
        } else {
            SP_NEXT[0][0] = 3.14;
            SP_NEXT[1][0] = -3;
        }
        if (i16iterSP < 300-MPC_HP_LEN+1) {
            i16iterSP++;
        } else {
            i16iterSP = 0;
        }
        for (int32_t _i = 0; _i < (MPC_HP_LEN-1); _i++) {
            SP = SP.InsertSubMatrix(SP, (_i*SS_Z_LEN), 0, ((_i+1)*SS_Z_LEN), 0, SS_Z_LEN, 1);
        }
        SP = SP.InsertSubMatrix(SP_NEXT, ((MPC_HP_LEN-1)*SS_Z_LEN), 0, 0, 0, SS_Z_LEN, 1);
        /* -------------------------------- Updating Set Point --------------------------------- */
        
        
        digitalWrite(LED_BUILTIN, HIGH);
        /* ===================================== MPC Update ==================================== */
        u64compuTime = micros();
        
        MPC_HIL.bUpdate(SP, x, u);
        
        u64compuTime = (micros() - u64compuTime);        
        
        _maxu64compuTime < u64compuTime? _maxu64compuTime = u64compuTime: 0;
        /* ------------------------------------- MPC Update ------------------------------------ */
        digitalWrite(LED_BUILTIN, LOW);
        
        
        
        /* ================================= Plant Simulation ================================== */
        x = A*x + B*u;
        z = C*x;
        /* --------------------------------- Plant Simulation ---------------------------------- */
        
        
        
        /* =========================== Print to serial (for plotting) ========================== */
        #if (1)
            /* Print: Computation time, number of active set iteration, Set-Point, z */
            snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "%.3f %i %.3f %.3f %.3f %.3f", ((float)_maxu64compuTime)/1000., MPC_HIL.cntIterActiveSet, SP[0][0], SP[1][0], z[0][0], z[1][0]);
        #else
            /* Print: Computation time, Set-Point, z, u */
            snprintf(bufferTxSer, sizeof(bufferTxSer)-1, "%.3f %i %.3f %.3f %.3f %.3f %.3f %.3f", ((float)_maxu64compuTime)/1000., MPC_HIL.cntIterActiveSet, SP[0][0], SP[1][0], z[0][0], z[1][0], u[0][0], u[1][0]);
        #endif
        Serial.print(bufferTxSer);
        Serial.print('\n');
        /* --------------------------- Print to serial (for plotting) -------------------------- */
        
        
        timerMPC = 0;
    }
}


void vCreateConstraintsLHS(Matrix &CnstLHS, const Matrix &CTHETA)
{
    #if (defined(MPC_USE_CONSTRAINT_DU) || defined(MPC_USE_CONSTRAINT_U) || defined(MPC_USE_CONSTRAINT_Z))
        int16_t _i16offsetCnstr = 0;
        CnstLHS = Matrix((INEQ_LEN_DU + INEQ_LEN_U + INEQ_LEN_Z), (MPC_HU_LEN*SS_U_LEN));
    #endif
    #if defined(MPC_USE_CONSTRAINT_DU)
        /* for dUmin <= dU(kj <= dUmax                                                      ...{MPC_5a} */
        {
            Matrix CnstDu_LHS   {(2*(MPC_HU_LEN*SS_U_LEN)), (MPC_HU_LEN*SS_U_LEN)};
            Matrix CnstDu_Identity((MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN));
            CnstDu_Identity.vSetIdentity();
            CnstDu_LHS = CnstDu_LHS.InsertSubMatrix(-CnstDu_Identity, 0, 0);
            CnstDu_LHS = CnstDu_LHS.InsertSubMatrix(CnstDu_Identity, (MPC_HU_LEN*SS_U_LEN), 0);

            /* Append to Constraint_LHS & Constraint_RHS */
            CnstLHS = CnstLHS.InsertSubMatrix(CnstDu_LHS, _i16offsetCnstr, 0);
            _i16offsetCnstr += CnstDu_LHS.i16getRow();
        }
    #endif
    #if defined(MPC_USE_CONSTRAINT_U)
        /* for Umin <= U(kj <= Umax                                                         ...{MPC_5b} */
        {
            Matrix CnstU_LHS   {(2*(MPC_HU_LEN*SS_U_LEN)), (MPC_HU_LEN*SS_U_LEN)};
            Matrix U_Identity(SS_U_LEN, SS_U_LEN);
            U_Identity.vSetIdentity();
            for (int16_t _i = 0; _i < MPC_HU_LEN; _i++) {
                for (int16_t _j = 0; _j < MPC_HU_LEN; _j++) {
                    if (_j <= _i) {
                        /* for Umin <= U(k) */
                        CnstU_LHS = CnstU_LHS.InsertSubMatrix(-U_Identity, (_i*SS_U_LEN), (_j*SS_U_LEN));
                        /* for U(k) <= Umax */
                        CnstU_LHS = CnstU_LHS.InsertSubMatrix(U_Identity,  (MPC_HU_LEN*SS_U_LEN)+(_i*SS_U_LEN), (_j*SS_U_LEN));
                    }
                }
            }

            /* Append to Constraint_LHS & Constraint_RHS */
            CnstLHS = CnstLHS.InsertSubMatrix(CnstU_LHS, _i16offsetCnstr, 0);
            _i16offsetCnstr += CnstU_LHS.i16getRow();
        }
    #endif
    #if defined(MPC_USE_CONSTRAINT_Z)
        /* for Zmin <= Z(kj <= Zmax                                                         ...{MPC_5c} */
        {
            Matrix CnstZ_LHS   {(2*(MPC_HP_LEN*SS_Z_LEN)), (MPC_HU_LEN*SS_U_LEN)};
            CnstZ_LHS = CnstZ_LHS.InsertSubMatrix(-CTHETA, 0, 0);
            CnstZ_LHS = CnstZ_LHS.InsertSubMatrix(CTHETA, (MPC_HP_LEN*SS_Z_LEN), 0);

            /* Append to Constraint_LHS & Constraint_RHS */
            CnstLHS = CnstLHS.InsertSubMatrix(CnstZ_LHS, _i16offsetCnstr, 0);
            _i16offsetCnstr += CnstZ_LHS.i16getRow();
        }
    #endif
}

#if defined(MPC_USE_CONSTRAINT_DU)
    Matrix CnstDu_RHS   {(2*(MPC_HU_LEN*SS_U_LEN)), 1};
    Matrix CnstDu_MaxMin((MPC_HU_LEN*SS_U_LEN), 1);
#endif
#if defined(MPC_USE_CONSTRAINT_U)
    Matrix CnstU_RHS   {(2*(MPC_HU_LEN*SS_U_LEN)), 1};
    Matrix CnstU_MaxMin((MPC_HU_LEN*SS_U_LEN), 1);
#endif
#if defined(MPC_USE_CONSTRAINT_Z)
    Matrix CnstZ_RHS   {(2*(MPC_HP_LEN*SS_Z_LEN)), 1};
    Matrix CnstZ_Min((MPC_HP_LEN*SS_Z_LEN), 1);
    Matrix CnstZ_Max((MPC_HP_LEN*SS_Z_LEN), 1);
#endif

void vCreateConstraintsRHS(Matrix &CnstRHS, const Matrix &COMEGA, const Matrix &CPSI, const Matrix &u_prev, const Matrix &x)
{
	#if (defined(MPC_USE_CONSTRAINT_DU) || defined(MPC_USE_CONSTRAINT_U) || defined(MPC_USE_CONSTRAINT_Z))
        int16_t _i16offsetCnstr = 0;
        CnstRHS = Matrix((INEQ_LEN_DU + INEQ_LEN_U + INEQ_LEN_Z), 1);
    #endif
    #if defined(MPC_USE_CONSTRAINT_DU)
        /* for dUmin <= dU(kj <= dUmax                                                      ...{MPC_5a} */
        {
            CnstDu_MaxMin.vSetHomogen(MPC_MIN_DU);
            CnstDu_RHS = CnstDu_RHS.InsertSubMatrix(-CnstDu_MaxMin, 0, 0);
            CnstDu_MaxMin.vSetHomogen(MPC_MAX_DU);
            CnstDu_RHS = CnstDu_RHS.InsertSubMatrix(CnstDu_MaxMin, (MPC_HU_LEN*SS_U_LEN), 0);

            /* Append to Constraint_LHS & Constraint_RHS */
            CnstRHS = CnstRHS.InsertSubMatrix(CnstDu_RHS, _i16offsetCnstr, 0);
            _i16offsetCnstr += CnstDu_RHS.i16getRow();
        }
    #endif
    #if defined(MPC_USE_CONSTRAINT_U)
        /* for Umin <= U(kj <= Umax                                                         ...{MPC_5b} */
        {
            for (int16_t _i = 0; _i < MPC_HU_LEN; _i++) {
                for (int16_t _j = 0; _j < SS_U_LEN; _j++) {
                    CnstU_MaxMin[(_i*SS_U_LEN)+_j][0] = -MPC_MIN_U + u_prev[_j][0];
                }
            }
            CnstU_RHS = CnstU_RHS.InsertSubMatrix(CnstU_MaxMin, 0, 0);
            for (int16_t _i = 0; _i < MPC_HU_LEN; _i++) {
                for (int16_t _j = 0; _j < SS_U_LEN; _j++) {
                    CnstU_MaxMin[(_i*SS_U_LEN)+_j][0] = MPC_MAX_U - u_prev[_j][0];
                }
            }
            CnstU_RHS = CnstU_RHS.InsertSubMatrix(CnstU_MaxMin, (MPC_HU_LEN*SS_U_LEN), 0);

            /* Append to Constraint_LHS & Constraint_RHS */
            CnstRHS = CnstRHS.InsertSubMatrix(CnstU_RHS, _i16offsetCnstr, 0);
            _i16offsetCnstr += CnstU_RHS.i16getRow();
        }
    #endif
    #if defined(MPC_USE_CONSTRAINT_Z)
        /* for Zmin <= Z(kj <= Zmax                                                         ...{MPC_5c} */
        {
            CnstZ_Min = CPSI*x + COMEGA*u_prev - MPC_MIN_Z;
            CnstZ_RHS = CnstZ_RHS.InsertSubMatrix(CnstZ_Min, 0, 0);

            CnstZ_Max = -CPSI*x + -COMEGA*u_prev + MPC_MAX_Z;
            CnstZ_RHS = CnstZ_RHS.InsertSubMatrix(CnstZ_Max, (MPC_HP_LEN*SS_Z_LEN), 0);

            /* Append to Constraint_LHS & Constraint_RHS */
            CnstRHS = CnstRHS.InsertSubMatrix(CnstZ_RHS, _i16offsetCnstr, 0);
            _i16offsetCnstr += CnstZ_RHS.i16getRow();
        }
    #endif
}



void SPEW_THE_ERROR(char const * str)
{
    #if (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_PC)
        cout << (str) << endl;
    #elif (SYSTEM_IMPLEMENTATION == SYSTEM_IMPLEMENTATION_EMBEDDED_ARDUINO)
        Serial.println(str);
    #else
        /* Silent function */
    #endif
    while(1);
}
