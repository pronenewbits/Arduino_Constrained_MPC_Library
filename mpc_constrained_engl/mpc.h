/**************************************************************************************************
 * Class for MPC without constraint.
 * 
 * 
 * See https://github.com/pronenewbits for more!
 *************************************************************************************************/
#ifndef MPC_H
#define MPC_H

#include "konfig.h"
#include "matrix.h"


#if (MPC_HP_LEN < MPC_HU_LEN)
    #error("The MPC_HP_LEN must be more than or equal MPC_HU_LEN!");
#endif
#if (((MPC_HP_LEN*SS_Z_LEN) > MATRIX_MAXIMUM_SIZE) || ((MPC_HP_LEN*SS_X_LEN) > MATRIX_MAXIMUM_SIZE) || ((MPC_HU_LEN*SS_U_LEN) > MATRIX_MAXIMUM_SIZE))
    #error("The MATRIX_MAXIMUM_SIZE is too small to do MPC calculation!");
#endif

#define MPC_MAXIMUM_ACTIVE_SET_ITERATION    (50)

#if defined(MPC_USE_CONSTRAINT_DU)
    #define INEQ_LEN_DU     (2*(MPC_HU_LEN*SS_U_LEN))       /* size DU = (MPC_HU_LEN*SS_U_LEN), times 2 for DU_MIN & DU_MAX */
#else
    #define INEQ_LEN_DU     (0)
#endif
#if defined(MPC_USE_CONSTRAINT_U)
    #define INEQ_LEN_U      (2*(MPC_HU_LEN*SS_U_LEN))       /* size U  = (MPC_HU_LEN*SS_U_LEN), times 2 for U_MIN & U_MAX */
#else
    #define INEQ_LEN_U      (0)
#endif
#if defined(MPC_USE_CONSTRAINT_Z)
    #define INEQ_LEN_Z      (2*(MPC_HP_LEN*SS_Z_LEN))       /* size Z  = (MPC_HP_LEN*SS_Z_LEN), times 2 for Z_MIN & Z_MAX */
#else
    #define INEQ_LEN_Z      (0)
#endif
#if (( (MPC_HP_LEN*SS_Z_LEN) + INEQ_LEN_DU + INEQ_LEN_U + INEQ_LEN_Z) > MATRIX_MAXIMUM_SIZE)
    /* Worst case for KKT Matrix size: Q + Cnst_DU + Cnst_U + Cnst_Z */
//     #error("The MATRIX_MAXIMUM_SIZE is too small to do MPC calculation!");
// WARNING! we disable this because Teensy has too little RAM if we enable this to do anything non-trivia, you __HAVE__ to enable this for mission critical implementation.
#endif



class MPC
{
public:
    MPC(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR);
    void vReInit(Matrix &A, Matrix &B, Matrix &C, float_prec _bobotQ, float_prec _bobotR);
    bool bUpdate(Matrix &SP, Matrix &x, Matrix &u);
    
    /* for analyzing purpose. If cntIterActiveSet = MPC_MAXIMUM_ACTIVE_SET_ITERATION, the solution might not be optimal (although feasible) */
    int16_t cntIterActiveSet;
protected:
    bool bActiveSet(Matrix &x, Matrix &Q, Matrix &c, Matrix &ineqLHS, Matrix &ineqRHS, int16_t &_i32iterActiveSet);

private:
    Matrix CPSI     {(MPC_HP_LEN*SS_Z_LEN), SS_X_LEN};
    Matrix COMEGA   {(MPC_HP_LEN*SS_Z_LEN), SS_U_LEN};
    Matrix CTHETA   {(MPC_HP_LEN*SS_Z_LEN), (MPC_HU_LEN*SS_U_LEN)};

    Matrix DU       {(MPC_HU_LEN*SS_U_LEN), 1};

    Matrix A        {SS_X_LEN, SS_X_LEN};
    Matrix B        {SS_X_LEN, SS_U_LEN};
    Matrix C        {SS_Z_LEN, SS_X_LEN};

    Matrix Q        {(MPC_HP_LEN*SS_Z_LEN), (MPC_HP_LEN*SS_Z_LEN)};
    Matrix R        {(MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN)};
    
    Matrix H        {(MPC_HU_LEN*SS_U_LEN), (MPC_HU_LEN*SS_U_LEN)};
};



#endif // MPC_H
