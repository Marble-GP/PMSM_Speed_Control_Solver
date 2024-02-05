/**
 * @file pmsm_speed_control_solver.h
 * @author S. Watanabe (marble321done@gmail.com)
 * 
 * @brief This module comtains PMSM id&iq target value solver algorithm using speed target value (elec. rad/s) and actual voltage amp.
 *        This algorithm provides id&iq target values in MTPA control or flux weakening control.
 *        The comtrol styles of MTPA control and flux weakening selected automatically determined by calculated motor voltage amp.
 *        
 * @version 0.1
 * @date 2024-01-31
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PMSM_CONTROL_SOLVER_H
#define PMSM_CONTROL_SOLVER_H

#include <math.h>

#define PSCS_MTPA_MAX_ITERATION (10)
#define PSCS_MTPA_CONVERGENCE_THRESHOLD (1e-2f)
#define PSCS_MTPA_ID0_TORQUE_RATIO_THRESHOLD (1e-2f)
#define PSCS_MTPA_ID0_I_THRESHOLD (1e-1f)
#define PSCS_FW_MAX_ITERATION (20)
#define PSCS_FW_CONVERGENCE_THRESHOLD (1e-2f)

#define PSCS_SIGN(x) ((x) >= 0.0f ? 1.0f : -1.0f)
#define PSCS_ABS(x) ((x) > 0.0f ? (x) : -(x))

/**
 * @brief This is solution structure for speed control solver algorithm
 */
typedef struct
{

    float id_ref;
    float iq_ref;
    float Ia_ref;
    float beta_ref;
    float vd_calc;
    float vq_calc;
    float Va_calc;
    int FW_flag;

} PSCS_PMSM_Solutions_t;

/**
 * @brief This is PMSM parameters structure that be used by speed control solver algorithm
 */
typedef struct
{

    float Rs;
    float Ld;
    float Lq;
    float Psi_a;
    float Poles;

} PSCS_PMSM_Params_t;

/**
 * @brief This is PMSM speed control solver algorithm structure
 */
typedef struct
{
    PSCS_PMSM_Params_t* pmsm;

    float Va_lim;
    float Ia_lim;
    float Kc_mtpa;

} PSCS_Condition_t;

/**
 * @brief Initialize condition variables for PMSM speed control solver.
 * 
 * @param condition pointer to a PSCS_Condition_t variable to be initialized
 * @param motor_params pointer to initialized PSCS_PMSM_Params_t variable
 * @param Va_lim limit value of voltage amplitude 
 * @param Ia_lim limit value of current amplitude
 */
void PSCS_Init(PSCS_Condition_t* condition, PSCS_PMSM_Params_t* motor_params, float Va_lim, float Ia_lim);

/**
 * @brief  Initialize PMSM parameter variables for PMSM speed control solver.
 * 
 * @param params  pointer to a PSCS_PMSM_Params_t variable to be initialized
 * @param Rs stator registance
 * @param Ld d-axis motor inductance
 * @param Lq q-axis motor inductance
 * @param Psi_a amplitude of magnetic flux linkage from field
 * @param Poles number of poles
 */
void PSCS_PMSM_Params_Init(PSCS_PMSM_Params_t* params, float Rs, float Ld, float Lq, float Psi_a, float Poles);

/**
 * @brief Calculate PMSM speed control solver
 * 
 * @param condition pointer to initialized PSCS_Condition_t variable
 * @param solution pointer to a PSCS_Solution_t variable, its members are assigned the solutions
 * @param torque_ref Target value of torque (Nm)
 * @param w_ref Target value of electric angular velocity (rad/s)
 * @return int Returns 0 if solved succesfully. Returns -1 if not converged.
 */
int PSCS_Calculate(PSCS_Condition_t* condition, PSCS_PMSM_Solutions_t* solution, float torque_ref, float w_ref);

#endif