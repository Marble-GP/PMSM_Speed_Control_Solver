/**
 * @file pmsm_speed_control_solver.c
 * @author S. Watanabe (marble321done@gmail.com)
 *      
 * @version 0.1
 * @date 2024-01-31
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "pmsm_speed_control_solver.h"


void PSCS_Init(PSCS_Condition_t* condition, PSCS_PMSM_Params_t* motor_params, float Va_lim, float Ia_lim)
{
    condition->pmsm = motor_params;
    condition->Va_lim = PSCS_ABS(Va_lim);
    condition->Ia_lim = PSCS_ABS(Ia_lim);
    condition->Kc_mtpa = 1.0f;
}

void PSCS_PMSM_Params_Init(PSCS_PMSM_Params_t* params, float Rs, float Ld, float Lq, float Psi_a, float Poles)
{
    params->Rs = Rs;
    params->Ld = Ld;
    params->Lq = Lq;
    params->Psi_a = Psi_a;
    params->Poles = Poles;
}

int PSCS_Calculate(PSCS_Condition_t* condition, PSCS_PMSM_Solutions_t* solution, float torque_ref, float w_ref)
{
    int iter = 0, errcode = -1;
    float torque_calc;

    solution->Ia_ref = condition->Kc_mtpa * torque_ref / condition->pmsm->Psi_a / condition->pmsm->Poles * 2.0f;
    solution->beta_ref = 0.0f;

    if(PSCS_ABS(solution->Ia_ref) > PSCS_MTPA_ID0_I_THRESHOLD && 0.5f * PSCS_ABS(condition->pmsm->Lq - condition->pmsm->Ld) * solution->Ia_ref / condition->pmsm->Psi_a > PSCS_MTPA_ID0_TORQUE_RATIO_THRESHOLD ) // when it's worth to calculate MTPA current equation
    {
        do
        {
            solution->beta_ref = asinf( (-condition->pmsm->Psi_a + sqrtf(condition->pmsm->Psi_a * condition->pmsm->Psi_a + 8.0f*(condition->pmsm->Lq - condition->pmsm->Ld) * (condition->pmsm->Lq - condition->pmsm->Ld) * solution->Ia_ref * solution->Ia_ref)) / (4.0f * (condition->pmsm->Lq - condition->pmsm->Ld) * solution->Ia_ref));

            torque_calc =  0.5f * condition->pmsm->Poles * ( condition->pmsm->Psi_a*solution->Ia_ref*cosf(solution->beta_ref) + 0.5f*(condition->pmsm->Lq - condition->pmsm->Ld)*solution->Ia_ref*solution->Ia_ref*sinf(2.0f*solution->beta_ref) );

            solution->Ia_ref *= torque_ref/torque_calc;

            ++iter;

        } while (iter < PSCS_MTPA_MAX_ITERATION && PSCS_ABS((torque_calc - torque_ref)/torque_ref) > PSCS_MTPA_CONVERGENCE_THRESHOLD);

        solution->id_ref = -solution->Ia_ref * sinf(solution->beta_ref);
		solution->iq_ref = solution->Ia_ref * cosf(solution->beta_ref);

        solution->vd_calc = condition->pmsm->Rs * solution->id_ref - w_ref * condition->pmsm->Lq * solution->iq_ref;
        solution->vq_calc = condition->pmsm->Rs * solution->iq_ref + w_ref * condition->pmsm->Ld * solution->id_ref + w_ref * condition->pmsm->Psi_a;
        solution->Va_calc = sqrtf(solution->vd_calc * solution->vd_calc  +  solution->vq_calc * solution->vq_calc);

        errcode = -(iter == PSCS_MTPA_MAX_ITERATION);
    }
    else //assume id=0 control is equivalent to MTPA control
    {
        solution->Ia_ref = solution->iq_ref = torque_ref / condition->pmsm->Psi_a / condition->pmsm->Poles * 2.0f;
		solution->id_ref = 0.0f;
        solution->beta_ref = 0.0f;
		torque_calc = 0.5f * condition->pmsm->Psi_a * condition->pmsm->Poles * solution->Ia_ref;

        solution->vd_calc = - w_ref * condition->pmsm->Lq * solution->iq_ref;
        solution->vq_calc = condition->pmsm->Rs * solution->iq_ref + w_ref * condition->pmsm->Psi_a;
        solution->Va_calc = sqrtf(solution->vd_calc * solution->vd_calc  +  solution->vq_calc * solution->vq_calc);

        errcode = 0;
    }

    solution->FW_flag = solution->Va_calc > condition->Va_lim;

    if(solution->FW_flag)
    {
        float id_ref_pre;
        solution->iq_ref = torque_ref / condition->pmsm->Poles * 2.0f;
        iter = 0;
        
        do
        {
            id_ref_pre = solution->id_ref;

			solution->id_ref = condition->Va_lim > (condition->pmsm->Ld * condition->pmsm->Lq * w_ref * w_ref * solution->iq_ref + condition->pmsm->Rs * condition->pmsm->Rs * solution->iq_ref + condition->pmsm->Rs * condition->pmsm->Psi_a * w_ref) * sqrtf(1.0f / (condition->pmsm->Ld * condition->pmsm->Ld * w_ref * w_ref + condition->pmsm->Rs * condition->pmsm->Rs))
            
            ?   (-w_ref * (condition->pmsm->Ld * condition->pmsm->Rs * solution->iq_ref + condition->pmsm->Ld * condition->pmsm->Psi_a * w_ref - condition->pmsm->Lq * condition->pmsm->Rs * solution->iq_ref)
                + sqrtf(   -condition->pmsm->Ld * condition->pmsm->Ld * condition->pmsm->Lq * condition->pmsm->Lq * w_ref * w_ref * w_ref * w_ref * solution->iq_ref * solution->iq_ref + condition->pmsm->Ld * condition->pmsm->Ld * condition->Va_lim * condition->Va_lim * w_ref * w_ref - 2.0f * condition->pmsm->Ld * condition->pmsm->Lq * condition->pmsm->Rs * condition->pmsm->Rs * w_ref * w_ref * solution->iq_ref * solution->iq_ref - 2.0f * condition->pmsm->Ld * condition->pmsm->Lq * condition->pmsm->Rs * condition->pmsm->Psi_a * w_ref * w_ref * w_ref * solution->iq_ref - condition->pmsm->Rs * condition->pmsm->Rs * condition->pmsm->Rs * condition->pmsm->Rs * solution->iq_ref * solution->iq_ref
                            - 2.0f * condition->pmsm->Rs * condition->pmsm->Rs * condition->pmsm->Rs * condition->pmsm->Psi_a * w_ref * solution->iq_ref + condition->pmsm->Rs * condition->pmsm->Rs * condition->Va_lim * condition->Va_lim - condition->pmsm->Rs * condition->pmsm->Rs * condition->pmsm->Psi_a * condition->pmsm->Psi_a * w_ref * w_ref))
                / (condition->pmsm->Ld * condition->pmsm->Ld * w_ref * w_ref + condition->pmsm->Rs * condition->pmsm->Rs)

            : -w_ref * (condition->pmsm->Ld * condition->pmsm->Rs * solution->iq_ref + w_ref * condition->pmsm->Ld * condition->pmsm->Psi_a - condition->pmsm->Lq * condition->pmsm->Rs * solution->iq_ref) / (w_ref * w_ref * condition->pmsm->Ld * condition->pmsm->Ld + condition->pmsm->Rs * condition->pmsm->Rs);
			
            // easier id equation
            // solution->id_ref = Vref_lim > w_ref * condition->pmsm->Lq * solution->iq_ref ? (-condition->pmsm->Psi_a + sqrtf((Vref_lim / w_ref) * (Vref_lim / w_ref) - condition->pmsm->Lq * condition->pmsm->Lq * solution->iq_ref * solution->iq_ref)) / condition->pmsm->Ld : -condition->pmsm->Psi_a / condition->pmsm->Ld;
			
            solution->iq_ref = 2.0f * torque_ref / condition->pmsm->Poles / (condition->pmsm->Psi_a + (condition->pmsm->Ld - condition->pmsm->Lq) * solution->id_ref); // solve torque equation at iq
			
            if (sqrtf(solution->id_ref * solution->id_ref + solution->iq_ref * solution->iq_ref) > condition->Ia_lim) // when over the current limiter
			{
				solution->iq_ref = condition->Ia_lim * condition->Ia_lim > solution->id_ref * solution->id_ref ? PSCS_SIGN(torque_ref) * sqrtf(condition->Ia_lim * condition->Ia_lim - solution->id_ref * solution->id_ref) : 0.0f;
			}

			++iter;

        } while (iter < PSCS_FW_MAX_ITERATION && PSCS_ABS((solution->id_ref - id_ref_pre) / solution->id_ref) > PSCS_FW_CONVERGENCE_THRESHOLD);
        

        solution->Ia_ref = sqrtf(solution->id_ref * solution->id_ref + solution->iq_ref * solution->iq_ref);
        solution->beta_ref = atan2f(solution->iq_ref, -solution->id_ref);

        solution->vd_calc = condition->pmsm->Rs * solution->id_ref - w_ref * condition->pmsm->Lq * solution->iq_ref;
        solution->vq_calc = condition->pmsm->Rs * solution->iq_ref + w_ref * condition->pmsm->Ld * solution->id_ref + w_ref * condition->pmsm->Psi_a;
        solution->Va_calc = sqrtf(solution->vd_calc * solution->vd_calc  +  solution->vq_calc * solution->vq_calc);

        errcode = -(iter == PSCS_FW_MAX_ITERATION);

    }
    

    return errcode;
}