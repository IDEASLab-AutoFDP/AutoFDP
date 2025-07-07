#include "layout.h"

void AutoFDP(graph &g, int t_max, int choose_metric)
{
    g.solveDij();
    g.preSolveSGD(0.1, t_max, 42, 10);
    float maxt = 20, mint = 5;
    float theta_decay = pow((0.2), 1 / (float)(t_max));

    float w1 = 0.5, w2 = 0.5;

    switch (choose_metric)
    {
    case 0:
        fprintf(stderr, "opti SE\n");
        g.b2 = -1;
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5);
        }
        break;

    case 1:
        fprintf(stderr, "opti IL\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_IL(iter, 0.05, 1);
        }
        break;

    case 2:
        fprintf(stderr, "opti NP\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_NP(iter, 0.5 * pow(theta_decay, iter), 2, 20);
        }
        break;

    case 3:
        fprintf(stderr, "opti CL\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_CL(iter, 0.05, 10);
        }
        break;

    case 4:
        fprintf(stderr, "opti CA\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_CA(iter, 0.1 * pow(theta_decay, iter), 10);
        }
        break;

    case 5:
        fprintf(stderr, "opti AR\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_AR(iter, 0.2 * pow(theta_decay, iter), 10);
        }
        break;

    case 6:
        fprintf(stderr, "opti MA\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_MA(iter, 0.1 , 10);
        }
        break;

    case 7:
        fprintf(stderr, "opti NR\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_NR(iter, 0.2 , 20);
        }
        break;

    case 8:
        fprintf(stderr, "opti GP\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_GP(iter, 0.1 * pow(theta_decay, iter) , 20);
        }
        break;

    case 9:
        fprintf(stderr, "opti SE+IL\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_IL_mult(iter, 0.2 * pow(theta_decay, iter), 10 * w2) * w2;
        }
        break;

    case 10:
        fprintf(stderr, "opti SE+NP\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_NP(iter, 0.5 * pow(theta_decay, iter), 2, 20 * w2) * w2;
        }
        break;

    case 11:
        fprintf(stderr, "opti SE+CL\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_CL(iter, 0.05, 10 * w2) * w2;
        }
        break;

    case 12:
        fprintf(stderr, "opti SE+CA\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_CA(iter, 0.1 * pow(theta_decay, iter), 10 * w2) * w2;
        }
        break;

    case 13:
        fprintf(stderr, "opti SE+AR\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_AR(iter, 0.2 * pow(theta_decay, iter), 10 * w2) * w2;
        }
        break;
    case 14:
        fprintf(stderr, "opti SE+NR\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_NR(iter, 0.2, 20 * w2) * w2;
        }
        break;
    case 15:
        fprintf(stderr, "opti SE+GP\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_GP(iter, 0.1 * pow(theta_decay, iter), 20 * w2) * w2;
        }
        break;
    
    case 16:
        fprintf(stderr, "opti NP+CL\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_NP(iter, 0.5 * pow(theta_decay, iter), 2, 20 * w2) * w2;
            g.opti_CL(iter, 0.05, 10 * w2) * w2;
        }
        break;

    // DeepGD
    case 17:
        fprintf(stderr, "opti SE+NR\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_NR_DeepGD(iter, 0.2, 20 * w2) * w2;
        }
    break;

    case 18:
        fprintf(stderr, "opti SE+MA\n");
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_MA_DeepGD(iter, 0.1, 10 * w2) * w2;

        }
    break;

    case 19:
        fprintf(stderr, "opti SE+TSNE\n");
        g.calP_DeepGD();
        for (int iter = 0; iter < t_max; iter++)
        {
            g.solve_bilevel(iter);
            g.opti_SE(iter, 0.05, 5 * w1) * w1;
            g.opti_tSNE_DeepGD(iter, 0.05, 5 * w2) * w2;
        }
    break;

    case -1:
        fprintf(stderr, "opti None\n");
        for (int iter = 0; iter < t_max; iter++)
            g.solve_bilevel(iter);
    }
}