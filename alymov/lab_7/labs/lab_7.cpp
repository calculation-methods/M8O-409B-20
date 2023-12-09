#include "lab_7.hpp"

Matrix libman(float acc=0.01){

    //init soliton and temporary variables
    Matrix sol(M, N);
    Matrix prev_sol(M, N);
    float max_difference = 10000.0;
    float dif = 0;
    float min = 0;

    //init boundary conditions
    for(int k = 0; k < M; k++){
        prev_sol[k][0] = U_X_0(0, H_Y * k);
        prev_sol[k][M - 1] = U_X_1(1, H_Y * k);
    }

    for(int k = 0; k < N; k++){
        prev_sol[0][k] = U_Y_0(H_X * k, 0);
        prev_sol[N-1][k] = U_Y_1(H_X * k, 1);
    }

    //solve
    for(int k=0; max_difference > acc; k++){
        min = 0;
        for(int i = 1; i < M - 1; i++){
            for(int j = 1; j < N - 1; j++){
                sol[i][j] = 0.25 * (prev_sol[i][j - 1] + prev_sol[i][j + 1] \
                + prev_sol[i - 1][j] + prev_sol[i + 1][j]);
                
                dif = fabs(sol[i][j] - prev_sol[i][j]);
                max_difference = (dif < max_difference && dif > 0) ? fabs(sol[i][j] - prev_sol[i][j]) : max_difference;
            }
        }
        for(int c = 0; c < M; c++){
            sol[c][0] = U_X_0(0, H_Y * c);
            sol[c][M - 1] = U_X_1(1, H_Y * c);
        }

        for(int c = 0; c < N; c++){
            sol[0][c] = U_Y_0(H_X * c, 0);
            sol[N-1][c] = U_Y_1(H_X * c, 1);
        }
        prev_sol = sol;
    }

    return sol;
}