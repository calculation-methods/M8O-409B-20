#include "lab_8.hpp"

std::vector<Matrix> alterning_direction_method(){
    std::vector<Matrix> sol;
    
    Matrix u(M, N);
    Matrix u_prev = init();

    sol.push_back(u_prev);

    for(int k = 1;k < K; k++){
        u = step(u_prev, k);
        sol.push_back(u);
        u_prev = u;
    }
    return sol;
}


Matrix init()
{
    Matrix initial(M, N);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            initial[i][j] = U_XY(i * H_X, j * H_Y);
        }
    }
    return initial;
}

Matrix step(Matrix& prev, int k)
{
    Matrix u(prev.cols, prev.rows);
    u = horizontal(prev, k);
    u = vertical(u, k);
    return u;
}


Matrix horizontal(Matrix& prev_sol, float k)
{
    Matrix sol(prev_sol.cols, prev_sol.rows);
    Vector row_sol(prev_sol.cols);

    Matrix A(prev_sol.cols, prev_sol.cols);
    Vector b(prev_sol.cols);

    float y;
    float t = k * TAU - TAU / 2;

    for(int i = 1; i < prev_sol.rows-1; i++){
        Vector prev_sol_row = prev_sol.get_row(i);
        y = i * H_Y; 
        A[0][0] = 1;
        b[0] = U_X_0(0, y, t);

        for(int j = 1; j < A.cols - 1; j++){
            A[j][j - 1] = ALPHA_SQ * TAU / (H_X * H_X);
            A[j][j] = - 1 - 2 * ALPHA_SQ * TAU / (H_X * H_X); 
            A[j][j + 1] =  ALPHA_SQ * TAU / (H_X * H_X);
            b[j] = -1 * prev_sol_row[j];
        }
        A[A.cols - 1][A.cols - 1] = 1;
        b[A.cols - 1] = U_X_1(L, y, t);
        
        row_sol = Tridiagonal_matrix_algorithm(A, b);
        embed_vector(sol, row_sol, false, i);
    }

    return sol;
}

Matrix vertical(Matrix& prev_sol, float k)
{
    Matrix sol(prev_sol.cols, prev_sol.rows);
    Vector col_sol(prev_sol.rows);

    Matrix A(prev_sol.rows, prev_sol.rows);
    Vector b(prev_sol.rows);

    float x;
    float t = k * TAU;

    for(int i = 1; i < prev_sol.cols - 1; i++){
        Vector prev_sol_col = prev_sol.get_col(i);
        x = i * H_X;

        A[0][0] = 1;
        b[0] = U_Y_0(x, 0, t);

        for(int j = 1; j < A.rows - 1; j++){
            A[j][j - 1] = ALPHA_SQ * TAU / (H_Y * H_Y);
            A[j][j] = - 1 - 2 * ALPHA_SQ * TAU / (H_Y * H_Y); 
            A[j][j + 1] =  ALPHA_SQ * TAU / (H_Y * H_Y);
            b[j] = -1 * prev_sol_col[j];
        }

        A[A.rows - 1][A.rows - 1] = 1;
        b[A.rows - 1] = U_Y_1(x, L, t);
        col_sol = Tridiagonal_matrix_algorithm(A, b);
        embed_vector(sol, col_sol, true, i);
    }
    
    for(int j = 0; j < prev_sol.rows; j++){
        sol[0][j] = U_Y_0(j*H_X, 0, t);
        sol[sol.rows-1][j] = U_Y_1(j*H_X, L, t);
    }

    for(int j = 0; j < prev_sol.cols; j++){
        sol[j][0] = U_X_0(0, j*H_Y, t);
        sol[j][sol.cols-1] = U_X_1(L, j*H_Y, t);
    }

    return sol;
}

void embed_vector(Matrix& A, Vector& b, bool cols, int ind)
{
    if((A.rows != b.len && ind > A.rows && cols == true) || (A.cols != b.len && ind > A.rows))
        return;

    if(cols == true){
        for(int i  = 0; i < A.cols; i++)
            A[i][ind] = b[i];
    } else {
        for(int i  = 0; i < A.rows; i++)
            A[ind][i] = b[i];
    }
    return;
}