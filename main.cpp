#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>



extern "C"
{
    double dnrm2_(const int *, const double *, const int *);
    void dgels_(const char *, const int *, const int *, const int *, double *, const int *, double *, const int *, double *, const int *, int *);
}


void
tensor_make(double** tensor, long long length)
{
    *tensor = new double[length];
    for (int i = 0; i < length; i++) {
        (*tensor)[i] = -1.0 + 2.0 * double(std::rand()) / RAND_MAX;
    }
}

double 
gen_func(int k, int *razm, int N)
{
    int *b = new int[N];
    int m = 0;
    b[N - 1] = k % razm[N - 1];
    k = (k - b[N - 1]) / razm[N - 1];
    m = b[N - 1];
    for (int i = N - 2; i >= 0; i--)
    {
        b[i] = k % razm[i];
        k = (k - b[i]) / razm[i];
        m = m + b[i];
    }
    return std::sin(m);
}
void
tensor_make_not_random(double** tensor, int *razm, int length, int N)
{
    *tensor = new double[length];
    for (int i = 0; i < length; i++) {
        (*tensor)[i] = gen_func(i, razm, N);
    }
}

void
create_matrix(double* &matrix, long long razm, int rank)
{
    matrix = new double[razm * rank];
    for(int r = 0; r < rank; r++) {
        for(int row = 0; row < razm; row++) {
            matrix[row + r * razm] = -1.0 + 2.0 * double(std::rand()) / RAND_MAX;
        }
    }
}


void
create_S_T(double** matrices, int busy_side, int rank, int N, int* razm, double* S, double* res, double* tensor)
{
    

    long long* multipliers = new long long[N + 1];

    multipliers[N] = 1;
    for(int i = N - 1; i >= 0; i--) {
        multipliers[i] = multipliers[i + 1] * razm[i + 1];
    }

    long long busy_multipliers = multipliers[busy_side];
    for(int i = busy_side; i < N; i++) {
        multipliers[i] = multipliers[i + 1];
    }

    long long* mas = new long long[N];
    long long* mas_razm = new long long[N];

    for(int i = 0; i < N + 1; i++) {
        if(i < busy_side) {
            mas_razm[i] = razm[i];
        } else if(i > busy_side) {
            mas_razm[i - 1] = razm[i];
        } else {
            continue;
        }
    }

    long long count = 0, tensor_val = 0;
    double val;

    long long max_r = std::max(rank, razm[busy_side]);

    for(long long r = 0; r < max_r; r++) {

        for(int pas = 0; pas < N; pas++) {
            mas[pas] = 0;
        }

        while(mas[0] < mas_razm[0]) {
            while(mas[N - 1] < mas_razm[N - 1]) {

                if(r < rank) {

                    val = 1;
                    
                    for(int j = 0; j < N + 1; j++) {
                        if(j < busy_side) {
                            val *= matrices[j][r * mas_razm[j] + mas[j]];
                        } else if(j > busy_side) {
                            val *= matrices[j][r * mas_razm[j - 1] + mas[j - 1]];
                        } else {
                            continue;
                        }
                    }

                    S[count] = val;
                }

                if(r < razm[busy_side]) {

                    tensor_val = r * busy_multipliers;

                    for(int i = 0; i < N; i++) {
                        tensor_val += mas[i] * multipliers[i];
                    }

                    res[count] = tensor[tensor_val];
                }

                count++;
                mas[N - 1]++;
            }

            for(int i = N - 1; i >= 1; i--) {
                if(mas[i] >= mas_razm[i]) {
                    mas[i - 1]++;
                    mas[i] = 0;
                }
            }
        }
    }
}

int
main(void)
{
    std::srand(std::time(nullptr)); 
    int N, rank;
    long long length = 1;
    std::cout << "Введите количество размерностей тензора: ";
    std::cin >> N;
    int* razm = new int[N];
    for(int i = 0; i < N; i++) {
        std::cout << "Введите " << i + 1 << "-ю" << " размерность: ";
        std::cin >> razm[i];
        length *= razm[i];
    }
    std::cout << "Введите ранг: ";
    std::cin >> rank;
    double* tensor;

    double** matrices = new double*[N];

    for(int i = 0; i < N; i++) {
        create_matrix(matrices[i], razm[i], rank);
    }

    // tensor_make_not_random(&tensor, razm, length, N);
    tensor_make(&tensor, length);
    
    for(int i = 0; i < N; i++) {
        delete[] matrices[i];
        create_matrix(matrices[i], razm[i], rank);
    }

    std::cout << std::endl << std::endl;

    int szfull = length;
    int ione = 1;
    int info = 654;
    int lwork = 64 * szfull;
    double * work = new double[lwork];
    char cN = 'N';
    int leftsize;
    int leftsize_trunc;

    double right_side_norm;
    double relative_residual = 1.0, end_relative_residual = 2.0;
    double diffrent = end_relative_residual - relative_residual;

    right_side_norm = dnrm2_(&szfull, tensor, &ione);

    std::cout << "Right side norm: " << right_side_norm << std::endl;
    double* S;
    double* T;
    int iteration = 0;



    while(diffrent  > 1.0e-8 && iteration < 100000)
    {
        end_relative_residual = relative_residual;
        iteration++;

        for(int i = 0; i < N; i++) {
            S = new double[length / razm[i] * rank];
            T = new double[length];
            create_S_T(matrices, i, rank, N - 1, razm, S, T, tensor);
            
            leftsize = length / razm[i];
            dgels_(&cN, &leftsize, &rank, &(razm[i]), S, &leftsize, T, &leftsize, work, &lwork, &info);

            relative_residual = 0.0;
            leftsize_trunc = leftsize - rank;

            for(int k = 0; k < razm[i]; k++)
            {
                double tmp = dnrm2_(&leftsize_trunc, T + k * leftsize + rank, &ione);

                relative_residual += tmp * tmp;
            }
            relative_residual = sqrt(relative_residual) / right_side_norm;

            for(int r = 0; r < rank; r++) {
                for(int j = 0; j < razm[i]; j++) {
                    matrices[i][r * razm[i] + j] = T[j * leftsize + r];
                }
            }

            delete[] S;
            delete[] T;

        }
        diffrent = end_relative_residual - relative_residual;
        if(iteration % 1000 == 0) {
            std::cout << diffrent << std::endl;
        }
    }
    std::cout << "Residual norm: " << relative_residual << std::endl;
    return 0;
}

