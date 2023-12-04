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
        (*tensor)[i] = 1.0 / (i + 1);
    }
}

void
create_matrix(double* &matrix, long long razm, int rank)
{
    std::srand(std::time(nullptr));
    matrix = new double[razm * rank];
    for(int r = 0; r < rank; r++) {
        for(int row = 0; row < razm; row++) {
            matrix[row + r * razm] = double(std::rand()) / RAND_MAX;
        }
    }
}

double*
create_S(long long length, double** matrices, int busy_side, int rank, int N, int* razm)
{
    double* S = new double[length * rank];
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

    long long count = 0;

    for(long long r = 0; r < rank; r++) {

        for(int pas = 0; pas < N; pas++) {
            mas[pas] = 0;
        }

        while(mas[0] < mas_razm[0]) {
            while(mas[N - 1] < mas_razm[N - 1]) {

                double val = 1;

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
    return S;
}

double*
create_T(double* tensor, long long length, int* razm, int busy_side, int N)
{
    
    double *res = new double[length];


    long long* multipliers = new long long[N + 1];

    multipliers[0] = 1;
    for(int i = 1; i < N + 1; i++) {
        multipliers[i] = multipliers[i - 1] * razm[i - 1];
    }

    long long busy_multipliers = multipliers[busy_side];
    for(int i = busy_side; i < N; i++) {
        multipliers[i] = multipliers[i + 1];
    }

    double* mas = new double[N];
    double* mas_razm = new double[N];

    for(int i = 0; i < N + 1; i++) {
        if(i < busy_side) {
            mas_razm[i] = razm[i];
        } else if(i > busy_side) {
            mas_razm[i - 1] = razm[i];
        } else {
            continue;
        }
    }

    for(int i = 0; i < N; i++) {
        std::cout << mas_razm[i] << " ";
    }
    std::cout << std::endl;

    long long count = 0, val;

    for(int r = 0; r < razm[busy_side]; r++) {

        for(int pas = 0; pas < N; pas++) {
            mas[pas] = 0;
        }

        while(mas[0] < mas_razm[0]) {
            while(mas[N - 1] < mas_razm[N - 1]) {
                val = r * busy_multipliers;

                for(int i = 0; i < N; i++) {
                    val += mas[i] * multipliers[i];
                }

                res[count] = tensor[val];
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

    return res;
}

int
main(void)
{
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
    tensor_make(&tensor, length);

    double** matrices = new double*[N];

    for(int i = 0; i < N; i++) {
        create_matrix(matrices[i], razm[i], rank);
    }

    int szfull = length;
    int ione = 1;
    int info = 654;
    int lwork = 64 * szfull;
    double * work = new double[lwork];
    char cN = 'N';
    int leftsize;
    int leftsize_trunc;

    double right_side_norm;
    double relative_residual = 1.0;

    right_side_norm = dnrm2_(&szfull, tensor, &ione);

    std::cout << "Right side norm: " << right_side_norm << std::endl;
    double* S;
    double* T;
    int iteration = 0;



    while(relative_residual > 1.0e-3 && iteration < 10000)
    {
        iteration++;

        for(int i = 0; i < N; i++) {
            S = create_S(length / razm[i], matrices, i, rank, N - 1, razm);
            T = create_T(tensor, length, razm, i, N - 1);
            
            leftsize = length / razm[i];
            dgels_(&cN, &leftsize, &rank, &(razm[i]), S, &leftsize, T, &leftsize, work, &lwork, &info);
            std::cout << "info: " << info << std::endl;

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

            std::cout << "Residual norm: " << relative_residual << std::endl;
        }
    }


    return 0;
}
