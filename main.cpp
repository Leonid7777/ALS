#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <random>
#include <vector>

extern "C" {
double dznrm2_(const int *, const std::complex<double> *, const int *);
void zgels_(const char *, const int *, const int *, const int *,
            std::complex<double> *, const int *, std::complex<double> *,
            const int *, std::complex<double> *, const int *, int *);
}

void tensor_make_not_random(std::complex<double> **tensor, int *razm,
                            int length, int N, int R, double &right_side_norm,
                            double noise, std::complex<double> **end_tensor) {
    int ione = 1;
    int szfull = length;

    *tensor = new std::complex<double>[length];
    *end_tensor = new std::complex<double>[length];
    for (int i = 0; i < length; i++) {
        std::complex<double> val(0, 0);
        (*tensor)[i] = val;
        (*end_tensor)[i] = val;
    }

    double phi;
    for (int r = 0; r < R; r++) {
        phi = -100.0 + 200.0 * double(std::rand()) / RAND_MAX;
        for (int i = 0; i < length; i++) {
            std::complex<double> val(cos(phi * i), sin(phi * i));
            (*tensor)[i] += val;
            (*end_tensor)[i] += val;
        }
    }

    right_side_norm = dznrm2_(&szfull, (*tensor), &ione);

    std::random_device rd;
    std::mt19937 gen(rd());
    double a = noise * right_side_norm / sqrt(length);
    std::normal_distribution<double> distribution(0.0, a / sqrt(2));

    std::complex<double> sample;

    for (int i = 0; i < length; i++) {
        sample.real(distribution(gen));
        sample.imag(distribution(gen));
        (*tensor)[i] += sample;
    }
}

void create_matrix(std::complex<double> *&matrix, long long razm, int rank) {
    matrix = new std::complex<double>[razm * rank];
    for (int r = 0; r < rank; r++) {
        for (int row = 0; row < razm; row++) {
            std::complex<double> val(-1.0 + 2.0 * double(std::rand()) / RAND_MAX,
                                    -1.0 + 2.0 * double(std::rand()) / RAND_MAX);
            matrix[row + r * razm] = val;
        }
    }
}

void create_S_T(std::complex<double> **matrices, int busy_side, int rank, int N,
                int *razm, std::complex<double> *S, std::complex<double> *res,
                std::complex<double> *tensor) {

    int *multipliers = new int[N + 1];

    multipliers[N] = 1;
    for (int i = N - 1; i >= 0; i--) {
        multipliers[i] = multipliers[i + 1] * razm[i + 1];
    }

    long long busy_multipliers = multipliers[busy_side];
    for (int i = busy_side; i < N; i++) {
        multipliers[i] = multipliers[i + 1];
    }

    int *mas = new int[N];
    int *mas_razm = new int[N];

    for (int i = 0; i < N + 1; i++) {
        if (i < busy_side) {
            mas_razm[i] = razm[i];
        } else if (i > busy_side) {
            mas_razm[i - 1] = razm[i];
        } else {
            continue;
        }
    }

    long long count = 0, tensor_val = 0;
    std::complex<double> val;

    long long max_r = std::max(rank, razm[busy_side]);

    for (long long r = 0; r < max_r; r++) {

        for (int pas = 0; pas < N; pas++) {
            mas[pas] = 0;
        }

        while (mas[0] < mas_razm[0]) {
            while (mas[N - 1] < mas_razm[N - 1]) {

                if (r < rank) {

                    val = 1;

                    for (int j = 0; j < N + 1; j++) {
                        if (j < busy_side) {
                            val *= matrices[j][r * mas_razm[j] + mas[j]];
                        } else if (j > busy_side) {
                            val *= matrices[j][r * mas_razm[j - 1] + mas[j - 1]];
                        } else {
                            continue;
                        }
                    }

                    S[count] = val;
                }

                if (r < razm[busy_side]) {

                    tensor_val = r * busy_multipliers;

                    for (int i = 0; i < N; i++) {
                        tensor_val += mas[i] * multipliers[i];
                    }

                    res[count] = tensor[tensor_val];
                }

                count++;
                mas[N - 1]++;
            }

            for (int i = N - 1; i >= 1; i--) {
                if (mas[i] >= mas_razm[i]) {
                    mas[i - 1]++;
                    mas[i] = 0;
                }
            }
        }
    }

    delete[] multipliers;
    delete[] mas;
    delete[] mas_razm;
}



//tensor - I J K
//S - JK R 
//T - JK I

double ALS(std::complex<double> *tensor, std::complex<double> **matrices, int N,
           int rank, int *razm, int length, double &right_side_norm, std::map<std::vector<int>, std::complex<double>> rand_ten) {
    int szfull = length;
    int ione = 1;
    int info = 654;
    int lwork = 64 * szfull;
    std::complex<double> *work = new std::complex<double>[lwork];
    char cN = 'N';
    int leftsize, leftsize_trunc;
    double relative_residual = 100.0, end_relative_residual = 2.0;
    double diffrent = abs(end_relative_residual - relative_residual);

    std::cout << "Right side norm: " << right_side_norm << std::endl;
    std::complex<double> *S;
    std::complex<double> *T = new std::complex<double>[length];
    int iteration = 0;

    int* multipliers = new int[N];

    while (diffrent > 1.0e-9 && iteration < 100000) {
        end_relative_residual = relative_residual;
        iteration++;

        for (int i = 0; i < N; i++) {

            multipliers[N - 1] = 1;

            for(int j = N - 2; j >= 0; j--) {
                multipliers[j] = 1;
                if(j != i) {
                    for(int k = j + 1; k < N; k++) {
                        if(k != i) {
                            multipliers[j] *= razm[k]; 
                        }
                    }
                }                
            }

            S = new std::complex<double>[length / razm[i] * rank];

            create_S_T(matrices, i, rank, N - 1, razm, S, T, tensor);

            leftsize = length / razm[i];

            std::complex<double> *S_l = new std::complex<double>[leftsize * rank];
            std::complex<double> *T_l = new std::complex<double>[length];

            relative_residual = 0.0;

            //перебор по столбцам тензора и перебора матрриц
            for(int j = 0; j < razm[i]; j++) {

                //заполнение массивов по стотбцу получившегося тензора и строкам произведения матриц
                int count_elem_i = 0;
                for(auto elem = rand_ten.begin(); elem != rand_ten.end(); elem++) {
                    if((elem->first)[i] == j) {
                        T_l[count_elem_i] = elem->second;
                        count_elem_i++;
                    }
                }

                int pas = 0;
                int count = 0;
                for(auto elem = rand_ten.begin(); elem != rand_ten.end(); elem++) {
                    if((elem->first)[i] == j) {
                        pas = 0;
                        for(int p = 0; p < N; p++) {
                            if(p != i) {
                                pas += (elem->first)[p] * multipliers[p];
                            }
                        }
                        for(int t = 0; t < rank; t++) {
                            S_l[count + t * count_elem_i] = S[pas + t * leftsize];
                        }
                        count++;
                    }
                }

                int rows_S_l = count_elem_i * rank;

                zgels_(&cN, &count_elem_i, &rank, &(ione), S_l, &count_elem_i, T_l, &count_elem_i,
                    work, &lwork, &info);


                leftsize_trunc = count_elem_i - rank;

                double tmp = dznrm2_(&leftsize_trunc, T_l + rank, &ione);

                relative_residual += tmp * tmp;

                for (int r = 0; r < rank; r++) {
                    matrices[i][r * razm[i] + j] = T_l[r];
                }

            }

            relative_residual = sqrt(relative_residual) / right_side_norm;

            delete[] S;
            delete[] S_l;
            delete[] T_l;
        }

        diffrent = end_relative_residual - relative_residual;
        if (iteration % 1000 == 0) {
            std::cout << diffrent << std::endl;
        }
    }

    std::cout << "iteration: " << iteration << std::endl;

    delete[] work;
    delete[] T;
    delete[] multipliers;

    return relative_residual;
}

void create_ALS_tensor(std::complex<double> **matrices, int rank, int N,
                       int *razm, std::complex<double> *tensor, int length,
                       double right_side_norm) {
    int ione = 1;
    int szfull = length;
    int *mas = new int[N];
    long long count = 0, tensor_val = 0;
    std::complex<double> val;
    long long max_r = rank;
    std::complex<double> sum;
    std::complex<double> *res = new std::complex<double>[length];

    for (int i = 0; i < N; i++) {
        mas[i] = 0;
    }

    while (mas[0] < razm[0]) {
        while (mas[N - 1] < razm[N - 1]) {

            res[count] = -tensor[count];

            for (int i = 0; i < rank; i++) {
                sum.real(1);
                sum.imag(0);
                for (int j = 0; j < N; j++) {
                    sum *= matrices[j][mas[j] + i * razm[j]];
                }
                res[count] += sum;
            }

            count++;
            mas[N - 1]++;
        }

        for (int i = N - 1; i >= 1; i--) {
            if (mas[i] >= razm[i]) {
                mas[i - 1]++;
                mas[i] = 0;
            }
        }
    }

    std::cout << "Real error: " << dznrm2_(&szfull, res, &ione) / right_side_norm
                << std::endl;

    delete[] res;
    delete[] mas;
}

void create_random_tensor(std::complex<double> *tensor,
                          std::map<std::vector<int>, std::complex<double>> *rand_ten,
                          int length, int len_rand_ten, int* razm, int N) {

    int* multipliers = new int[N];
    multipliers[N - 1] = 1;
    for(int i = N - 2; i >= 0; i--) {
        multipliers[i] = razm[i + 1] * multipliers[i + 1];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distrib(0, length - 1);
    int count = 0;

    while (count < len_rand_ten) {
        std::vector<int> val;
        int num = distrib(gen);
        int pas = num;

        for(int i = 0; i < N; i++) {
            val.push_back(num / multipliers[i]);
            num %= multipliers[i];
        }

        if ((*rand_ten).find(val) == (*rand_ten).end()) {
            (*rand_ten)[val] = tensor[pas];
            count++;
        }
    }
    delete[] multipliers;
}

int main(void) {
    std::srand(std::time(nullptr));
    int N, rank, R;
    int length = 1;
    double right_side_norm, noise, pro_inc_it;
    std::cout << "Введите количество размерностей тензора: ";
    std::cin >> N;

    int *razm = new int[N];
    for (int i = 0; i < N; i++) {
        std::cout << "Введите " << i + 1 << "-ю"
                << " размерность: ";
        std::cin >> razm[i];
        length *= razm[i];
    }

    std::cout << "Введите желаемый ранг тензора (без учёта шума): ";
    std::cin >> R;

    std::cout << "Введите долю шума: ";
    std::cin >> noise;

    std::cout << "Введите ранг: ";
    std::cin >> rank;

    std::cout << "Введите долю полученных элементов: ";
    std::cin >> pro_inc_it;

    std::complex<double> *tensor;
    std::complex<double> **matrices = new std::complex<double> *[N];
    std::map<std::vector<int>, std::complex<double>> rand_ten;

    int len_rand_ten = length * pro_inc_it;

    std::complex<double> *end_tensor;


    for (int i = 0; i < N; i++) {
        create_matrix(matrices[i], razm[i], rank);
    }

    tensor_make_not_random(&tensor, razm, length, N, R, right_side_norm, noise,
                            &end_tensor);

    create_random_tensor(tensor, &rand_ten, length, len_rand_ten, razm, N);

    std::cout << std::endl << std::endl;

    double relative_residual =
        ALS(tensor, matrices, N, rank, razm, length, right_side_norm, rand_ten);

    std::cout << "Residual norm: " << relative_residual << std::endl;

    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += razm[i];
    }

    sum = sqrt((rank * sum) / length);

    std::cout << "Expected error: " << noise * sum << std::endl;

    create_ALS_tensor(matrices, rank, N, razm, end_tensor, length,
                        right_side_norm);

    // for(auto i = rand_ten.begin(); i != rand_ten.end(); i++) {
    //     std::cout << i->first << " " << i->second << std::endl;
    // }

    delete[] razm;
    delete[] tensor;
    delete[] end_tensor;
    for (int i = 0; i < N; i++) {
        delete[] matrices[i];
    }
    delete[] matrices;

    return 0;
}
