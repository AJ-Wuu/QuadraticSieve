#include "Quadratic-Sieve.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <vector>

using namespace std;

const uint64_t SIEVE_CHUNK = 65536;  // 2^16 -- the size of the chunk we sieve at a time

// modular exponentiation using the right-to-left binary method
uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t r = 1;
    while (b > 0) {
        if (b & 1) {
            r = r * a % m;
        }
        b >>= 1;
        a = a * a % m;
    }
    return r;
}

int64_t legendre_symbol(uint64_t a, uint64_t p) {
    unsigned long t = pow_mod(a, (p - 1) / 2, p);
    return t > 1 ? -1 : t;
}

// solve the congruence x^2 = n (mod p)
// reference: https://rosettacode.org/wiki/Tonelli-Shanks_algorithm
void tonelli_shanks(uint64_t n, uint64_t p, uint64_t *result) {
    if (p == 2) {
        result[0] = n;
        result[1] = n;
        return;
    }

    uint64_t S = 0;
    uint64_t Q = p - 1;
    while (Q % 2 == 0) {
        Q /= 2;
        ++S;
    }

    uint64_t z = 2;
    while (legendre_symbol(z, p) != -1) {
        ++z;
    }

    uint64_t c = pow_mod(z, Q, p);
    uint64_t R = pow_mod(n, (Q + 1) / 2, p);
    uint64_t t = pow_mod(n, Q, p);
    uint64_t M = S;

    while (t % p != 1) {
        int64_t i = 1;
        while (pow_mod(t, pow(2, i), p) != 1) {
            ++i;
        }
        uint64_t b = pow_mod(c, pow(2, M - i - 1), p);
        R = R * b % p;
        t = t * b * b % p;
        c = b * b % p;
        M = i;
    }

    result[0] = R;
    result[1] = p - R;
}

// get the ith bit in row
inline int64_t get_bit(uint64_t i, uint64_t *row) {
    return (row[i / sizeof(uint64_t)] & (1 << (i % sizeof(uint64_t)))) != 0;
}

// set the ith bit in row to 1
inline void set_bit(uint64_t i, uint64_t *row) {
    row[i / sizeof(uint64_t)] |= (1 << (i % sizeof(uint64_t)));
}

// set the ith bit in row to 0
inline void unset_bit(uint64_t i, uint64_t *row) {
    row[i / sizeof(uint64_t)] &= ~(1 << (i % sizeof(uint64_t)));
}

// toggle the ith bit in row
inline void toggle_bit(uint64_t i, uint64_t *row) {
    row[i / sizeof(uint64_t)] ^= (1 << (i % sizeof(uint64_t)));
}

// reference: http://www.personal.psu.edu/rcv4/467qs1.pdf
mpz_class quadratic_sieve(mpz_class &N) {
    vector<uint64_t> factor_base;
    mpz_class sqrt_N = sqrt(N);

    // 1 - Initialization
    // 1.1 pick the smoothness bound
    double log_N = mpz_sizeinbase(N.get_mpz_t(), 2) * log(2);
    uint64_t B = (uint64_t)ceil(exp(0.56 * sqrt(log_N * log(log_N)))) + 300;

    // 1.2 generate the factor base using a sieve: set p0 = 2, find odd primes p2 < p3 < ... < pK = p <= B
    char *sieve = new char[B + 1];
    memset(sieve, 1, B + 1);
    for (unsigned long p = 2; p <= B; ++p) {
        if (!sieve[p]) {
            continue;
        }

        if (mpz_legendre(N.get_mpz_t(), mpz_class(p).get_mpz_t()) == 1) {
            factor_base.push_back(p);
        }

        for (unsigned long i = p; i <= B; i += p) {
            sieve[i] = 0;
        }
    }
    delete[] sieve;

    // 2 - Sieving
    vector<uint64_t> X;
    double *Y = new double[SIEVE_CHUNK];
    vector<vector<uint64_t>> smooth;
    int fails = 0;
    uint64_t min_x = 0;  // sieve boundary
    uint64_t max_x = SIEVE_CHUNK;

    // calculate sieve index -- where to start the sieve -- for each factor base number
    uint64_t **fb_indexes = new uint64_t *[2];
    fb_indexes[0] = new uint64_t[factor_base.size()];
    fb_indexes[1] = new uint64_t[factor_base.size()];
    for (uint64_t p = 0; p < factor_base.size(); ++p) {
        // 1.3 find the solutions to x^2 = n (mod p)
        uint64_t idxs[2];
        mpz_class temp = N % mpz_class(factor_base[p]);
        tonelli_shanks(temp.get_ui(), factor_base[p], idxs);

        temp = idxs[0] - sqrt_N;
        temp = ((temp % factor_base[p]) + factor_base[p]) % factor_base[p];
        fb_indexes[0][p] = temp.get_ui();

        temp = idxs[1] - sqrt_N;
        temp = ((temp % factor_base[p]) + factor_base[p]) % factor_base[p];
        fb_indexes[1][p] = temp.get_ui();
    }

    double last_estimate = 0;
    uint64_t next_estimate = 1;
    // sieve new chunks until having enough smooth numbers
    while (smooth.size() < (factor_base.size() + 20)) {
        // generate Y vector for the sieve, containing log approximations that fit in machine words
        for (uint64_t t = 1; t < SIEVE_CHUNK; ++t) {
            // calculating a log estimate is expensive, so don't do it for every Y[t]
            if (next_estimate <= (t + min_x)) {
                mpz_class y = (sqrt_N + t + min_x) * (sqrt_N + t + min_x) - N;

                // to estimate the 2 logarithm, just count the number of bits that v takes up
                last_estimate = mpz_sizeinbase(y.get_mpz_t(), 2);

                // the higher t gets, the less the logarithm of Y[t] changes
                next_estimate = next_estimate * 1.8 + 1;
            }

            Y[t] = last_estimate;
        }

        // perform the actual sieve
        for (uint64_t p = 0; p < factor_base.size(); ++p) {
            double lg = log(factor_base[p]) / log(2);

            for (uint64_t t = 0; t < 2; ++t) {
                while (fb_indexes[t][p] < max_x) {
                    Y[fb_indexes[t][p] - min_x] -= lg;
                    fb_indexes[t][p] += factor_base[p];
                }

                // p = 2 only has one modular root.
                if (factor_base[p] == 2) {
                    break;
                }
            }
        }

        // factor all values whose logarithms were reduced to approximately zero using trial division
        double threshold = log(factor_base.back()) / log(2);
        for (uint64_t i = 0; i < SIEVE_CHUNK; ++i) {
            if (fabs(Y[i]) < threshold) {
                mpz_class y = (sqrt_N + i + min_x) * (sqrt_N + i + min_x) - N;
                smooth.push_back(vector<uint64_t>());

                for (uint64_t p = 0; p < factor_base.size(); ++p) {
                    while (mpz_divisible_ui_p(y.get_mpz_t(), factor_base[p])) {
                        mpz_divexact_ui(y.get_mpz_t(), y.get_mpz_t(), factor_base[p]);
                        smooth.back().push_back(p);
                    }
                }

                if (y == 1) {
                    // This V was indeed B-smooth.
                    X.push_back(i + min_x);

                    // Break out of trial division loop if we've found enough smooth numbers.
                    if (smooth.size() >= (factor_base.size() + 20)) {
                        break;
                    }
                } else {
                    // This V was apparently not B-smooth, remove it.
                    smooth.pop_back();
                    ++fails;
                }
            }
        }

        min_x += SIEVE_CHUNK;
        max_x += SIEVE_CHUNK;
    }

    // 3 - Linear Algebra
    // 3.1 form the matrix
    uint64_t **matrix = new uint64_t *[factor_base.size()];
    int row_words = (smooth.size() + sizeof(uint64_t)) / sizeof(uint64_t);  // the amount of words needed to accomodate a row in the augmented matrix
    for (uint64_t i = 0; i < factor_base.size(); ++i) {
        matrix[i] = new uint64_t[row_words];
        memset(matrix[i], 0, row_words * sizeof(uint64_t));
    }
    for (uint64_t s = 0; s < smooth.size(); ++s) {
        // for each factor in the smooth number, add the factor to the corresponding element in the matrix
        for (uint64_t p = 0; p < smooth[s].size(); ++p) {
            toggle_bit(s, matrix[smooth[s][p]]);
        }
    }

    // 3.2 solve λM = 0 (mod 2) where λ = (λ1, λ2, ..., λ(K+2))
    // Gauss elimination, the dimension of the augmented matrix is factor_base.size() * (smooth.size() + 1)
    uint64_t i = 0, j = 0;
    while (i < factor_base.size() && j < (smooth.size() + 1)) {
        uint64_t maxi = i;

        // Find pivot element.
        for (uint64_t k = i + 1; k < factor_base.size(); ++k) {
            if (get_bit(j, matrix[k]) == 1) {
                maxi = k;
                break;
            }
        }
        if (get_bit(j, matrix[maxi]) == 1) {
            swap(matrix[i], matrix[maxi]);

            for (uint64_t u = i + 1; u < factor_base.size(); ++u) {
                if (get_bit(j, matrix[u]) == 1) {
                    for (int64_t w = 0; w < row_words; ++w)
                        matrix[u][w] ^= matrix[i][w];
                }
            }
            ++i;
        }
        ++j;
    }

    mpz_class a, b;
    uint64_t **back_matrix = new uint64_t *[factor_base.size()];  // copy of matrix that to perform back-substitution on
    for (uint64_t i = 0; i < factor_base.size(); ++i) {
        back_matrix[i] = new uint64_t[row_words];
    }

    // 4 - Factorization
    // 4.1 compute x = x1^λ1 * x2^λ2 * ... * x(K+2)^λ(K+2) mod n and y mod n
    uint64_t *x = new uint64_t[smooth.size()];
    uint64_t *combination = new uint64_t[factor_base.size()];
    // loop until found a non-trivial factor
    do {
        // copy the Gauss eliminated matrix
        for (uint64_t i = 0; i < factor_base.size(); ++i) {
            memcpy(back_matrix[i], matrix[i], row_words * sizeof(uint64_t));
        }

        // clear the x vector
        memset(x, 0, smooth.size() * sizeof(uint64_t));

        // perform back-substitution on our matrix that's now in row echelon form to get x
        int64_t i = factor_base.size() - 1;
        while (i >= 0) {
            // count non-zero elements in current row
            int64_t count = 0;
            int64_t current = -1;
            for (uint64_t c = 0; c < smooth.size(); ++c) {
                count += get_bit(c, back_matrix[i]);
                current = get_bit(c, back_matrix[i]) ? c : current;
            }

            // empty row -> advance to next
            if (count == 0) {
                --i;
                continue;
            }

            // the system is underdetermined and we can choose x[current] freely
            // to avoid the trivial solution -> avoid always setting it to 0
            uint64_t val = count > 1 ? rand() % 2 : get_bit(smooth.size(), back_matrix[i]);

            x[current] = val;
            for (int64_t u = 0; u <= i; ++u) {
                if (get_bit(current, back_matrix[u]) == 1) {
                    if (val == 1) {
                        toggle_bit(smooth.size(), back_matrix[u]);
                    }
                    unset_bit(current, back_matrix[u]);
                }
            }

            if (count == 1) {
                --i;
            }
        }

        a = 1, b = 1;
        memset(combination, 0, sizeof(uint64_t) * factor_base.size());  // combine the factor base to get square
        for (uint64_t i = 0; i < smooth.size(); ++i) {
            if (x[i] == 1) {
                for (uint64_t p = 0; p < smooth[i].size(); ++p) {
                    ++combination[smooth[i][p]];
                }
                b *= (X[i] + sqrt_N);
            }
        }

        for (uint64_t p = 0; p < factor_base.size(); ++p) {
            for (uint64_t i = 0; i < (combination[p] / 2); ++i) {
                a *= factor_base[p];
            }
        }
    } while (a % N == b % N || a % N == (-b) % N + N);  // if a = +/- b (mod N) -> found a trivial factor, run the loop again to find a new a and b
    b -= a;

    // 4.2 compute m = gcd(x - y, n)
    mpz_class factor;
    mpz_gcd(factor.get_mpz_t(), b.get_mpz_t(), N.get_mpz_t());
    for (uint64_t i = 0; i < factor_base.size(); ++i) {
        delete[] matrix[i];
        delete[] back_matrix[i];
    }

    delete[] combination;
    delete[] Y;
    delete[] fb_indexes[0];
    delete[] fb_indexes[1];
    delete[] fb_indexes;
    delete[] matrix;
    delete[] back_matrix;
    delete[] x;

    return factor;
}
