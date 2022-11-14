#include "Quadratic-Sieve.h"
#include <gmpxx.h>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
    mpz_class N(argv[1]);  // read from command line

    // Quadratic Sieve is the fastest for integers up to 100 digits
    // hence, this implementation only takes care of the input less than or equal to 100 digits
    if (mpz_sizeinbase(N.get_mpz_t(), 10) > 100) {
        cerr << N << " is over 100 digits\n" << endl;
        return 0;
    }

    stack<mpz_class> factors;
    vector<mpz_class> result;
    factors.push(N);
    result.push_back(N);

    while (!factors.empty()) {
        mpz_class factor = factors.top();
        factors.pop();

        // if the factor is prime
        if (mpz_probab_prime_p(factor.get_mpz_t(), 10)) {
            result.push_back(factor);
            continue;
        }
        else {
            // perfect power case
            // Quadratic Sieve cannot deal with perfect power very well, so we split it out
            if (mpz_perfect_power_p(factor.get_mpz_t())) {
                mpz_class root, rem;

                // check root remainders up to half of the amount of bits in factor
                uint32_t max = mpz_sizeinbase(factor.get_mpz_t(), 2) / 2;
                for (uint32_t n = 2; n < max; ++n) {
                    mpz_rootrem(root.get_mpz_t(), rem.get_mpz_t(), factor.get_mpz_t(), n);
                    if (rem == 0) {
                        // push the n root factors
                        for (uint32_t i = 0; i < n; ++i) {
                            factors.push(root);
                        }
                        break;
                    }
                }
            }
            else {
                mpz_class f = quadratic_sieve(factor);
                factors.push(f);
                factors.push(factor / f);
           }
        }
    }

    cout << result.at(0) << " = ";
    for (unsigned int i = 1; i < result.size() - 1; i++) {
        cout << result.at(i) << " * ";
    }
    cout << result.at(result.size() - 1) << endl;
}


