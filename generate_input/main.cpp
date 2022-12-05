#include <cstdio>
#include <iostream>
#include <algorithm>
#include <random>

using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

void generate_number(int digits, char* string) {
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 10.0); //range is [0, 10)

    int random = 0;
    for (int i = 0; i < digits; i++) {
        random = abs((int)dist(mt) % 10);
        while (i == 0 && random == 0) { // make sure the first bit is non-zero
            random = abs((int)dist(mt) % 10);
        }
        while (i == digits - 1 && random % 2 == 0) { // make sure the number is not even
            random = abs((int)dist(mt) % 10);
        }
        *(string + i) = random + '0';
    }
    *(string + digits) = '\0'; // add the end mark of this string
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
    	cerr << "Not enough inputs" << endl;
    }
    int min = atoi(argv[1]);
    int max = atoi(argv[2]);
    for (int i = min; i <= max; ++i) {
        char* n_string = (char*)malloc((i + 1) * sizeof(char));
        generate_number(i, n_string);
	cout << n_string << endl;
    }
}
