# Quadratic Sieve
## Environmental Settings (in Linux System)
* GMP is required to run the files (see [details](https://gmplib.org/)). To install GMP, please do:
```
$ sudo apt-get update
$ sudo apt-get install g++
$ sudo apt-get install m4
$ sudo apt-get install libgmp-dev
$ wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz
$ tar -xJvf gmp-6.2.1.tar.xz
$ cd gmp-6.2.1
$ make
$ ./configure
```
* Boost Library is also required (see [documents](https://www.boost.org/doc/libs/1_39_0/more/getting_started/unix-variants.html) and [downloads](https://www.boost.org/users/download/)). To install it, please do:
```
$ wget https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz
$ tar -xvzf boost_1_80_0.tar.gz
$ ./bootstrap.sh    # note that I directly write this library under QuadraticSieve directory
                    # if you do it somewhere else, please change cpp/Makefile accordingly
$ ./b2 install
```
