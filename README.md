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
* NTL Library is also required (see [official page](https://libntl.org/))
