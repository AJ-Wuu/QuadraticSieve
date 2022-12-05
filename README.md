# Quadratic Sieve
## Environmental Settings (in Linux System)
* GMP is needed to run the files (see [details](https://gmplib.org/))
### (OPTIONAL) Updates -- if you are on your personal computer (i.e., you can act as admin), you may do so
```bash
$ sudo apt-get update
$ sudo apt-get install g++
$ sudo apt-get install m4
$ sudo apt-get install libgmp-dev
```
### Get things ready
#### To keep the compatibility of all Makefiles, we install everything under the current directory `QuadraticSieve`
#### I'm working on a lab machine, so you may want to change the following absolute directories to yours
```bash
$ cd QuadraticSieve                                                # make sure you are under this code directory
$ wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz
$ tar -xJvf gmp-6.2.1.tar.xz
$ cd gmp-6.2.1
$ ./configure --prefix=/srv/home/awu53/QuadraticSieve/usr/local    # note: 1. prefix must be an absolute path
```

## Code Usage
```bash
$ cd QuadraticSieve                                                # skip this line if you're already inside the dir
$ ls                                                               # showing four dirs: cpp, generate_input, gmp-6.2.1, openmp
$ 
$ cd generate_input
$ make
$ ./generate_input 18 40 >> test.in                                # syntax: ./generate_input min-length max-length >> output-file
$                                                                  # note that the min and max length are both INCLUSIVE
$                                                                  # this code doesn't support over 40 digits
$ 
$ cd ../cpp
$ make
$ ./quadratic_sieve_cpp < ../generate_input/test.in
$ 
$ cd ../openmp
$ make
$ ./quadratic_sieve_openmp < ../generate_input/test.in
```
