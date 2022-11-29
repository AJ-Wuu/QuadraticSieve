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
# install GMP first -- this order matters
$ cd QuadraticSieve                                                # make sure you are under this code directory
$ wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz
$ tar -xJvf gmp-6.2.1.tar.xz
$ cd gmp-6.2.1
$ ./configure --prefix=/srv/home/awu53/QuadraticSieve/usr/local    # note: 1. prefix must be an absolute path
```
