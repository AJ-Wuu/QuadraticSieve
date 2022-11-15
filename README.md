# Quadratic Sieve
## Environmental Settings (in Linux System)
* GMP is needed to run the files (see [details](https://gmplib.org/)), and NTL Library is also needed (see [official page](https://libntl.org/))
### (OPTIONAL) if you are on your personal computer (i.e., you can act as admin), you may do the updates
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
                                                                   #       2. create .../usr/local first if it yells on this step
$ make
$ make check
$ make install
# install NTL now
$ cd ..                                                            # get back to QuadraticSieve
$ wget https://libntl.org/ntl-11.5.1.tar.gz
$ tar xf ntl-11.5.1.tar.gz
$ cd ntl-11.5.1/src
$ ./configure PREFIX=/srv/home/awu53/QuadraticSieve/usr/local GMP_PREFIX=/srv/home/awu53/QuadraticSieve/usr/local
$ make
$ make check
$ make install
```
### (OPTIONAL) Test the environment
1. Under `QuadraticSieve`, create [foo.cpp](https://libntl.org/doc/tour-ex2.html)
```cpp
#include <NTL/ZZ.h>
#include <NTL/vector.h>

using namespace std;
using namespace NTL;

int main() {
   Vec<ZZ> v;
   cin >> v;
   long n = v.length();
   v.SetLength(2*n);
   for (long i = 0 ; i < n; i++)
      v[n+i] = v[n-1-i];
   cout << v << "\n";
}
```
2. Execute in command line: ```g++ -I ./usr/local/include -L ./usr/local/lib -L ./usr/local/lib  foo.cpp -o foo -lntl -lgmp -lm```
   * This is of the form ```g++ -I<prefix>/include -L<prefix>/lib -L<gmp_prefix>/lib  foo.c -lntl -lgmp -lm```
   * Note: the paths are all relative -- indicating that we are still under `QuadraticSieve`
3. If it builds successfully, then try ```./foo [1 -2 3]``` and you shall get ```[1 -2 3 3 -2 1]```
