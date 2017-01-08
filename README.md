#smdlt-iterative
This code requires [OpenBLAS](https://github.com/xianyi/OpenBLAS). [OpenMP](https://computing.llnl.gov/tutorials/openMP/) is optional but recommanded.
This code can run iterative method of semidilute system in **planar extensional flow**
I only test it with gcc, may not be compatible with icc or clang. The code is *only* compativble with C99. 

#Usage
I use Makefile to compile the code, Some common input is also in the Makefile.
##Makefile
* `CC=gcc-5` change this to the compiler you are using I am using gcc/5.3 here.
* `OMPGLAG=-fopenmp` the OpenMP flag, for icc user use `-qopenmp`. Clang does not support OpenMP therefore leave it blank
* `HPATH=` and `LIBPARH=` theses are the library and header path for OpenBLAS. See it's [wiki](https://computing.llnl.gov/tutorials/openMP/) for detail.

* Other Inputs: there are some inputs in the Makefile, the meaning of them are described in it. 
* `TAEGETS=` the target program name, make sure to change the name before its rule when changing the target name.

To compile the code just run the following in terminal:
````
make clean && make
````

The inputs that are not included in the Makefile need to be specified through arguments
##Input arguments
* Running mode use flog `-m`, 0 - FD, 1 - HI. Default 0. 
* Chain length use flag `-n`. Default value 20
* Box size use flag `-l`. Default value 20. Number of chains will be calculated based on the box size and overlap comcentration specifed in Makefile
* Normalized concentration use flag `-c`. Default value 1 (overlap comcentration)
* Flow rate (hencky strain rate) use flag `-r`. *must* be specified. Default value 0.
* Number of OpenMP threads use flag `-t`. Default 1 (no parallelism)

Example:
````
./program -t 4 -n 20 -l 40.0 -c 0.1 -r 0.01 -m 0
````
#Others
* To automatically start the next iteration after the first iteration finished I included one method using a `run.sh` file. 

* This is an old code and I did not use it for a long time. It may contain typo.