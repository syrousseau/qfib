<p align="center"><img src="sources/Documentation/html/qfib.png" width="200" /></p>

# Overview
This software compress and decompress dMRI brain tractograms. Those tractograms should have been generated using a method that ensure a constant stepsize along each individual fiber. This property can be checked using our software. We currently can compress .tck files, and decompressed .qfib files. This is the C++ reference implementation for the paper: 

>**QFib: Fast and Efficient Brain Tractogram Compression** <br/>
*Corentin Mercier\*, Sylvain Rousseau\*, Pietro Gori, Isabelle Bloch and Tamy Boubekeur.*<br/>
(\* C. Mercier · S. Rousseau contributed equally to this work)<br/>
NeuroInformatics 2020<br/>
DOI: 10.1007/s12021-020-09452-0<br/>

This is not the exact same version of the source code that is used to measure performance for the corresponding paper. Performance might have been slightly affected during refactoring.

Copyright(C) 2019-2020 Corentin Mercier, Sylvain Rousseau, Pietro Gori, Isabelle Bloch and Tamy Boubekeur

All right reserved

## Release Notes ##
### v1.0 ###
Initial version

# Building and Running

This program uses CMake. It has been tested on Linux using Ubuntu 18.04 (LTS) with gcc 7.3 and Ubuntu 16.04 (LTS) with gcc 5.4. Building on Windows 10 or MacOS has not been tested. Direct compatibility with visual c++ is not ensured as their version of OpenMP (that we used to parallelize our code) is too old, and we are using getopt.h. The easiest way to get a working version *under* windows is to use WSL. To build the source code, you can use the following commands: 

```
mkdir qfib/build
cd qfib/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```

Then, you can simply compress a file with the default parameters (8 bits spherical Fibonacci, incore, with data validity verification) using the following command line:

```
./qfib path/file.tck
```

and decompress a file in the same way
```
./qfib path/file.qfib
```
It is not necessary to indicate the parameters of the compression during decompression as these informations are stored in the compressed file. Some option can be activated to improve the performances (remove checkings), modify the quantization method or precision, or to display more informations. 
Those option can be activated using the following flags: 

```
-u or --out_of_core: force out-of-core, in case there is not enough RAM (default is in core)
-o or --octahedral: force octahedral quantification, only for compression (default is fibonacci)
-b or --16bits: compress using 16 bits (default is 8 bits)
-t or --timings: display compression and decompression times
-e or --error: display the error when compressing fibers
-d or --detailed-error: display a detailed error when compressing fibers
-v or --verbose: display additional details
-n or --no-verif: no verification of data constant stepsize for better performance
-f or --fa + filename: compare fa values after compression and decompression
-c or --compression-ratio: display compression ratio when compressing fibers
-h or --help: access to this help
```

# Docker

To use docker with the code, you can use the following commands: 
```
sudo docker build --rm -f C++Base -t cpp-build-base:1.0 .
sudo docker build --rm -f QfibDocker -t qfib:1.0 .
sudo docker run -v /path/to/folder/containing/the/fibers:/src/data qfib:1.0 ../data/file -o -t
```

# Authors

* [**Corentin Mercier**](https://perso.telecom-paristech.fr/comercier/)
* [**Sylvain Rousseau**](https://perso.telecom-paristech.fr/srousseau/) 
* [**Isabelle Bloch**](https://perso.telecom-paristech.fr/bloch/)
* [**Pietro Gori**](https://perso.telecom-paristech.fr/pgori/)
* [**Tamy Boubekeur**](https://perso.telecom-paristech.fr/boubek)

# License

This project is licensed under the MIT license - see the [LICENSE](LICENSE) file for details.

All right reserved. The Authors
# Data 
Data can be generated using [MrTrix](https://www.mrtrix.org).

Some examples are available here: https://gitlab.com/comercier/qfib-data
