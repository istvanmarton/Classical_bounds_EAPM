# Classical_bounds_EAPM
This repository contains the implementation of calculating the $L_d$ norm of an $W_{xyz}$ witness tensor standing for the entanglement assisted prepare-and-measure scenario. The code was used to calculate the numerical results in the article ["Bound entanglement-assisted prepare-and-measure scenarios"](https://doi.org/10.1088/2058-9565/ae095f) by I. Márton, E. Bene, and T. Vértesi. The code is written in C++ using OpenMP.

## Usage:
The code can be compiled with the

    g++ L_OMP_org.c -fopenmp -o L_OMP_org -lm

command. The input of the program should be az $W_{xyz}$ tensor of size ($n \times n \times n$), where the tensor is given in a file as a single row or column having $n \times n \times n$ entries, where the given entry stands for the $(x + y \times n + z \times n \times n)$-th element of the tensor. The program can be invoked with the

    ./L_EAPM number_of_threads filename_of_tensor d

where
 + 'number_of_threads' is the number of threads the code utilizes,
 + 'filename_of_tensor' is the name of the file containing the input tensor,
 + 'd' is an integer ($d \ge 2$) indicating the order of the L norm one intends to calculate.

The code calculates the L norm of the tensor considering only the strategies corresponding to the Stirling numbers of the second kind and writes the result to the screen. Apart from that, the code writes the strategies belonging to the L norm to a file 'strategy_filename_of_tensor'. For example if the name containing the input tensor has the name 'something.txt', the output file name containing the strategies will have 'strategy_something.txt'.
