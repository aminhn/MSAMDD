MSAMDD Version 1.0 - Exact Multiple Sequence Alignment by Synchronized Decision Diagrams

Copyright (C) 2019 Amin Hosseininasab
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

DESCRIPTION
    MSAMDD is an multiple sequence alignment algorithm designed to generate optimal sequence alignments.
    	
REFERENCE
    Please cite the below paper when using MSAMDD. Technical and algorithm specifics may also be found in this paper.
    Hosseininasab, A. and van Hoeve, W.J., 2019. 
    Exact Multiple Sequence Alignment by Synchronized Decision Diagrams. 
    INFORMS Journal on Computing (to appear).

CONTACT
    amin.hosseininasab@warrington.ufl.edu

INSTALLATION
    Designed for Linux systems. Installation from Makefile. Simply execute make from the terminal at the MSAMDD directory. The output is two executables: msa_aff and msa_cnv. 

EXTERNAL SOFTWARE REQUIREMENTS FOR INSTALLATION
    MSAMDD requires the installation of IBM ILOG CPLEX Optimization Studio to solve its mathematical programming models. 
    CPLEX may be downloaded from https://www.ibm.com/products/ilog-cplex-optimization-studio
    CPLEX should be installed in the opt/ibm/ILOG directory at the root folder of the Linux system. 
    MSAMDD automatically detects and uses the latest version at this directory. In case CPLEX is installed in a different directory, 
    the user must change the default directory in the Makefile at line 5.
    MSAMDD further requires the MUSCLE program for its heuristic alignments. Although MSAMDD can be executed without MUSCLE, its solution time will be significantly slowed.
    MUSCLE may be downloaded from https://www.drive5.com/muscle/
    All credit for the usage of the MUSCLE software goes to Robert C. Edgar, and should be properly cited as 
    Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput Nucleic Acids Res. 32(5):1792-1797 
    The MSAMDD package includes MUSCLE for convenience, but claims no right to its ownership.

DOWNLOAD
    MSAMDD can be downloaded from https://github.com/aminhn/MSAMDD

PENALTY FUNCTIONS
    MSAMDD supports linear, concave, and convex penalty functions for optimal sequence alignments through two programs: msa_aff and msa_cnv. 
    msa_aff is intended for constant, linear, or affine gap penalty functions. The affine penalty function is of the form a + bx. Here, x is the gap length and a and b are 
    constants specified by the user. In sequence alignment terminology, constant a corresponds to the opening gap penalty, and constant b corresponds to the extension gap penalty. 
    Constant and linear penalty functions are a special case of the affine gap penalty. A constant penalty function is achieved when b=0, and a linear gap penalty function is achieved when a=0. 
    msa_cnv is intended for convex gap penalty functions. The penalty functions are of the form a + bx + cSQRT(x) or a + bx + cLn(x), where c is a constant specified by the user
    SQRT(x) is the square root of the gap length, and Ln(x) is the natural logarithm of the gap length. In sequence alignment terminology, constant c corresponds to the convex gap penalty.
    Note that the affine gap penalty is a special case of the convex gap penalty with c=0. In that regard, although msa_cnv may be used to accommodate affine, linear or constant gap penalty functions
    the runtime of msa_cnv is significantly slower than msa_aff. Therefore, the user is advised to use msa_cnv only when convex gap penalty functions are required, and use msa_aff otherwise.

COMMAND LINE COMMANDS
	
	-in		path of input .fa file containing the sequences to be aligned. The input file MUST be of the fasta .fa format
	-out		path of output file containing the aligned sequences. This file may be of .fa or .txt format
	-time		time limit (in seconds) for the algorithm to find the optimal solution. MSAMDD terminates after the time limit and outputs the best found solution (not necessarily optimal).
	-start 		path of file containing a warm-start alignment of sequences. A high quality warm-start solution may help reduce the solution time of MSAMDD. This file may be of .fa or .txt format
	-submat		path containing the substitution matrix for alignment rewards (default is ./Data/blosum.ncbi)	
	-op 		opening gap penalty (double precision, default is 12 in msa_aff, and 8 in msa_cnv)
	-ep 		extension gap penalty (double precision, default is 2.22 in msa_aff, and 2 in msa_cnv)
	-cp 		convex gap penalty (double precision, default is 2 in msa_cnv) (available only for msa_cnv)
	-cfun 		convex gap penalty function. 1: SQRT, or 2: Ln (default is 1 and SQRT) (available only for msa_cnv)

COMMAND LINE USAGE EXAMPLE
    The below command uses msa_cnv for an alignment of sequences in file BB11001.fa found at directory ./Data/bali3/ 
    The algorithm uses a convex gap penalty function of 8 + 2x + 2SQRT(x), and substitution matrix of blosum.ncbi found at directory ./Data/
    The algorithm terminates at most after 3600 seconds with the best found solution, which is outputted to file BB11001_opt.txt in directory ./Data/results/

    ./msa_cnv -in ./Data/bali4/BB11001.fa -out ./Data/results/BB11001_opt.txt -time 3600 -op 8 -ep 2 -cp 2 -cfun 1 -submat ./Data/blosum.ncbi
    

DATABASE
    The MSAMDD package includes sequences for the BAliBASE versions 1 (./Data/bali1/) and BAliBASE version 4 (./Data/bali4). These databases may be downloaded from http://www.lbgi.fr/balibase/
    The MSAMDD package includes the BAliBASE databases for convenience, but claims no right to their ownership.


