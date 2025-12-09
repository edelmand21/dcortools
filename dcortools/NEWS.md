### dcortools 0.1.7

9.12.2025

Dominic Edelmann

- Removed Specification to C++11
- Solved a bug for distance covariance version for survival with type.X = distance (thanks to Na You, Essex!!)

### dcortools 0.1.6

9.12.2022

Dominic Edelmann

- Added #define NDEBUG 1 to RcppExports.cpp after #include <RcppArmadillo.h>
- Checked that no assert was left in any(!) files 

### dcortools 0.1.5

7.12.2022

Dominic Edelmann

- Added #define NDEBUG 1 to matrix_vector_operations.cpp after #include <RcppArmadillo.h>
-Checked that no assert was left in compiled code.

## dcortools 0.1.4

1.12.2022

Dominic Edelmann

- Updated DESCRIPTION file (adding a Description field)
- Added "NDEBUG 1" to header of RcppExports.cpp to avoids problems with "assert_rtn" 


### dcortools 0.1.3

30.11.2022

Dominic Edelmann

- Native routines registered in src/init.c


### dcortools 0.1.2

9.7.2021

Dominic Edelmann

- Removed "rlc" from "Imports"
- Fixed a bug in "distcov.test" when metr.X = "discrete", metr.Y = "discrete",
  method = "permutation"



### dcortools 0.1.1

7.7.2021

Dominic Edelmann

Deleted all functions containing "<<-"


### dcortools 0.1.0

8.3.2021

Dominic Edelmann

Initial submission of dcortools to CRAN.

