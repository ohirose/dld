
# Dependent Landmark Drift

This is an implementation of a point matching algorithm, dependent landmark drift (DLD). 
Details of the algorithm are available [here](https://arxiv.org/abs/1711.06588).

## Demo

Demo movies are available [here](https://youtu.be/3sDNes4n_RY). If you are a MATLAB (Windows) user,
demo codes can be executed by following the instructions included in `Body`, `Hand`, and
`Face` folders.

## Preparation: shape model construction

Statistical shape models used as shape prior information MUST be constructed before execution.
The construction of a statistical shape model might be technical especially if source datasets
are images or point clouds without point-by-point correspondence. For these cases,
see [one](http://www.sciencedirect.com/science/article/pii/S1361841509000425) of review papers,
for example. If source datasets are set of points of correspondence, the construction of the
statistical shape model is simple. For example, the PCA-based statistical shape model is
constructed as follows.

1. Pre-align shapes for removing shape variations originated from the rigid transformation of objects.
2. Convert a point set corresponding to an object shape into a single vector.
3. Compute the mean shape as the sample average of shape vectors.
4. Compute the covariance matrix of shape vectors.
5. Output the mean shape and eigenvectors (+ eigenvalues if needed).
   Shape variations of the object geometry are represented as eigenvectors.

Please see `trainingHand.m` in the `Hand` folder, for example, to construct a statistical shape model
from point sets with point-by-point correspondence.

## Compilation

### MacOS and Linux

1. Install the LAPACK library if not installed.
2. Download the zip file that includes source codes and uncompress it.
3. Move into the root directory of the uncompressed folder using the terminal window.
4. Type `make` in the terminal window.

OpenMP can be used if your MacOS is quite old, e.g., 10.11. Install `gcc-mp-6` and
replace `make` with `make CFLAGS=-DUSE_OPENMP CC=gcc-mp-6`.

### Windows

Not required. Use the binary file in this repository.

## Usage 

For MacOS and Linux, type the following command in the terminal window.

` ./dld <target: X> <model: shape variations> <model: mean shape> (+options) `

For Windows, type the following command in the DOS prompt.

` dld <target: X> <model: shape variations> <model: mean shape> (+options) `

Some instructions are printed by typing `./dld` in the terminal window, or `dld` in the DOS prompt.

### Terms and symbols

- Target point set: the point set corresponding to the reference shape.
- N: The number of the target point set.
- M: The number of points in the mean shape.
- K: The number of shape variations.
- D: Dimension of the space in which the mean shape and the target point set are embedded.

### Input data (required)

- 1st argument: The target shape represented as a matrix of size N x D.
- 2nd argument: The matrix of size MD x K corresponding to shape variations.
- 3rd argument: The mean shape represented as a vector of size MD.

Only tab-separated files are accepted. A column of the matrix required for the 3rd argument
MUST be the eigenvector multiplied by the square root of the corresponding eigenvalue.
The multiplication of the eigenvalue are required for imposing the Tikhonov regularizer
which plays a role in the smooth displacement field.

### Parameters

- `-w [real]`: Noise probability in (0,1).
- `-g [real]`: Regularization constant. Positive. A value in (0,1) is recommended.

### Options

- `-i`: Pre-alignment by estimating the best rigid transformation.
- `-z`: Estimation of shape deformation and translation. Scale and rotation are not estimated.
- `-l`: Toggle layout. Read the 2nd and 3rd arguments in the order of xyzxyz instead of xxyyzz.
- `-q`: Quiet mode. Print nothing.
- `-s`: Save the optimization path.
- `-v`: Print the version of this software.

### Configurations

- `-c [real]`: Convergence tolerance. 
- `-n [int ]`: The maximum number of EM loops.
- `-o [name]`: File name of the output shape.
- `-p [name]`: File name of the parameter file. If turned on, parameters specified by command-line options are ignored.
- `-f [real]`: IFGT error bound, positive. A value between 1e-2 and 1e-4 is recommended.
- `-y [int ]`: The number of points to be sampled for Nystrom method, positive.

Nystrom method and IFGT are options for speeding up the execution. Nystrom method is usually
faster than the direct computation if M and N are moderately large and the number of points
to be sampled is set to much smaller than both N and M. If the number of sampled points in
the Nystrom method is not enough, the optimization will become unstable. IFGT is recommended
to be turned on if the number of cores in a CPU is greater than or equal to four and if N and M
are more than at least several thousands. <!-- Otherwise, execution will be slower if IFGT is turned on. -->
Nystrom method and IFGT can't be turned on at the same time.

