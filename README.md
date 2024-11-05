## Numerical Methods for Partial Differential Equations: laboratory material

### Required software

The laboratories rely on the `deal.II` finite element library ([official website](https://dealii.org), [GitHub repository](https://github.com/dealii/dealii)), which is already available in mk modules (`module load gcc-glibc dealii`).

For the visualization of numerical solutions, we will use Paraview ([official website](https://www.paraview.org/)). Please download and install it in advance from their [download page](https://www.paraview.org/download/). Notice that Paraview should be installed on the host system, not inside the Docker container.

### Laboratory folder structure
Each laboratory has its own folder `lab-XY`, containing:
- a sub-folder `lab-XY/text`, with the text of the exercises;
- a sub-folder `lab-XY/src`, with the starting source code;
- a file `CMakeLists.txt`, used by CMake for the configuration and compilation;
- a sub-folder `lab-XY/solution` (added after the lecture), with the solution to the exercises (both theory and code).

### Compiling
To build the executable for a laboratory, make sure you have loaded the needed modules with
```bash
$ module load gcc-glibc dealii
```
Then, from the folder `lab-XY`, run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```
The executable for the laboratory will be created into `lab-XY/build`, and can be executed through
```bash
$ ./lab-XY
```

### Examples folder

The repository contains a folder `examples`, that will contain short examples of code. Examples will be added during the course: refer to the `README.md` file in each example subfolder for information on each example and building instructions.