# PolyDG Library

**_Politecnico di Milano (ITALY)_**

**Author :** Andrea Vescovini

**Mailto :** andrea.vescovini@mail.polimi.it

## The project

#### PolyDG
PolyDG is a library for the approximation of elliptic problems with Discontinuous
Galerkin finite element methods on polyhedral grids.  
It makes use of the technique of *expression templates* in order to speed up the
assembling of the matrix and of the right-hand-side of the linear system, moreover
this approach allows a flexible implementation of different models starting from
the variational formulation of the problem.  
It is written in C++11.  
Developed with gcc 5.3.0, tested with clang???

#### convergence_analisys
The folder `convergence_analisys` contains examples of use of the library concerning
the calculation of the convergence order of the method, both with respect to the
refinment of the mesh and with respect to the degree of the polynomials used for
the basis of the finite element space.

## How to install

#### Requirments
* The library **Eigen** is used for linear algebra, it does not have to be compiled
  since it is a header only library.  
  The project has been developed with the **version 3.3.3**  
  See http://eigen.tuxfamily.org/index.php?title=Main_Page for more information.

* The library **GetPot** is used to parse the tests commaind line arguments, it is not
  necessary for building the library but it is needed if you want to run the tests
  or the examples provided about the convergence analisys. It is only one header
  file, so you do not need to compile anything.  
  See http://getpot.sourceforge.net/ for more information.

* The tool **Doxygen** is used to generate the documentation.  
  See http://www.stack.nl/~dimitri/doxygen/index.html for more information.

#### Installation of the library
First of all you have to edit the `Makefile.inc` in this folder and insert in
the first three variables the directories where Eigen (`EIGEN_INC`) and GetPot (`GETPOT_INC`) are in your system
and the directory where you want to install the library (`POLYDG_PATH`).

Then enter the folder `libPolyDG` and type:
```shell
make all RELEASE=yes
```
This builds both a static and a dynamic version of the library, builds the executables
of the tests and the documentation. Then typing:
```shell
make install
```
you can install the library in the directory specified with `POLYDG_PATH`.

The option `RELEASE=yes` is recommended in order to enable all the optimizations
and get the full performance from the code. Without it you compile in the debug mode
and you enable some output messages during the execution.

Typing `make help` you can get some information about other kinds of commands.

If the library has been successfully built, you should find in the folder `libPolyDG/lib`
the two libraries, static and dynamic.  
In the folder `libPolyDG/doc` there should be the documentation, in HTML and LaTeX format. In the main page of the documentation you can find a tutorial for the use of the library.  
In the folder `libPolyDG/bin` there should be
the executables of the tests, that run with the configuration file `libPolyDG/data.pot`.
If you have build the dynamic version of the library before building the tests (for
example if you have used `make all`), then the executables have been linked with the
dynamic one, so remember to tell the loader where the library is editing
the environmental variable `LD_LIBRARY_PATH`.

If you want to uninstall the library removing all the files that had been copied
in `POLYDG_PATH`, type:
```shell
make uninstall
```

#### Examples
Before building the examples, be sure to have installed the library.  
Enter the directory `convergence_analisys` and type:
```shell
make all RELEASE=yes
```
in order to compile the two examples. Again remember to specify `RELEASE=yes` for the optimization.

Then you can run the exacutables and you can modify some parameters from the
input files `h_convergence.pot` and `p_convergence.pot`.
