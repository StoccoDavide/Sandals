# Sandals

`Sandals` is a library for the numerical integration of Ordinary Differential Equations (ODEs) and Differential Algebraic Equations (DAEs) written in C++17.

Are you looking for the online documentation? Visit [this link](https://stoccodavide.github.io/Sandals/)!

## Installation

### Quick and dirty

`Sandals` is a header-only library that depends only on [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) (version >= 5.0.1), so the quick and dirty way of installing it is by simply copying the `include` directory to your project and make sure to have [`Eigen`](https://eigen.tuxfamily.org/index.php?title=Main_Page) available however you see fit. Alternatively, you can do things properly and use `CMake` (version >= 3.14).

You also need to install [`Optimist`](https://stoccodavide.github.io/Optimist/) to use root-finding algorithms. But don't worry, it is already fetched by `CMake` if you use it.

Optionally, you can also install CERN's [`ROOT`](https://root.cern.ch) (version >= 6.34.08) to plot the results of your simulations.

### CMake

If you are using CMake, you can add the library as a subdirectory in your project.

```cmake
add_subdirectory(path/to/Sandals)
target_link_libraries(your_target PRIVATE Sandals::Sandals)
```

If you already have `Sandals` somewhere on your system, you can use `find_package` directly.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  Sandals
  ${YOUR_DESIRED_SANDALS_VERSION}
  NO_MODULE
)

target_link_libraries(your_target PRIVATE Sandals::Sandals)
```

## Authors

- Davide Stocco <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: davide.stocco@unitn.it

- Enrico Bertolazzi <br>
  University of Trento <br>
  Department of Industrial Engineering <br>
  email: enrico.bertolazzi@unitn.it

Aka...

```
▗▄▄▄  ▄   ▄  ▐▌    ▗▞▀▜▌▄▄▄▄     ▐▌    ▗▄▄▖ ▗▞▀▚▖ ▄▄▄ ▄   ▄
▐▌  █ █   █  ▐▌    ▝▚▄▟▌█   █    ▐▌    ▐▌ ▐▌▐▛▀▀▘█    █   █
▐▌  █  ▀▄▀▗▞▀▜▌         █   █ ▗▞▀▜▌    ▐▛▀▚▖▝▚▄▄▖█     ▀▀▀█
▐▙▄▄▀     ▝▚▄▟▌               ▝▚▄▟▌    ▐▙▄▞▘          ▄   █
                                                       ▀▀▀
```

## License

The `Sandals` project is distributed under the BSD 2-Clause License - see the [LICENSE](https://StoccoDavide.github.io/Sandals/LICENSE) file for details.

Here's what the license entails:

1. Anyone can copy, modify and distribute this software.
2. You have to include the license and copyright notice with each and every distribution.
3. You can use this software privately.
4. You can use this software for commercial purposes.
5. This software is provided without warranty.
6. The software author or license can not be held liable for any damages inflicted by the software.

## Acknowledgements

The name `Sandals` has been suggested by [Mattia Piazza](https://github.com/Mattiapzz) as a poor imitation of the Lawrence Livermore National Laboratory's (LLNL) [`Sundials`](https://computing.llnl.gov/projects/sundials) (SUite of Nonlinear and DIfferential/ALgebraic Equation Solvers).
