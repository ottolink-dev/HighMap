# HighMap

A C++ library to generate two-dimensional terrain heightmaps for software rendering or video games. HighMap is the core solver of the [Hesiod](https://github.com/otto-link/Hesiod) GUI.

<p align="center"><img src="https://github.com/otto-link/HighMap/assets/121820229/4451f45b-c74e-4a06-9431-9302477a12c5" width="512"></p>

>[!WARNING] 
> Work in progress, use at your own risk!

## License

This project is licensed under the GNU General Public License v3.0.

## Python binding

[pyHighMap](https://github.com/otto-link/pyHighMap) allows you to use a few functionalities (for now) of HighMap directly in your Python projects.

```bash
pip3 install pyhighmap
```

## Getting started

### Getting the sources

Use `git` to retrieve the sources: 
``` bash
git clone git@github.com:otto-link/HighMap.git
cd HighMap
git submodule update --init --recursive
```

### Building

#### Linux

Build by making a build directory (i.e. `build/`), run `cmake` in that dir, and then use `make` to build the desired target.

Example:
``` bash
mkdir build && cd build
cmake ..
make
```

Simple examples are provided for most of the functions, for instance:
```
bin/./ex_noise_fbm
eog ex_noise_fbm_*.png
```

#### Windows

Use `PowerShell` to clone the repository (because Visual Studio won't clone the submodules) using `git` command lines (see above).

Install the missing OpenSource dependencies using `vcpkg`:
```
vcpkg install libpng glm opencl assimp opencv[openexr]
```

You should then be able to build the sources using Visual Studio.

#### macOS

Install the dependencies using [Homebrew](https://brew.sh):
``` bash
brew install gsl assimp glm opencv libomp spdlog opencl-headers opencl-clhpp-headers
```

Then build as on Linux:
``` bash
mkdir build && cd build
cmake ..
make
```

>[!NOTE]
> Apple Clang does not bundle OpenMP, the `libomp` package provides it. If
> `libomp` is missing the library still builds, but the few OpenMP-based
> kernels (flow accumulation, drainage basins, local metrics) run
> single-threaded. Apple Silicon builds are tuned with `-mcpu=apple-m1`.


### CMake Integration

To integrate HighMap into your CMake-based project, follow these steps:

1. Add the library to your project's `CMakeLists.txt`:
   ```cmake
   add_subdirectory(HighMap)
   target_link_libraries(your_project_target highmap)
   ```
2. Link the `highmap` library to your target, as shown above.
   
### Usage examples

See the `examples` folder.

## Third Party Code

This project includes code from external open-source projects.
See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).