cmake_minimum_required(VERSION 3.5)
project(bsplines VERSION 3.1.0)

# Find python and Boost - both are required dependencies
find_package(PythonLibs 3.7 REQUIRED)
find_package(Boost COMPONENTS python REQUIRED)


# Without this, any build libraries automatically have names "lib{x}.so"
set(CMAKE_SHARED_MODULE_PREFIX "")

set(bsplines_VERSION_MAJOR 3)
set(bsplines_VERSION_MINOR 0)
set(CMAKE_CXX_STANDARD 11)
include_directories(/home/ashish/intel/oneapi/intelpython/python3.7/include/python3.7m)
add_executable(bspline_basis examples/bspline_basis.cpp)
