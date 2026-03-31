cmake_minimum_required(VERSION 3.10)
project(lbimproved
    VERSION 1.0.0
    DESCRIPTION "Fast Dynamic Time Warping nearest-neighbor retrieval with LB_Improved"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library target
add_library(lbimproved INTERFACE)
target_include_directories(lbimproved INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Options
option(LBIMPROVED_BUILD_TESTS "Build unit tests" ON)
option(LBIMPROVED_BUILD_BENCHMARKS "Build benchmarks" ON)
option(LBIMPROVED_BUILD_EXAMPLES "Build examples" ON)

if(LBIMPROVED_BUILD_TESTS)
    enable_testing()
    add_executable(unittesting tests/unittesting.cpp)
    target_link_libraries(unittesting PRIVATE lbimproved)
    add_test(NAME unittesting COMMAND unittesting)
endif()

if(LBIMPROVED_BUILD_BENCHMARKS)
    add_executable(benchmark benchmarks/benchmark.cpp)
    target_link_libraries(benchmark PRIVATE lbimproved)
    target_compile_options(benchmark PRIVATE -O2)
endif()

if(LBIMPROVED_BUILD_EXAMPLES)
    add_executable(example examples/example.cpp)
    target_link_libraries(example PRIVATE lbimproved)
    target_compile_options(example PRIVATE -O2)
endif()
