cmake_minimum_required(VERSION 3.20)
project(epmc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(HDF5 REQUIRED COMPONENTS CXX)

add_executable(epmc
    src/main.cpp
    src/Args.cpp
    src/Mesh.cpp
    src/Geometry.cpp
    src/Phonon.cpp
    src/Simulation.cpp
)

target_include_directories(epmc PRIVATE include ${HDF5_INCLUDE_DIRS})
target_link_libraries(epmc PRIVATE ${HDF5_CXX_LIBRARIES})

option(EPMC_ENABLE_OPENMP "Enable OpenMP parallel loops" ON)
if(EPMC_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(epmc PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(epmc PRIVATE EPMC_USE_OPENMP=1)
        message(STATUS "OpenMP enabled for epmc")
    elseif(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # Fallback for AppleClang + Homebrew libomp where CMake may fail to auto-detect OpenMP.
        set(EPMC_LIBOMP_ROOT "/opt/homebrew/opt/libomp")
        if(EXISTS "${EPMC_LIBOMP_ROOT}/include/omp.h" AND EXISTS "${EPMC_LIBOMP_ROOT}/lib/libomp.dylib")
            target_compile_options(epmc PRIVATE -Xpreprocessor -fopenmp -I${EPMC_LIBOMP_ROOT}/include)
            target_link_options(epmc PRIVATE -L${EPMC_LIBOMP_ROOT}/lib -Wl,-rpath,${EPMC_LIBOMP_ROOT}/lib)
            target_link_libraries(epmc PRIVATE omp)
            target_compile_definitions(epmc PRIVATE EPMC_USE_OPENMP=1)
            message(STATUS "OpenMP enabled via Homebrew libomp: ${EPMC_LIBOMP_ROOT}")
        else()
            message(WARNING "EPMC_ENABLE_OPENMP=ON but OpenMP was not found (auto and Homebrew fallback failed). Building without OpenMP.")
        endif()
    else()
        message(WARNING "EPMC_ENABLE_OPENMP=ON but OpenMP was not found. Building without OpenMP.")
    endif()
endif()
