cmake_minimum_required(VERSION 3.20)
project(PhonoMC LANGUAGES CXX)

option(BUILD_TESTING "Build the optional local regression suite" OFF)
include(CTest)

# C++ language requirements
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(HDF5 REQUIRED COMPONENTS CXX)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

# File/schema parsing has no solver, SimulationConfig or OpenMP dependency.
add_library(phonomc_material_data STATIC src/material/MaterialData.cpp src/material/ScatteringMatrix.cpp src/material/CallawayOperator.cpp)
target_include_directories(phonomc_material_data PUBLIC include)
add_library(phonomc_material_io STATIC src/material/MaterialDataReader.cpp src/material/ScatteringMatrixReader.cpp)
target_include_directories(phonomc_material_io PUBLIC include PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(phonomc_material_io PUBLIC phonomc_material_data ${HDF5_CXX_LIBRARIES})

# In-memory physics model: no HDF5, filesystem, configuration or loading policy.
add_library(phonomc_material STATIC src/PhononMaterial.cpp src/material/CallawayMaterial.cpp)
target_include_directories(phonomc_material PUBLIC include)
target_link_libraries(phonomc_material PUBLIC phonomc_material_data)

# Compile the core once; all consumers inherit the same dependency configuration.
add_library(phonomc_core STATIC
    src/SimulationConfig.cpp
    src/SurfaceMesh.cpp
    src/SimulationDomain.cpp
    src/material/MaterialFactory.cpp
    src/MonteCarloSolver.cpp
    src/solver/ParticleStorage.cpp
    src/solver/CellParticleIndex.cpp
    src/solver/RtaCollisionOperator.cpp
    src/solver/CarrierCollisionOperator.cpp
    src/solver/NonlinearCallawayCollision.cpp
    src/solver/FullMatrixCollisionOperator.cpp
    src/solver/FullMatrixSource.cpp
    src/solver/TemperatureGradientDrive.cpp
    src/solver/Initialization.cpp
    src/solver/Transport.cpp
    src/solver/Interfaces.cpp
    src/solver/InterfaceModeSampler.cpp
    src/solver/AcousticInterfaceSampler.cpp
    src/solver/ChannelCell.cpp
    src/solver/ReservoirSampler.cpp
    src/solver/Reservoirs.cpp
    src/solver/RoughBoundaries.cpp
    src/solver/RoughBoundarySampler.cpp
    src/solver/HeatSource.cpp
    src/solver/PrescribedHeatSource.cpp
    src/solver/PhononSourceOperator.cpp
    src/solver/Diagnostics.cpp
    src/solver/Observables.cpp
    src/solver/GridStatistics.cpp
    src/solver/ModalResampler.cpp
)

# File formatting/I/O has no physics, HDF5 or OpenMP dependency.
add_library(phonomc_output STATIC src/solver/ResultWriter.cpp)
target_include_directories(phonomc_output PUBLIC include)
# Accounting is independently testable without HDF5 or OpenMP linkage.
add_library(phonomc_accounting STATIC src/solver/BoundaryLedger.cpp src/solver/EnergyLedger.cpp)
target_include_directories(phonomc_accounting PUBLIC include)
target_link_libraries(phonomc_core PUBLIC phonomc_output phonomc_accounting)

# Public include paths
target_include_directories(phonomc_core PUBLIC include ${HDF5_INCLUDE_DIRS})

# Core library dependencies
target_link_libraries(phonomc_core PUBLIC phonomc_material phonomc_material_io Threads::Threads)

# Optional shared-memory parallelism
option(PHONOMC_ENABLE_OPENMP "Enable OpenMP parallel loops" ON)
if(PHONOMC_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(phonomc_core PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(phonomc_core PUBLIC PHONOMC_USE_OPENMP=1)
        target_link_libraries(phonomc_material PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(phonomc_material PRIVATE PHONOMC_USE_OPENMP=1)
        message(STATUS "OpenMP enabled for PhonoMC")
    else()
        message(WARNING "OpenMP was not found. Building serial PhonoMC.")
    endif()
endif()

add_executable(PhonoMC src/main.cpp)
target_link_libraries(PhonoMC PRIVATE phonomc_core)

# The public source checkout does not require the local validation directory.
if(BUILD_TESTING)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
        message(FATAL_ERROR
            "BUILD_TESTING requires the local tests/ suite. "
            "Configure this source checkout with -DBUILD_TESTING=OFF.")
    endif()
    include("${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
endif()
