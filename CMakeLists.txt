cmake_minimum_required(VERSION 3.20)
project(PhonoMC LANGUAGES CXX)

include(CTest)

# 基础编译设置
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 查找 HDF5 库
find_package(HDF5 REQUIRED COMPONENTS CXX)

# 添加可执行文件及其源码
add_executable(PhonoMC
    src/main.cpp
    src/SimulationConfig.cpp
    src/SurfaceMesh.cpp
    src/SimulationDomain.cpp
    src/PhononMaterial.cpp
    src/MonteCarloSolver.cpp
)

# 设置包含目录
target_include_directories(PhonoMC PRIVATE include ${HDF5_INCLUDE_DIRS})

# 链接 HDF5 库
# Linux 下推荐使用 HDF5::HDF5_CPP_LIBRARIES 或 ${HDF5_LIBRARIES}
target_link_libraries(PhonoMC PRIVATE ${HDF5_CXX_LIBRARIES})

# OpenMP 配置
option(PHONOMC_ENABLE_OPENMP "Enable OpenMP parallel loops" ON)
if(PHONOMC_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(PhonoMC PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(PhonoMC PRIVATE PHONOMC_USE_OPENMP=1)
        message(STATUS "OpenMP enabled for PhonoMC")
    else()
        message(WARNING "OpenMP was not found. Building without OpenMP.")
    endif()
endif()

if(BUILD_TESTING)
    add_executable(PhonoMCConfigTests
        tests/SimulationConfigTests.cpp
        src/SimulationConfig.cpp
    )
    target_include_directories(PhonoMCConfigTests PRIVATE include)
    target_compile_definitions(PhonoMCConfigTests PRIVATE
        PHONOMC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(NAME simulation_config COMMAND PhonoMCConfigTests)

    add_executable(PhonoMCDomainTests
        tests/SimulationDomainTests.cpp
        src/SimulationConfig.cpp
        src/SurfaceMesh.cpp
        src/SimulationDomain.cpp
    )
    target_include_directories(PhonoMCDomainTests PRIVATE include)
    target_compile_definitions(PhonoMCDomainTests PRIVATE
        PHONOMC_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
    )
    add_test(NAME simulation_domain COMMAND PhonoMCDomainTests)

    add_executable(PhonoMCMaterialTests
        tests/PhononMaterialTests.cpp
        src/SimulationConfig.cpp
        src/PhononMaterial.cpp
    )
    target_include_directories(PhonoMCMaterialTests PRIVATE include ${HDF5_INCLUDE_DIRS})
    target_link_libraries(PhonoMCMaterialTests PRIVATE ${HDF5_CXX_LIBRARIES})
    add_test(NAME phonon_material COMMAND PhonoMCMaterialTests)
endif()
