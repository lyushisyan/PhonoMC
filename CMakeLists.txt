cmake_minimum_required(VERSION 3.20)
project(epmc LANGUAGES CXX)

# 基础编译设置
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 查找 HDF5 库
find_package(HDF5 REQUIRED COMPONENTS CXX)

# 添加可执行文件及其源码
add_executable(epmc
    src/main.cpp
    src/Args.cpp
    src/Mesh.cpp
    src/Geometry.cpp
    src/Phonon.cpp
    src/Simulation.cpp
)

# 设置包含目录
target_include_directories(epmc PRIVATE include ${HDF5_INCLUDE_DIRS})

# 链接 HDF5 库
# Linux 下推荐使用 HDF5::HDF5_CPP_LIBRARIES 或 ${HDF5_LIBRARIES}
target_link_libraries(epmc PRIVATE ${HDF5_CXX_LIBRARIES})

# OpenMP 配置
option(EPMC_ENABLE_OPENMP "Enable OpenMP parallel loops" ON)
if(EPMC_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(epmc PRIVATE OpenMP::OpenMP_CXX)
        target_compile_definitions(epmc PRIVATE EPMC_USE_OPENMP=1)
        message(STATUS "OpenMP enabled for epmc")
    else()
        message(WARNING "OpenMP was not found. Building without OpenMP.")
    endif()
endif()