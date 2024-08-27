
cmake_minimum_required(VERSION 3.10)
project(ParabolaPlotter)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/out/${CMAKE_BUILD_TYPE})

# 添加VCPKG的CMake工具链
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "")
endif()

# 查找必要的包
find_package(CGAL CONFIG REQUIRED)
# 查找必要的包
find_package(Eigen3 REQUIRED)
#VTK库
find_package(VTK CONFIG REQUIRED)

# 包含头文件
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include) 

# 添加源文件
set(SOURCE_DIR src)
set(SOURCE_FILES
    ${SOURCE_DIR}/demo.cpp 
    ${SOURCE_DIR}/cdt.cpp 
    ${SOURCE_DIR}/Parabola.cpp 
    ${SOURCE_DIR}/uniformGrid.cpp 
    ${SOURCE_DIR}/kdTreee.cpp 
   )

# Create an executable with the specified sources
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# 链接必要的库
# Link the CGAL target to the executable
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)
target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES})

