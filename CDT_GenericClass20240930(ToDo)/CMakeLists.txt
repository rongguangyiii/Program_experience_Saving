
cmake_minimum_required(VERSION 3.20)
set(CMAKE_TOOLCHAIN_FILE "D:/dev_1st/vcpkg/scripts/buildsystems/vcpkg.cmake")
project (T1 VERSION 3.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 添加VCPKG的CMake工具链
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
        CACHE STRING "")
endif()


# 指定构建的子目录
add_subdirectory(${PROJECT_SOURCE_DIR}/src)

#生成构建目标
add_executable(${PROJECT_NAME} "${PROJECT_SOURCE_DIR}/src/demo.cpp" ${SRC_FILE})

#程序输出路径
# set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/out/${CMAKE_BUILD_TYPE})

#可执行文件后缀
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX ${PROJECT_VERSION})
set_target_properties(${PROJECT_NAME} PROPERTIES RELEASE_POSTFIX ${PROJECT_VERSION}) 

#指定头文件搜索路径
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
# 查找必要的包
#Eigen库
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

#VTK库
find_package(CGAL CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)

