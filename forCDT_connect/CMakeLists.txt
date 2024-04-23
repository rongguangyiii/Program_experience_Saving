cmake_minimum_required(VERSION 3.0)
set(CMAKE_TOOLCHAIN_FILE "D:/dev_VCPKG/vcpkg/scripts/buildsystems/vcpkg.cmake")
project(newCDT02)

# Find CGAL package
find_package(CGAL CONFIG REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include) 

# Organize source files in a dedicated directory variable
set(SOURCE_DIR src)
set(SOURCE_FILES
    ${SOURCE_DIR}/demo.cpp 
    ${SOURCE_DIR}/filter_Tri.cpp 
    
)

# Create an executable with the specified sources
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Link the CGAL target to the executable
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)
