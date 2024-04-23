cmake_minimum_required(VERSION 3.0)
set(CMAKE_TOOLCHAIN_FILE "C:/DEVELOP/vcpkg/scripts/buildsystems/vcpkg.cmake")
project(AlphaShapeExample)

# Find CGAL package
find_package(CGAL CONFIG REQUIRED)

# Set source files
set(SOURCE_FILES ./src/ConnectGrid.cpp)

# Add executable
add_executable(AlphaShapeExample ${SOURCE_FILES})

# Link CGAL
target_link_libraries(AlphaShapeExample PRIVATE CGAL::CGAL)
