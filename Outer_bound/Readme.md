> ### This is a program that uses VCPKG package management tool to install CGAL library and test the output of 2D point cloud outer edge point set.
> ## Author： Liu Guangying
> ## Email:   liugy36@mail2.sysu.edu.cn 

# 概览
- 该项目主要是测试使用VCPKG经典模式安装CGAL库，查找2D点云模型的外边缘点集，注意不是凸包点集。使用cgal库中的delunay算法和alpha-shape算法。
- 测试点云模型放在了./model/points.txt中。
- 文件夹./src下有两个测试代码，都是可以用的，使用时只需要在CMakeLists.txt文件中修改重新build即可：
```cpp
set(SOURCE_FILES ./src/***.cpp)
```

20240407
这两天因为需要对点集进行一系列操作，使用VCPKG包管理工具安装了CGAL库。特此记录一下。

## VCPKG

* [VCPKG](https://github.com/microsoft/vcpkg)是一个由微软开发的跨平台开源软件包管理器。它简化了在Windows、Linux和macOS上获取和管理C/C++库的过程。使用VCPKG，您可以轻松地将第三方库安装并集成到您的C++项目中。
* 要安装VCPKG，您可以按照[VCPKG GitHub repository](https://github.com/microsoft/vcpkg)中提供的说明进行操作. 安装后，可以使用VCPKG命令行界面搜索和安装库，包括CGAL库。

在一个文件夹下，打开bash，或者powershell输入
```bash  
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
```
就这么简单！ vcpkg 已安装并可供使用。

## CGAL
 * CGAL（Computational Geometry Algorithms Library）是一个用于计算几何的开源软件库。
 * 它提供了一系列的算法和数据结构，用于解决各种计算几何问题，如点、线、多边形、曲线、曲面等的计算和操作。
 * CGAL的设计目标是提供高效、可靠和易于使用的计算几何工具，以满足各种应用领域的需求。
 * 该库支持多种编程语言，包括C++、Python和Java，并且可以与其他计算几何库和软件集成使用。
 * CGAL的代码质量和可移植性得到了广泛认可，并且在学术界和工业界都得到了广泛应用。
 * 更多关于CGAL的信息可以在官方网站上找到：https://www.cgal.org/


由于vcpkg for windows中的gmp存在错误，需要安装32位的yasm工具才能正确构建cgal所需的gmp_64位：(当前所用VCPKG的安装方式是经典模式，若要使用清单模式(VCPKG.json)，请仔细阅读[Microsoft vcpkg manifest-mode](https://learn.microsoft.com/zh-cn/vcpkg/consume/manifest-mode))
```bash  
.\vcpkg.exe install yasm-tool:x86-windows
.\vcpkg.exe install cgal:x64-windows
.\vcpkg.exe integrate install 
```

## CMakeLists.txt

这里提供一个我安装时用的测试文件。
```CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
set(CMAKE_TOOLCHAIN_FILE "C:/DEVELOP/vcpkg/scripts/buildsystems/vcpkg.cmake")
project(AlphaShapeExample)

# Find CGAL package
find_package(CGAL CONFIG REQUIRED)

# Set source files
set(SOURCE_FILES ./src/boundaryNode.cpp)

# Add executable
add_executable(AlphaShapeExample ${SOURCE_FILES})

# Link CGAL
target_link_libraries(AlphaShapeExample PRIVATE CGAL::CGAL)
```

## 遇到的问题
* 安装CGAL失败，多次尝试失败的原因是因为[tukaani-project](https://github.com/tukaani-project/xz/archive/v5.4.4.tar.gz)库因为系统环境安全问题在github.com被禁用了。目前的解决方法是：

在如下网址(https://github.com/bminor/xz/archive/refs/tags/v5.4.4.tar.gz)下载会得到一个v5.4.4.tar.gz，然后将其重命名为tukaani-project-xz-v5.4.4.tar.gz，并放到VCPKG的downloads文件夹下。这样vcpkg就不会尝试下载它了。

> * 环境变量 
```bash
VCPKG_ROOT="C:\path\to\vcpkg"
PATH=%VCPKG_ROOT%;%PATH%
```
> * yasm-tool:x86-windows

由于vcpkg for windows中的gmp存在错误，需要安装32位的yasm工具才能正确构建cgal所需的gmp_64位
```bash
.\vcpkg.exe install yasm-tool:x86-windows
```
这样应该可以在Windows，cmake，visual studio2022中使用了。