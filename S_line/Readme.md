> ## Name :CDT(Constrain Delaunay Triangle)约束delaunay三角形
> ## Author :Guangying Liu 
> ## Email: EatShark1208@163.com  
> ## Date: 2024.04.23

## CDT(Constrain Delaunay Triangle)约束delaunay三角形

> **该项目，用于生成约束delaunay三角形，并进行三角形的处理，挑选和删除。在头文件中定义了三个类分别是TriEdge，TriEle，TriBase，用于存储和输出数据。结果可以通过输出的tecplot网格文件查看。主要用初始三角形网格，筛选三角形网格，和重构三角形三个过程的网格。**

> **后续增加了三角形和四边形混合的网格，主要用于适应多种情况**

* 需要链接库：CGAL
## Windows 下编译
- 编译环境：Win 10，Win 11系统； Visual Studio 2019/2022
### [0] 环境变量设置，方便后续库文件的查找链接
> * 在windows系统环境变量中添加以下环境变量：  
        VCPKG_ROOT="安装vcpkg所在目录"  
        VCPKG_DEFAULT_TRIPLET="x64-windows" 
> * 在变量为***path***的环境变量中，添加以下变量(此步骤也可以省略)：  
        PATH=%VCPKG_ROOT%;%PATH%
### [1] VCPKG
* [VCPKG](https://github.com/microsoft/vcpkg)是一个由微软开发的跨平台开源软件包管理器。它简化了在Windows、Linux和macOS上获取和管理C/C++库的过程。使用VCPKG，您可以轻松地将第三方库安装并集成到您的C++项目中。
* 要安装VCPKG，您可以按照[VCPKG GitHub repository](https://github.com/microsoft/vcpkg)中提供的说明进行操作. 安装后，可以使用VCPKG命令行界面搜索和安装库，包括CGAL库。

在任意目录下，打开bash或者powershell输入：
```bash  
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
```
就这么简单！ vcpkg 已安装并可供使用。
### [2] CGAL
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
### [3] CMakeLists.txt

这里提供一个我安装时用的***CMakeLists.txt***测试文件。
```CMake
cmake_minimum_required(VERSION 3.0)
set(CMAKE_TOOLCHAIN_FILE "D:/dev_VCPKG/vcpkg/scripts/buildsystems/vcpkg.cmake")
project(TestDemo)

# Find CGAL package
find_package(CGAL CONFIG REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include) 

# Organize source files in a dedicated directory variable
set(SOURCE_DIR src)
set(SOURCE_FILES
    ${SOURCE_DIR}/demo.cpp 
    ${SOURCE_DIR}/filter_Tri.cpp 
    # other cpp file....
)

# Create an executable with the specified sources
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Link the CGAL target to the executable
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)

```
### [3] Problem
如果上面步骤都安装成功，以下部分就不用看了！
* 在windows系统下,安装的库都要集成到Visual Studio中。这里应该注意的是：
    > VS的安装路径不能有中文字符，空格或连字符。否则会导致库查找MSVC编译器配置对应链接库时出错。
    
这样应该可以在Windows系统下的visual studio2022中使用了。