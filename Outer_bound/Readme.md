> ### This is a program that uses VCPKG package management tool to install CGAL library and test the output of 2D point cloud outer edge point set.
> ## Author�� Liu Guangying
> ## Email:   liugy36@mail2.sysu.edu.cn 

# ����
- ����Ŀ��Ҫ�ǲ���ʹ��VCPKG����ģʽ��װCGAL�⣬����2D����ģ�͵����Ե�㼯��ע�ⲻ��͹���㼯��ʹ��cgal���е�delunay�㷨��alpha-shape�㷨��
- ���Ե���ģ�ͷ�����./model/points.txt�С�
- �ļ���./src�����������Դ��룬���ǿ����õģ�ʹ��ʱֻ��Ҫ��CMakeLists.txt�ļ����޸�����build���ɣ�
```cpp
set(SOURCE_FILES ./src/***.cpp)
```

20240407
��������Ϊ��Ҫ�Ե㼯����һϵ�в�����ʹ��VCPKG�������߰�װ��CGAL�⡣�ش˼�¼һ�¡�

## VCPKG

* [VCPKG](https://github.com/microsoft/vcpkg)��һ����΢�����Ŀ�ƽ̨��Դ���������������������Windows��Linux��macOS�ϻ�ȡ�͹���C/C++��Ĺ��̡�ʹ��VCPKG�����������ɵؽ��������ⰲװ�����ɵ�����C++��Ŀ�С�
* Ҫ��װVCPKG�������԰���[VCPKG GitHub repository](https://github.com/microsoft/vcpkg)���ṩ��˵�����в���. ��װ�󣬿���ʹ��VCPKG�����н��������Ͱ�װ�⣬����CGAL�⡣

��һ���ļ����£���bash������powershell����
```bash  
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
```
����ô�򵥣� vcpkg �Ѱ�װ���ɹ�ʹ�á�

## CGAL
 * CGAL��Computational Geometry Algorithms Library����һ�����ڼ��㼸�εĿ�Դ����⡣
 * ���ṩ��һϵ�е��㷨�����ݽṹ�����ڽ�����ּ��㼸�����⣬��㡢�ߡ�����Ρ����ߡ�����ȵļ���Ͳ�����
 * CGAL�����Ŀ�����ṩ��Ч���ɿ�������ʹ�õļ��㼸�ι��ߣ����������Ӧ�����������
 * �ÿ�֧�ֶ��ֱ�����ԣ�����C++��Python��Java�����ҿ������������㼸�ο���������ʹ�á�
 * CGAL�Ĵ��������Ϳ���ֲ�Եõ��˹㷺�Ͽɣ�������ѧ����͹�ҵ�綼�õ��˹㷺Ӧ�á�
 * �������CGAL����Ϣ�����ڹٷ���վ���ҵ���https://www.cgal.org/


����vcpkg for windows�е�gmp���ڴ�����Ҫ��װ32λ��yasm���߲�����ȷ����cgal�����gmp_64λ��(��ǰ����VCPKG�İ�װ��ʽ�Ǿ���ģʽ����Ҫʹ���嵥ģʽ(VCPKG.json)������ϸ�Ķ�[Microsoft vcpkg manifest-mode](https://learn.microsoft.com/zh-cn/vcpkg/consume/manifest-mode))
```bash  
.\vcpkg.exe install yasm-tool:x86-windows
.\vcpkg.exe install cgal:x64-windows
.\vcpkg.exe integrate install 
```

## CMakeLists.txt

�����ṩһ���Ұ�װʱ�õĲ����ļ���
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

## ����������
* ��װCGALʧ�ܣ���γ���ʧ�ܵ�ԭ������Ϊ[tukaani-project](https://github.com/tukaani-project/xz/archive/v5.4.4.tar.gz)����Ϊϵͳ������ȫ������github.com�������ˡ�Ŀǰ�Ľ�������ǣ�

��������ַ(https://github.com/bminor/xz/archive/refs/tags/v5.4.4.tar.gz)���ػ�õ�һ��v5.4.4.tar.gz��Ȼ����������Ϊtukaani-project-xz-v5.4.4.tar.gz�����ŵ�VCPKG��downloads�ļ����¡�����vcpkg�Ͳ��᳢���������ˡ�

> * �������� 
```bash
VCPKG_ROOT="C:\path\to\vcpkg"
PATH=%VCPKG_ROOT%;%PATH%
```
> * yasm-tool:x86-windows

����vcpkg for windows�е�gmp���ڴ�����Ҫ��װ32λ��yasm���߲�����ȷ����cgal�����gmp_64λ
```bash
.\vcpkg.exe install yasm-tool:x86-windows
```
����Ӧ�ÿ�����Windows��cmake��visual studio2022��ʹ���ˡ�