> ## Name :CDT(Constrain Delaunay Triangle)Լ��delaunay������
> ## Author :Guangying Liu 
> ## Email: EatShark1208@163.com  
> ## Date: 2024.04.23

## CDT(Constrain Delaunay Triangle)Լ��delaunay������

> **����Ŀ����������Լ��delaunay�����Σ������������εĴ�����ѡ��ɾ������ͷ�ļ��ж�����������ֱ���TriEdge��TriEle��TriBase�����ڴ洢��������ݡ��������ͨ�������tecplot�����ļ��鿴����Ҫ�ó�ʼ����������ɸѡ���������񣬺��ع��������������̵�����**

> **���������������κ��ı��λ�ϵ�������Ҫ������Ӧ�������**

* ��Ҫ���ӿ⣺CGAL
## Windows �±���
- ���뻷����Win 10��Win 11ϵͳ�� Visual Studio 2019/2022
### [0] �����������ã�����������ļ��Ĳ�������
> * ��windowsϵͳ����������������»���������  
        VCPKG_ROOT="��װvcpkg����Ŀ¼"  
        VCPKG_DEFAULT_TRIPLET="x64-windows" 
> * �ڱ���Ϊ***path***�Ļ��������У�������±���(�˲���Ҳ����ʡ��)��  
        PATH=%VCPKG_ROOT%;%PATH%
### [1] VCPKG
* [VCPKG](https://github.com/microsoft/vcpkg)��һ����΢�����Ŀ�ƽ̨��Դ���������������������Windows��Linux��macOS�ϻ�ȡ�͹���C/C++��Ĺ��̡�ʹ��VCPKG�����������ɵؽ��������ⰲװ�����ɵ�����C++��Ŀ�С�
* Ҫ��װVCPKG�������԰���[VCPKG GitHub repository](https://github.com/microsoft/vcpkg)���ṩ��˵�����в���. ��װ�󣬿���ʹ��VCPKG�����н��������Ͱ�װ�⣬����CGAL�⡣

������Ŀ¼�£���bash����powershell���룺
```bash  
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
```
����ô�򵥣� vcpkg �Ѱ�װ���ɹ�ʹ�á�
### [2] CGAL
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
### [3] CMakeLists.txt

�����ṩһ���Ұ�װʱ�õ�***CMakeLists.txt***�����ļ���
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
������沽�趼��װ�ɹ������²��־Ͳ��ÿ��ˣ�
* ��windowsϵͳ��,��װ�ĿⶼҪ���ɵ�Visual Studio�С�����Ӧ��ע����ǣ�
    > VS�İ�װ·�������������ַ����ո�����ַ�������ᵼ�¿����MSVC���������ö�Ӧ���ӿ�ʱ����
    
����Ӧ�ÿ�����Windowsϵͳ�µ�visual studio2022��ʹ���ˡ�