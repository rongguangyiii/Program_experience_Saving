#pragma once
#include <string>

namespace GlobalData
{
	const std::string meshType = "Uniform";//�������� Uniform; Semicircular;
	const std::string calmethod = "Common";//����RHS���� Deer(����ϵ��); Common(ʹ�ü�������ϵ��); Basic(��ͨ������)
	const size_t nx = 101;
	const size_t ny = 51;
	const double cfl = 0.9;      // CFL��
	const double ctrltime = 100;      // ��ʱ��
	const size_t ctrlstep = 20000;      // ʱ�䲽��
	const double tolerance = 1e-15;  // ������ֵ
	const size_t ctrlout = 100; //�������
	const size_t movestep = 1;
	const double PI = 3.141592653589793238462643383279502884197169399;//Բ����pi
	void createFolder(const std::string& basePath);
	void outGrid();
}