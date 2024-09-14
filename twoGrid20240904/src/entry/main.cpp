/*----------------------------------------------------------------------------------------------------
* ���������ڼ����ά�򵥶����񣬳��߽�����⣬ȫ���ɲ�ֵõ����޲�ֵ��,����˼��㾫�ȡ�
* Author: Liu Guangying 
* Develop: 2024.9.4-2024.9.10
----------------------------------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <memory>
#include "gridGenerate/include/gridBase.h"
#include "gridGenerate/include/gridBoundS.h"
#include "gridGenerate/include/gridUniform.h"
#include "flowSolve/include/ctrlFlow.h"
#include "flowSolve/include/Solver.h"
#include "triangle/include/triBase.h"

int main() {

	GlobalData::createFolder("result");
	std::cout << "��ʼ���� .... " << std::endl;

	// ʹ�ö�̬�Դ�����ͬ���͵��������
	// -------------------------------------------------------------------------------------------------debug��
	//std::shared_ptr<GridBase> back = std::make_unique<GridUniform>(0.0, 0.6, 0.0, 0.4, 7, 5);
	//std::shared_ptr<GridBase> mgrid = std::make_unique<GridBoundS>(0.208, 0.408, 0.0, 0.4, 3, 5, 1, 1);
	// -------------------------------------------------------------------------------------------------debug��
	//std::shared_ptr<GridBase> back = std::make_unique<GridUniform>(0.0, 1, 0.0, 0.4, 11, 5);
	//std::shared_ptr<GridBase> mgrid = std::make_unique<GridBoundS>(0.208, 0.408, 0.0, 0.4, 3, 5, 1, 1);
	// -------------------------------------------------------------------------------------------------debug��
	std::shared_ptr<GridBase> back = std::make_unique<GridUniform>(0, 2, 0,1, 101, 51);
	std::shared_ptr<GridBase> mgrid = std::make_unique<GridBoundS>(0.307, 0.427, 0, 1, 7, 51, 3, 5);

	CtrlFlow ctrlflow;
	ctrlflow.addGrid(std::move(back));  // ʹ�� std::move ��������ָ��
	ctrlflow.addGrid(std::move(mgrid));
	ctrlflow.solve();
	std::shared_ptr<CtrlFlow> tosolve = std::make_shared<CtrlFlow>(ctrlflow);
	FlowSolver flowsolver(tosolve);
	flowsolver.solve();

	std::cout << "������� .... " << std::endl;

	return 0;
}