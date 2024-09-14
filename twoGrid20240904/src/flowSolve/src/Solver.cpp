#include "flowSolve/include/Solver.h"
#include "flowSolve/include/globalData.h"
#include <iostream>
#include <fstream>


void FlowSolver::solve()
{
	initialflow();
	timeAdvance();
}

void FlowSolver::initialflow()
{
	std::vector<std::shared_ptr<GridBase>>& gridVec = ctrlflow_->getGridVec();
	for (const auto& curgrid : gridVec)
	{ 
		size_t xpoints = curgrid->getXnum();
		size_t ypoints = curgrid->getYnum();
		auto& points = curgrid->getPoints();
		auto& neighbors = curgrid->getReneighbors();
		CoordPtr cp = nullptr;
		for (size_t i = 0; i < xpoints; i++)
		{
			for (size_t j = 0; j < ypoints; j++)
			{
				cp = points[i * ypoints + j];
				if (cp->type_ == PointType::b_notcal)
					continue;
				u[cp->id_] = 1.0;
			}
		}
	}
	u_new = u;
}

void FlowSolver::FilterPointsToCal(const bool isfirst)
{
	if(!isfirst)
	{
		calPoints_.clear();
		AllNeighbors_.clear();
	}
	std::vector<std::shared_ptr<GridBase>>& gridVec = ctrlflow_->getGridVec();
	//step 1. 从网格0中找到要计算的点
	auto& curgrid = gridVec.at(0);
	size_t xpoints = curgrid->getXnum();
	size_t ypoints = curgrid->getYnum();
	auto& points = curgrid->getPoints();
	auto& neighbors = curgrid->getReneighbors();
	CoordPtr cp = nullptr;
	size_t n1 = 0;
	for (size_t i = 1; i < xpoints - 1; i++)
	{
		for (size_t j = 1; j < ypoints - 1; j++)
		{
			cp = points[i * ypoints + j];
			if (cp->type_ == PointType::b_notcal)
				continue;
			calPoints_.push_back(cp);
			n1++;
			if (neighbors.find(cp->id_) == neighbors.end())
				std::cout << "errr 110\n";
			AllNeighbors_[cp->id_] = neighbors[cp->id_];
		}
	}

	//step 1. 从网格1中提取要计算的点
	auto& mgrid = gridVec.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();
	auto& m_neighbors = mgrid->getReneighbors();
	size_t n2 = 0;
	for (size_t i = 0; i < mx; i++)
	{
		for (size_t j = 1; j < my - 1; j++)
		{
			cp = mPointsvec[i * my + j];
			if (cp->type_ == PointType::b_notcal)
				continue;
			calPoints_.push_back(cp);
			n2++;
			if (m_neighbors.find(cp->id_) == m_neighbors.end())
				std::cout << "errr 110\n";
			AllNeighbors_[cp->id_] = m_neighbors[cp->id_];
		}
	}

}

void FlowSolver::FilterPointsToCal()
{
	if (isfirst_)
	{
		FilterPointsToCal(isfirst_);
		saveOldNeighbors();//因为第一次计算需要使用所以，这里复制一次。
		isfirst_ = false;
	}
	else
	{
		if (ctrlflow_->isSpit_ || ctrlflow_->isSwallow_)
		{
			FilterPointsToCal(isfirst_);
		}
		else
		{
			//没有吞吐，do Nothing!
		}
	}

}

void FlowSolver::saveOldNeighbors()
{
	OldNeighbors_.clear();
	std::vector<Coord>temp;
	for (const auto& perVec : AllNeighbors_)
	{
		for (const auto& perit : perVec.second)
		{
			temp.push_back(*perit);
		}
		OldNeighbors_[perVec.first] = temp;
		temp.clear();
	}

}

void FlowSolver::timeAdvance()
{
	size_t timestep = 0;
	writeTecplotFile(timestep);
	while (timestep < 150/*residual > GlobalData::tolerance*/)
	{
		++timestep;
		FilterPointsToCal();
		if (caltime > GlobalData::ctrltime /*||*/ /*residual < GlobalData::tolerance*/ /*timestep > GlobalData::ctrlstep*/)
			break;
		//1.欧拉时间推进
		solverEuler();
		//2.当前步残差计算
		computeResidual();
		//3.将u更新为u_new，准备进入下一步
		caltime += dt;
		u = u_new;
		writeTecplotFile(timestep);
		std::cout << "Timestep: " << timestep << " Residual: " << residual << std::endl;
		//4.移动网格
		if (ctrlflow_->isSpit_ || ctrlflow_->isSwallow_)
			saveOldNeighbors();
		ctrlflow_->move(1.0, dt);
	}
}


void FlowSolver::solverEuler()
{
	//step 1. 计算右端项
	computeRHS_deer();
	//computeRHS_common();

	//step 2. 计算u_new
	auto& coef_n = ctrlflow_->getallcoefTrans();//使用移动后新网格计算的coef
	for (const auto& it : calPoints_) {
		if (it->tag_ == PointTag::spit)
			continue;
		Tools::checkKeyExists(u, it->id_, "u");
		Tools::checkKeyExists(rhs, it->id_, "rhs");
		Tools::checkKeyExists(coef_n, it->id_, "coef_n");
		u_new[it->id_] = u[it->id_] + dt * coef_n[it->id_].jacob() * rhs[it->id_];
	}

	//step 3. 处理spit网格点
	if (ctrlflow_->isSpit_)
		calSpitPoints();

	//step 4. 处理边界点
	flowBoundary();

}

void FlowSolver::calSpitPoints()
{
	//step 1. 计算右端项
	computeSpitPointsRHS_deer();
	//computeRHS_common();

	//step 2. 计算u_new
	auto& coef_n = ctrlflow_->getSpitcoefTrans();//n时刻度量系数
	auto& spitPoint = ctrlflow_->getSpitoutPoints();
	auto& conPoints = ctrlflow_->getConnectPoints();//n时刻的点，未移动之前的点
	/*这里本质上其实是用n时刻与spit点直接相连的moveGrid上的点的值更新n+1时刻的点*/
	for (size_t i = 0; i < spitPoint.size(); i++)
	{
		const auto& it = conPoints[i];
		Tools::checkKeyExists(u, it->id_, "u");
		Tools::checkKeyExists(rhs, it->id_, "rhs");
		Tools::checkKeyExists(coef_n, it->id_, "coef_n");
		u_new[spitPoint[i]->id_] = u[it->id_] + dt * coef_n[it->id_].jacob() * spit_rhs[it->id_];
		//spitPoint[i]->tag_ = PointTag::surface;
	}
	
}

void FlowSolver::flowBoundary() 
{
	//step 3. 边界外推处理
	extrapBoundary_s();
	extrapBoundary_b();
}

void FlowSolver::extrapBoundary_s()
{
	std::vector<std::shared_ptr<GridBase>>& gridVec = ctrlflow_->getGridVec();
	const auto& it = gridVec.at(1);
	size_t numx = it->getXnum();
	size_t numy = it->getYnum();
	const auto& points = it->getPoints();

	CoordPtr p1 = nullptr;
	CoordPtr p2 = nullptr;
	//对于上边界的处理
	for (size_t i = 0; i < numx; i++) {
		p1 = points[i * numy + numy - 1];
		p2 = points[i * numy + numy - 2];
		if (p1->type_ == PointType::b_notcal)
			continue;
		if (u_new.find(p2->id_) == u_new.end())
			std::cout << "errr 110\n";
		u_new[p1->id_] = u_new[p2->id_];
		//u_new[i][numy - 1] = u_new[i][numy - 2];
	}
	//均匀网格左侧入口边界条件，对于下边界的处理
	for (size_t i = 0; i < numx; i++) {
		p1 = points[i * numy + 0];
		p2 = points[i * numy + 1];
		if (p1->type_ == PointType::b_notcal)
			continue;
		if (u_new.find(p2->id_) == u_new.end())
			std::cout << "errr 110\n";
		u_new[p1->id_] = u_new[p2->id_];
		//u_new[i][0] = u_new[i][1];
	}

}

void FlowSolver::extrapBoundary_b()
{
	std::vector<std::shared_ptr<GridBase>>& gridVec = ctrlflow_->getGridVec();
	const auto& it = gridVec.at(0);

	size_t numx = it->getXnum();
	size_t numy = it->getYnum();
	const auto& points = it->getPoints();

	CoordPtr p1 = nullptr;
	CoordPtr p2 = nullptr;
	//右边界
	for (size_t j = 1; j < numy - 1; j++) {
		p1 = points[(numx - 1) * numy + j];
		p2 = points[(numx - 2) * numy + j];
		if (p1->type_ == PointType::b_notcal)
			continue;
		if (u_new.find(p2->id_) == u_new.end())
			std::cout << "errr 110\n";
		u_new[p1->id_] = u_new[p2->id_];
		//u_new[nx - 1][j] = u_new[nx - 2][j];

	}
	//均匀网格左侧入口边界条件，对于z左边界的处理
	for (size_t j = 1; j < numy - 1; j++) {
		p1 = points[0 * numy + j];
		p2 = points[1 * numy + j];
		if (p1->type_ == PointType::b_notcal)
			continue;
		if (u_new.find(p2->id_) == u_new.end())
			std::cout << "errr 110\n";
		u_new[p1->id_] = u_new[p2->id_];
		//u_new[0][j] = u_new[1][j];
	}
	//均匀网格左侧入口边界条件，对于上边界的处理
	for (size_t i = 0; i < numx; i++) {
		p1 = points[i * numy + numy - 1];
		p2 = points[i * numy + numy - 2];
		if (p1->type_ == PointType::b_notcal)
			continue;
		if (u_new.find(p2->id_) == u_new.end())
			std::cout << "errr 110\n";
		u_new[p1->id_] = u_new[p2->id_];
		//u_new[i][numy - 1] = u_new[i][numy - 2];
	}
	//均匀网格左侧入口边界条件，对于下边界的处理
	for (size_t i = 0; i < numx; i++) {
		p1 = points[i * numy + 0];
		p2 = points[i * numy + 1];
		if (p1->type_ == PointType::b_notcal)
			continue;
		if (u_new.find(p2->id_) == u_new.end())
			std::cout << "errr 110\n";
		u_new[p1->id_] = u_new[p2->id_];
		//u_new[i][0] = u_new[i][1];
	}

}

void FlowSolver::computeRHS_deer()
{
	auto& coef_n = ctrlflow_->getallcoefTrans();//n时刻度量系数

	//采用DEER冻结系数和strger-warming方法计算流场
	rhs.clear();

	for (const auto& it : calPoints_) {
		if (it->tag_ == PointTag::spit)
			continue;
		Tools::checkKeyExists(coef_n, it->id_, "coef_n");
		const auto& jacobi = coef_n[it->id_].jacob();
		const auto& ksi_x = coef_n[it->id_].ksi_x();
		const auto& ksi_y = coef_n[it->id_].ksi_y();
		const auto& eta_x = coef_n[it->id_].eta_x();
		const auto& eta_y = coef_n[it->id_].eta_y();
		const auto& ksi_t = coef_n[it->id_].ksi_t();
		const auto& eta_t = coef_n[it->id_].eta_t();
		auto perneighbors = Tools::cood2CoodPtr(OldNeighbors_[it->id_]);//n时刻的邻居关系。

		CoordPtr& pR = GlobalData::PointAtAngle(perneighbors, 0);
		CoordPtr& pU = GlobalData::PointAtAngle(perneighbors, GlobalData::PI / 2);
		CoordPtr& pL = GlobalData::PointAtAngle(perneighbors, GlobalData::PI);
		CoordPtr& pD = GlobalData::PointAtAngle(perneighbors, GlobalData::PI * 3 / 2);

		// x 方向通量分裂
		double Fx_plus = SW_positive_ksi(u[it->id_], jacobi, ksi_x, ksi_y, ksi_t);
		double Fx_minus = SW_negative_ksi(u[it->id_], jacobi, ksi_x, ksi_y, ksi_t);
		double Fx_plus_left = SW_positive_ksi(u[pL->id_], jacobi, ksi_x, ksi_y, ksi_t);
		double Fx_minus_right = SW_negative_ksi(u[pR->id_], jacobi, ksi_x, ksi_y, ksi_t);

		// y 方向通量分裂
		double Fy_plus = SW_positive_eta(u[it->id_], jacobi, eta_x, eta_y, eta_t);
		double Fy_minus = SW_negative_eta(u[it->id_], jacobi, eta_x, eta_y, eta_t);
		double Fy_plus_down = SW_positive_eta(u[pD->id_], jacobi, eta_x, eta_y, eta_t);
		double Fy_minus_up = SW_negative_eta(u[pU->id_], jacobi, eta_x, eta_y, eta_t);

		rhs[it->id_] = (-1.0) * ((Fx_plus - Fx_plus_left + Fx_minus_right - Fx_minus) + (Fy_plus - Fy_plus_down + Fy_minus_up - Fy_minus));
	}

}

void FlowSolver::computeSpitPointsRHS_deer()
{
	auto& coef_n = ctrlflow_->getSpitcoefTrans();//n时刻度量系数
	//auto& spitPoint = ctrlflow_->getSpitoutPoints();
	auto& conPoints = ctrlflow_->getConnectPoints();//moveGrid上的n时刻的点，

	spit_rhs.clear();

	for (const auto& it : conPoints) {//moveGrid上的n时刻的点，是把边界点反向推进到吐出的点上
		Tools::checkKeyExists(coef_n, it->id_, "coef_n");
		const auto& jacobi = coef_n[it->id_].jacob();
		const auto& ksi_x = coef_n[it->id_].ksi_x();
		const auto& ksi_y = coef_n[it->id_].ksi_y();
		const auto& eta_x = coef_n[it->id_].eta_x();
		const auto& eta_y = coef_n[it->id_].eta_y();
		const auto& ksi_t = coef_n[it->id_].ksi_t();
		const auto& eta_t = coef_n[it->id_].eta_t();
		auto perneighbors = Tools::cood2CoodPtr(OldNeighbors_[it->id_]);//n时刻的邻居关系。

		CoordPtr& pR = GlobalData::PointAtAngle(perneighbors, 0);
		CoordPtr& pU = GlobalData::PointAtAngle(perneighbors, GlobalData::PI / 2);
		CoordPtr& pL = GlobalData::PointAtAngle(perneighbors, GlobalData::PI);
		CoordPtr& pD = GlobalData::PointAtAngle(perneighbors, GlobalData::PI * 3 / 2);
												
		// x 方向通量分裂
		double Fx_plus = SW_positive_ksi(u[it->id_], jacobi, ksi_x, ksi_y, ksi_t);
		double Fx_minus = SW_negative_ksi(u[it->id_], jacobi, ksi_x, ksi_y, ksi_t);
		double Fx_plus_left = SW_positive_ksi(u[pL->id_], jacobi, ksi_x, ksi_y, ksi_t);
		double Fx_minus_right = SW_negative_ksi(u[pR->id_], jacobi, ksi_x, ksi_y, ksi_t);

		// y 方向通量分裂
		double Fy_plus = SW_positive_eta(u[it->id_], jacobi, eta_x, eta_y, eta_t);
		double Fy_minus = SW_negative_eta(u[it->id_], jacobi, eta_x, eta_y, eta_t);
		double Fy_plus_down = SW_positive_eta(u[pD->id_], jacobi, eta_x, eta_y, eta_t);
		double Fy_minus_up = SW_negative_eta(u[pU->id_], jacobi, eta_x, eta_y, eta_t);

		spit_rhs[it->id_] = (-1.0) * ((Fx_plus - Fx_plus_left + Fx_minus_right - Fx_minus) + (Fy_plus - Fy_plus_down + Fy_minus_up - Fy_minus));
	}

}

//void FlowSolver::computeRHS_common()
//{
//	const auto& jacobi = grid_->getJacobi_old();
//	const auto& ksi_x = grid_->getKsix_old();
//	const auto& ksi_y = grid_->getKsiy_old();
//	const auto& eta_x = grid_->getEtax_old();
//	const auto& eta_y = grid_->getEtay_old();
//	const auto& ksi_t = grid_->getKsit_old();
//	const auto& eta_t = grid_->getEtat_old();
//
//	//采用DEER冻结系数和strger-warming方法计算流场
//	for (auto& row : this->rhs) {
//		std::fill(row.begin(), row.end(), 0.0);
//	}
//
//	for (size_t i = 1; i < nx - 1; ++i) {
//		for (size_t j = 1; j < ny - 1; ++j) {
//			// x 方向通量分裂
//			double Fx_plus = SW_positive_ksi(u[i][j], jacobi[i][j], ksi_x[i][j], ksi_y[i][j], ksi_t[i][j]);
//			double Fx_minus = SW_negative_ksi(u[i][j], jacobi[i][j], ksi_x[i][j], ksi_y[i][j], ksi_t[i][j]);
//			double Fx_plus_left = SW_positive_ksi(u[i - 1][j], jacobi[i][j], ksi_x[i][j], ksi_y[i][j], ksi_t[i][j]);
//			double Fx_minus_right = SW_negative_ksi(u[i + 1][j], jacobi[i][j], ksi_x[i][j], ksi_y[i][j], ksi_t[i][j]);
//
//			// y 方向通量分裂
//			double Fy_plus = SW_positive_eta(u[i][j], jacobi[i][j], eta_x[i][j], eta_y[i][j], eta_t[i][j]);
//			double Fy_minus = SW_negative_eta(u[i][j], jacobi[i][j], eta_x[i][j], eta_y[i][j], eta_t[i][j]);
//			double Fy_plus_down = SW_positive_eta(u[i][j - 1], jacobi[i][j], eta_x[i][j], eta_y[i][j], eta_t[i][j]);
//			double Fy_minus_up = SW_negative_eta(u[i][j + 1], jacobi[i][j], eta_x[i][j], eta_y[i][j], eta_t[i][j]);
//
//			rhs[i][j] = (-1.0) * ((Fx_plus - Fx_plus_left + Fx_minus_right - Fx_minus) + (Fy_plus - Fy_plus_down + Fy_minus_up - Fy_minus));
//
//		}
//	}
//}

void FlowSolver::computeResidual()
{
	residual = 0.0;
	for (const auto& it : calPoints_) {
		if (it->tag_ == PointTag::spit){
			it->tag_ = PointTag::surface;
			continue;
		}
		Tools::checkKeyExists(rhs, it->id_, "rhs");
		residual += std::pow(rhs[it->id_], 2);
	}

	if (ctrlflow_->isSpit_) {
		auto& conPoints = ctrlflow_->getConnectPoints();
		for (const auto& it : conPoints) {
			Tools::checkKeyExists(rhs, it->id_, "rhs");
			residual += std::pow(spit_rhs[it->id_], 2);
		}
	}
	residual = std::sqrt(residual);
}



void FlowSolver::writeTecplotFile(size_t timestep) {
	//if (timestep % GlobalData::ctrlout != 0)
	//	return;
	std::ofstream outFile("./result/solution_" + std::to_string(timestep) + ".dat");
	outplot_B(outFile, "zone1", timestep);
	outplot_M(outFile, "zone2", timestep);
	std::cout << "第 " << timestep << " 步结果输出完成;" << std::endl;
	outFile.close();
}

void FlowSolver::outplot_B(std::ofstream& out, const std::string& zone_title, const size_t timestep)
{
	std::vector<std::shared_ptr<GridBase>>& gridVec_ = ctrlflow_->getGridVec();

	//step1 先把要输出的点删选出来
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	size_t mx = mgrid->getXnum();
	size_t my = mgrid->getYnum();

	auto& bgrid = gridVec_.at(0);
	size_t bx = bgrid->getXnum();
	size_t by = bgrid->getYnum();
	auto& bPointsVec = bgrid->getPoints();
	std::vector<CoordPtr> outPoints;

	//step 1 确定编号
	size_t LB = 0;
	size_t RB = 0;
	for (size_t i = 0; i < bx; i++) {
		if (bPointsVec[i * by + 0]->type_ == PointType::b_Lbound)
			LB = i;
		if (bPointsVec[i * by + 0]->type_ == PointType::b_Rbound)
			RB = i;
	}
	//step 2 左半边
	for (size_t i = 0; i <= LB; i++) {
		for (size_t j = 0; j < by; j++)
		{
			const auto& it = bPointsVec[i * by + j];
			outPoints.emplace_back(it);
		}
	}
	//step 3 添加movegrid中的左边界，
	for (size_t j = 0; j < my; j++)
	{
		const auto& it = mPointsvec[0 * my + j];
		outPoints.emplace_back(it);
	}
	//step 4 添加movegrid中的右边界，
	for (size_t j = 0; j < my; j++)
	{
		const auto& it = mPointsvec[(mx - 1) * my + j];
		outPoints.emplace_back(it);
	}
	//step 4 右半边，
	for (size_t i = RB; i < bx; i++) {
		for (size_t j = 0; j < by; j++)
		{
			const auto& it = bPointsVec[i * by + j];
			outPoints.emplace_back(it);
		}
	}

	//step 5 形成单元
	std::vector<std::vector<size_t>> outelements;
	for (size_t i = 0; i < LB + 1; i++)
	{
		for (size_t j = 0; j < by - 1; j++)
		{
			size_t p1 = i * by + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + by;
			size_t p4 = p3 + 1;
			outelements.push_back({ p1, p3, p4, p2 });
		}

	}

	size_t outRow = 0;
	size_t row_num = bx - 1 - RB;
	for (size_t i = 0; i < row_num + 1; i++)
	{
		outRow = LB + 2 + i;
		for (size_t j = 0; j < by - 1; j++)
		{
			size_t p1 = outRow * by + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + by;
			size_t p4 = p3 + 1;
			outelements.push_back({ p1, p3, p4, p2 });
		}

	}

	out << "TITLE = \"" << zone_title << "\"\n";
	out << "Variables=\"X\", \"Y\", \"U\", \"ref\", \"err\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << outPoints.size() << ", E=" << outelements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << ", SOLUTIONTIME=" << timestep << std::endl;

	for (const auto& point : outPoints) {
		out << std::fixed << std::setprecision(6) << point->x << " " << point->y << " " << u[point->id_] << " " << 0 << " " << 0 << "\n";
	}

	for (const auto& elem : outelements) {
		out << elem[0] << " " << elem[1] << " " << elem[2] << " " << elem[3] << "\n";
	}

	//std::cout << zone_title << " output complete!\n";
}

void FlowSolver::outplot_M(std::ofstream& out, const std::string& zone_title, const size_t timestep) {
	std::vector<std::shared_ptr<GridBase>>& gridVec_ = ctrlflow_->getGridVec();

	//step1 先把要输出的点删选出来
	auto& mgrid = gridVec_.at(1);
	const auto& mPointsvec = mgrid->getPoints();
	const auto& outelements = mgrid->getelements();

	out << "TITLE = \"" << zone_title << "\"\n";
	out << "Variables=\"X\", \"Y\", \"U\", \"ref\", \"err\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << mPointsvec.size() << ", E=" << outelements.size() << ", F=FEPOINT, ET=QUADRILATERAL" << ", SOLUTIONTIME=" << timestep << std::endl;

	for (const auto& point : mPointsvec) {
		out << std::fixed << std::setprecision(6) << point->x << " " << point->y << " " << u[point->id_] << " " << 0 << " " << 0 << "\n";
	}

	for (const auto& elem : outelements) {
		const auto& eleconn = elem.getConnectivty();
		if (eleconn.size() == 3)
		{
			out << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[2] << "\n";
		}
		else if (eleconn.size() == 4)
		{
			out << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[3] << "\n";
		}
		else
		{
			std::cout << "error in outplot!\n";
		}
	}
	//std::cout << zone_title << " output complete!\n";
}


//=========================================================================================================
double FlowSolver::SW_positive_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t)
{
	//return 0.5 * (f + fabs(f));
	double f_plus_i = 0.5 * (1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)
		+ std::fabs(1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)));
	return f_plus_i;
}

double FlowSolver::SW_negative_ksi(double u, double jacobi, double ksi_x, double ksi_y, double ksi_t)
{
	//return 0.5 * (f - fabs(f));
	double f_minus_i = 0.5 * (1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)
		- std::abs(1 / jacobi * (u * ksi_t + u * ksi_x + (4.0 / 3.0) * u * ksi_y)));
	return f_minus_i;
}

double FlowSolver::SW_positive_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t)
{
	//return 0.5 * (g + fabs(g));
	double g_plus_j = 0.5 * (1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)
		+ std::fabs(1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)));
	return g_plus_j;
}

double FlowSolver::SW_negative_eta(double u, double jacobi, double eta_x, double eta_y, double eta_t)
{
	//return 0.5 * (g - fabs(g));
	double g_minus_j = 0.5 * (1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)
		- std::fabs(1 / jacobi * (u * eta_t + u * eta_x + (4.0 / 3.0) * u * eta_y)));
	return g_minus_j;
}
