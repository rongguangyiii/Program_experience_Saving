#include "grid/include/gridDisc.h"
#include "grid/include/globalData.h"
#include <iostream>
#include <filesystem>  

void DiscGrid::generatePoints()
{
	const double R0 = 1;
	const double R1 = 5;
	size_t count = 0;
	for (size_t i = 0; i < x_points_; i++)
	{
		for (size_t j = 0; j < y_points_; j++)
		{
			double coef = (std::exp(2.0 * i / (x_points_ - 1)) - 1.0) / (std::exp(2) - 1.0);//圆盘非均匀网格
			//double coef = (double)i / (x_points_ - 1);//圆盘均匀网格
			double R = R1 - j * (R1 - R0) / (y_points_ - 1);
			double theta = GlobalData::PI * coef;
			double x = -R * sin(theta);
			double y = R * cos(theta);
			count = i * y_points_ + j;
			points_.emplace_back(std::make_shared<Coord>(PointId(level_, count), x, y));
			this->x[i][j] = x;
			this->y[i][j] = y;
		}
	}

}

void DiscGrid::generateElements() {
	for (size_t i = 0; i < x_points_ - 1; ++i) {
		for (size_t j = 0; j < y_points_ - 1; ++j) {
			size_t p1 = i * y_points_ + j + 1;
			size_t p2 = p1 + 1;
			size_t p3 = p1 + y_points_;
			size_t p4 = p3 + 1;
			elements_.push_back({ p1, p3, p4, p2 });
		}
	}
}

void DiscGrid::genNeiborNode()
{
	size_t i = 0, j = 0;
	size_t xNodeNum = x_points_, yNodeNum = y_points_;
	for (size_t iNode = 0; iNode < points_.size(); ++iNode)
	{
		auto& currentNode = points_.at(iNode);
		PointType currentNodeType = PointType::inner;
		std::vector<CoordPtr> neiborNode;
		i = iNode / yNodeNum;
		j = iNode % yNodeNum;
		std::vector<size_t> neiborNodeIndex = { iNode + yNodeNum, iNode + 1, iNode - yNodeNum, iNode - 1 };//设置邻居节点的索引，右，上，左，下
		if (j == 0)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 3);
			currentNodeType = PointType::b_filedbound;
		}
		if (i == 0)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 2);
			currentNodeType = PointType::b_filedbound;
		}
		if (j == yNodeNum - 1)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin() + 1);
			currentNodeType = PointType::b_filedbound;;
		}
		if (i == xNodeNum - 1)
		{
			neiborNodeIndex.erase(neiborNodeIndex.begin());
			currentNodeType = PointType::b_filedbound;
		}
		for (const auto it : neiborNodeIndex)
			neiborNode.emplace_back(points_[it]);
		currentNode->type_ = currentNodeType;
		Reneighbors_[currentNode->id_] = neiborNode;
	}
}

void DiscGrid::outputToTecplot(std::ofstream& out, const std::string& zone_title) const {
	out << "TITLE = \"" << zone_title << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << points_.size() << ", E=" << elements_.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";

	for (const auto& point : points_) {
		out << std::fixed << std::setprecision(6) << point->x << " " << point->y << "\n";
	}

	for (size_t i = 0; i < elements_.size(); i++)
	{
		const auto& eleconn = elements_[i].getConnectivty();
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


	std::cout << zone_title << " output complete!\n";
}

void DiscGrid::outputToTecplot_VCFEM(std::ofstream& out, const std::string& zone_title) const {
	out << "TITLE = \"" << zone_title << "\"\n";
	out << "VARIABLES = \"X\", \"Y\"\n";
	out << "ZONE T=\"" << zone_title << "\", N=" << points_.size() << ", E=" << elements_.size() << ", F=FEPOINT, ET=QUADRILATERAL\n";

	//1. CONNECTIVITY
	out << "CONNECTIVITY \n";
	for (size_t i = 0; i < elements_.size(); i++)
	{
		const auto& eleconn = elements_[i].getConnectivty();
		const size_t num = eleconn.size();
		if (eleconn.size() == 3) {
			std::cout << "error!!!\n";
		}
		else if (eleconn.size() == 4) {
			out << i + 1 << " " << num << " " << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[3] << " " << -1 << " " << 1 << "\n";
		}
		else {
			std::cout << "error in outplot!\n";
		}
	}
	//2. COORDINATES
	out << "COORDINATES \n";
	size_t nn = 0;
	for (const auto& point : points_) {
		nn++;
		out << std::fixed << std::setprecision(6) << nn << " " << point->x << " " << point->y << "\n";
	}
	//3. INCLUSIONS
	out << "INCLUSIONS \n";
	for (size_t i = 0; i < elements_.size(); i++) {

		out << i + 1 << " " << 0 << "\n";
	}

	out << "INFORMATION AND VALUES OF PRESCRIBED INCREMENT DISPLACEMENT \n";

	//左边x方向固定u=0-------------------------------------------------------------
	out << 1 << ",    " << 1 << ",    " << y_points_ << "\n";
	out << 0 << "\n";
	for (size_t i = 0; i < y_points_; i++) {
		out << i + 1 << ",    ";
	}
	out << std::endl;
	//左边y方向固定一个点v=0-------------------------------------------------------------
	out << 2 << ",    " << 2 << ",    " << 1 << "\n";
	out << 0 << "\n";
	out << 1 << "\n";
	//左边x方向均匀拉伸u=0.01-------------------------------------------------------------
	out << 3 << ",    " << 1 << ",    " << y_points_ << "\n";
	out << 0.01 << "\n";
	size_t sta = (x_points_ - 1) * y_points_;
	for (size_t i = 0; i < y_points_; i++) {
		out << sta + i + 1 << ",    ";
	}

	std::cout << zone_title << " output complete!\n";
}

//重载操作符=
std::shared_ptr<GridBase> DiscGrid::operator=(const std::shared_ptr<GridBase>& other) {
	points_ = other->getPoints();
	elements_ = other->getelements();
	Reneighbors_ = other->getReneighbors();
	x_points_ = other->getXnum();
	y_points_ = other->getYnum();
	x_start_ = other->getXstart();
	x_end_ = other->getXend();
	y_start_ = other->getYstart();
	y_end_ = other->getYend();
	return std::make_shared<DiscGrid>(*this);
}

void DiscGrid::coefTrans() {
	coeftransVec_.clear();
	double ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, jacob = 0;
	double XK = 0, XI = 0, YK = 0, YI = 0, JJ = 0;
	const size_t mx = x_points_;
	const size_t my = x_points_;
	size_t index = 0;
	for (size_t i = 0; i < mx; ++i) {
		for (size_t j = 0; j < my; ++j) {
			if (i == 0 && j != 0 && j != my - 1) {  // 1. 左边界 (i == 0) 非角点 (0 < j < my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == mx - 1 && j != 0 && j != my - 1) {  // 2. 右边界 (i == mx) 非角点 (0 < j < my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == 0 && j == 0) {  // 3. 左下角 (i == 0, j == 0)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (i == 0 && j == my - 1) {  // 4. 左上角 (i == 0, j == my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == my - 1) {  // 5. 右上角 (i == mx, j == my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == 0) {  // 6. 右下角 (i == mx, j == 0)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (j == 0) {  // 7. 下边界 (j == 0) 非角点 (0 < i < mx)
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else if (j == my - 1) {  // 8. 上边界 (j == my) 非角点 (0 < i < mx)
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else {  // 9. 内部点
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			index=i*y_points_+j;
			JJ = XK * YI - XI * YK;  // 雅可比行列式
			ksi_x = YI / JJ;
			ksi_y = -XI / JJ;
			eta_x = -YK / JJ;
			eta_y = XK / JJ;
			jacob = 1.0 / JJ;
			coeftransVec_.emplace_back(CoefTrans(ksi_x,ksi_y,eta_x,eta_y,jacob));
		}
	}
}