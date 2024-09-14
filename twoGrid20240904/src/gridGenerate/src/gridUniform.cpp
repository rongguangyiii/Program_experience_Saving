
#include "gridGenerate/include/gridUniform.h"


void GridUniform::generatePoints() {
	double x_step = (x_end_ - x_start_) / (x_points_ - 1);
	double y_step = (y_end_ - y_start_) / (y_points_ - 1);
	for (size_t i = 0; i < x_points_; ++i) {
		for (size_t j = 0; j < y_points_; ++j) {
			size_t num = i * y_points_ + j;
			x[i][j] = x_start_ + i * x_step;
			y[i][j] = y_start_ + j * y_step;
			points_.emplace_back(std::make_shared<Coord>(PointId(level_, num), x_start_ + i * x_step, y_start_ + j * y_step,0 ));
		}
	}
}

void GridUniform::generateElements() {
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
void GridUniform::genNeiborNode()
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
		std::vector<size_t> neiborNodeIndex = { iNode + yNodeNum, iNode + 1, iNode - yNodeNum, iNode - 1 };//�����ھӽڵ���������ң��ϣ�����
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

		//if ((i == xNodeNum - 1 || i == 0) && (j == 0 || j == yNodeNum - 1))
		//{
		//	currentNodeType = Node::Nodetype::corner;
		//}
		//currentNode->setType(currentNodeType);
		for (const auto it : neiborNodeIndex)
			neiborNode.emplace_back(points_[it]);
		currentNode->type_ = currentNodeType;
		Reneighbors_[currentNode->id_] = neiborNode;

	}
}

void GridUniform::outputToTecplot(std::ofstream& out, const std::string& zone_title) const {
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
	//for (const auto& elem : elements_) {
	//	const auto& eleconn = elem.getConnectivty();
	//	if (eleconn.size() == 3)
	//	{
	//		out << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[2] << "\n";
	//	}
	//	else if (eleconn.size() == 4)
	//	{
	//		out << eleconn[0] << " " << eleconn[1] << " " << eleconn[2] << " " << eleconn[3] << "\n";
	//	}
	//	else
	//	{
	//		std::cout << "error in outplot!\n";
	//	}
	//}

	std::cout << zone_title << " output complete!\n";
}

//���ز�����=
std::shared_ptr<GridBase> GridUniform::operator=(const std::shared_ptr<GridBase>& other) {
	points_ = other->getPoints();
	elements_ = other->getelements();
	//neighbors_ = other->getneighbors();
	Reneighbors_ = other->getReneighbors();
	x_points_ = other->getXnum();
	y_points_ = other->getYnum();
	x_start_ = other->getXstart();
	x_end_ = other->getXend();
	y_start_ = other->getYstart();
	y_end_ = other->getYend();
	return std::make_shared<GridUniform>(*this);
}

void GridUniform::coordTrans() {
	double XK, XI, YK, YI, JJ;
	double ksi_t = 0, ksi_x = 0, ksi_y = 0, eta_x = 0, eta_y = 0, eta_t = 0;
	double jacob = 0.0;
	const size_t mx = x_points_;
	const size_t my = y_points_;
	for (size_t i = 0; i < mx; ++i) {
		for (size_t j = 0; j < my; ++j) {
			if (i == 0 && j != 0 && j != my - 1) {  // 1. ��߽� (i == 0) �ǽǵ� (0 < j < my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == mx - 1 && j != 0 && j != my - 1) {  // 2. �ұ߽� (i == mx) �ǽǵ� (0 < j < my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}
			else if (i == 0 && j == 0) {  // 3. ���½� (i == 0, j == 0)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (i == 0 && j == my - 1) {  // 4. ���Ͻ� (i == 0, j == my)
				XK = 0.5 * (4.0 * x[i + 1][j] - x[i + 2][j] - 3.0 * x[i][j]);
				YK = 0.5 * (4.0 * y[i + 1][j] - y[i + 2][j] - 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == my - 1) {  // 5. ���Ͻ� (i == mx, j == my)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
			}
			else if (i == mx - 1 && j == 0) {  // 6. ���½� (i == mx, j == 0)
				XK = 0.5 * (-4.0 * x[i - 1][j] + x[i - 2][j] + 3.0 * x[i][j]);
				YK = 0.5 * (-4.0 * y[i - 1][j] + y[i - 2][j] + 3.0 * y[i][j]);
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
			}
			else if (j == 0) {  // 7. �±߽� (j == 0) �ǽǵ� (0 < i < mx)
				XI = 0.5 * (4.0 * x[i][j + 1] - x[i][j + 2] - 3.0 * x[i][j]);
				YI = 0.5 * (4.0 * y[i][j + 1] - y[i][j + 2] - 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else if (j == my - 1) {  // 8. �ϱ߽� (j == my) �ǽǵ� (0 < i < mx)
				XI = 0.5 * (-4.0 * x[i][j - 1] + x[i][j - 2] + 3.0 * x[i][j]);
				YI = 0.5 * (-4.0 * y[i][j - 1] + y[i][j - 2] + 3.0 * y[i][j]);
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
			}
			else {  // 9. �ڲ���
				XK = (x[i + 1][j] - x[i - 1][j]) / 2.0;
				YK = (y[i + 1][j] - y[i - 1][j]) / 2.0;
				XI = (x[i][j + 1] - x[i][j - 1]) / 2.0;
				YI = (y[i][j + 1] - y[i][j - 1]) / 2.0;
			}

			JJ = XK * YI - XI * YK;  // �ſɱ�����ʽ

			ksi_x = YI / JJ;
			ksi_y = -XI / JJ;
			eta_x = -YK / JJ;
			eta_y = XK / JJ;
			jacob = 1.0 / JJ;
			coefTrans_[points_[i*y_points_+j]->id_]= CoefTrans(ksi_x, ksi_y, eta_x, eta_y, jacob);
		}
	}
}
