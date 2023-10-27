#pragma once
#include "../CDTAdaptor.hpp"
#include <Config.hpp>
#include <geometry/CGALGeometry.hpp>

NAMESPACE_BEGIN(mscut)
namespace unit_test
{
	using namespace core;

	inline bool testComputeCDT(/*const std::string& mesh_file*/)
	{
		struct CDT_BD {
			using _Point = Point_2;
			using _Edge = Segment_2;

			std::vector<_Point> bdPoints;
			std::vector<_Edge> bdEdges;

			CDT_BD() noexcept = default;

			CDT_BD(const std::vector<_Edge>& _bdEdges) :bdEdges(_bdEdges) {}
		};

		std::vector<typename CDT_BD::_Point> _bdPoints;
		_bdPoints.emplace_back(Point_2(0, 0));
		_bdPoints.emplace_back(Point_2(1, 0));
		_bdPoints.emplace_back(Point_2(0, 1));
		_bdPoints.emplace_back(Point_2(1, 1));

		//保证逆时针
		if (CGAL::orientation(_bdPoints[0], _bdPoints[1], _bdPoints[2]) != CGAL::POSITIVE)
		{
			std::reverse(_bdPoints.begin(), _bdPoints.end());
		}

		// 设置边界边
		int numBDPoints = _bdPoints.size();
		std::vector<typename CDT_BD::_Edge> _bdEdges;
		for (int i = 0; i < _bdPoints.size(); ++i)
			_bdEdges.emplace_back(CDT_BD::_Edge(_bdPoints[i], _bdPoints[(i + 1) % numBDPoints]));

		CDT_BD cdtBD(_bdEdges);

		CDTAdaptor cdtAdaptor;

		const auto& cdt_res = cdtAdaptor.computeCDForBoundary(cdtBD);

		const std::string vis_file = R"(.\src\core\test\cdt_test.obj)";
		str_util::checkDir(vis_file);
		std::ofstream vis_out(vis_file);
		if (!vis_out) { LOG::qpError("I/O: File ", vis_file.c_str(), " could not be opened!"); return false; }

		int faceIdx = 1;
		for (const Triangle_2& res_triangle_2 : cdt_res)
		{
			geometry::gvis::write_triangle_2_to_obj(vis_out, faceIdx, res_triangle_2);
		}

		vis_out.close();

		return true;
	}

} // namespace unit_test
NAMESPACE_END(mscut)