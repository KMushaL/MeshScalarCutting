#pragma once
#include "Config.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <map>
#include <list>
#include <type_traits>
// the number type
#include <CGAL/MP_Float.h>
// example that uses an exact number type
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>

// typedefs for the traits and the algorithm
struct FaceInfo_2
{
	FaceInfo_2() noexcept = default;

	int level;
	bool isInDomain()
	{
		return level % 2 == 1; // 非边界面
	}
};
typedef  CGAL::Exact_predicates_inexact_constructions_kernel                        CDT_Kernel;
typedef  CGAL::Triangulation_vertex_base_2<CDT_Kernel>                              CDT_Vb;
typedef  CGAL::Triangulation_face_base_with_info_2<FaceInfo_2, CDT_Kernel>          CDT_Fbb;
typedef  CGAL::Constrained_triangulation_face_base_2<CDT_Kernel, CDT_Fbb>           CDT_Fb;
typedef  CGAL::Triangulation_data_structure_2<CDT_Vb, CDT_Fb>                       CDT_TDS;
typedef  CGAL::Exact_predicates_tag                                                 CDT_Itag;
typedef  CGAL::Constrained_Delaunay_triangulation_2<CDT_Kernel, CDT_TDS, CDT_Itag>  CDT;
typedef  CDT::Point																    Point_2;
typedef  CDT::Edge																    Edge;
typedef  CDT::Face_handle														    Face_handle;
typedef  CDT_Kernel::Segment_2														Segment_2;
typedef  CGAL::Triangle_2<CDT_Kernel>												Triangle_2;


NAMESPACE_BEGIN(mscut)
namespace core
{
	class CDTAdaptor
	{
	public:
		/* Member type */
		using BoundaryPoint = Point_2;
		using BoundaryEdge = Segment_2;

		using CDTPoint = typename BoundaryPoint;
		using CDTEdge = typename BoundaryEdge;

	private:
		/* Data */
		CDT cdt;

	public:
		/* Constructor and Destructor */
		CDTAdaptor() noexcept = default;

		CDTAdaptor(const CDTAdaptor& _other) { cdt = _other.cdt; }

	private:
		void markDomains(CDT::Face_handle start,
			int index,
			std::list<CDT::Edge>& border)
		{
			if (start->info().level != -1) return;

			std::list<CDT::Face_handle> faceList;
			faceList.push_back(start);
			while (!faceList.empty())
			{
				CDT::Face_handle fh = faceList.front();
				faceList.pop_front();
				if (fh->info().level == -1)
				{
					fh->info().level = index;
					for (int i = 0; i < 3; i++)
					{
						CDT::Edge i_edge(fh, i); // 当前face的第i条边
						CDT::Face_handle i_neighbor_face = fh->neighbor(i);

						if (i_neighbor_face->info().level == -1)
						{
							if (cdt.is_constrained(i_edge)) border.push_back(i_edge);
							else faceList.push_back(i_neighbor_face);
						}
					}
				}
			}
		}

		// 通过给定的constrained point计算出CDT后，利用以下函数得到CDT的border edge(可以理解为constrained edge)，即变量cdtBorderEdge
		void markDomains()
		{
			for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it)
			{
				it->info().level = -1;
			}

			std::list<CDT::Edge> cdtBorderEdge;
			markDomains(cdt.infinite_face(), 0, cdtBorderEdge); // 从infinite face开始蔓延搜寻border edge
			while (!cdtBorderEdge.empty())
			{
				CDT::Edge e = cdtBorderEdge.front();
				cdtBorderEdge.pop_front();
				CDT::Face_handle neighborFace = e.first->neighbor(e.second);
				if (neighborFace->info().level == -1)
					markDomains(neighborFace, e.first->info().level + 1, cdtBorderEdge);
			}
		}

	public:
		template<class BD>
		std::vector<Triangle_2> computeCDForBoundary(const BD& boundary)
		{
			static_assert(std::is_same_v<typename BD::_Point, BoundaryPoint>, "TYPE_OF_BOUNDARY_POINT_FOR_CDT_MUST_BE_Point_2");
			static_assert(std::is_same_v<typename BD::_Edge, BoundaryEdge>, "TYPE_OF_BOUNDARY_EDGE_FOR_CDT_MUST_BE_Segment_2");

			cdt.clear();
			std::vector<Triangle_2> cdt_results; // 保存最后的结果

			// insert constrained points
			// make sure the boundary performs a closed-form behavior
			for (const auto& bdEdge : boundary.bdEdges)
				cdt.insert_constraint(bdEdge.vertex(0), bdEdge.vertex(1));
			//cdt.insert_constraint(boundary.bdPoints.begin(), boundary.bdPoints.end(), true);
			cdt.is_valid();
			//markDomains();

			for (const Face_handle& f : cdt.finite_face_handles()) // 遍历所有的三角形
			{
				//if (f->info().isInDomain()) // 如果不是边界面
				{
					CDTPoint a = f->vertex(0)->point();
					CDTPoint b = f->vertex(1)->point();
					CDTPoint c = f->vertex(2)->point();

					Triangle_2 tri_2(a, b, c);

					cdt_results.emplace_back(tri_2);
				}
			}

			return cdt_results;
		}
	};

} // namespace core
NAMESPACE_END(mscut)