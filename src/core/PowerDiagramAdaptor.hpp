#pragma once
#include "Config.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <map>
#include <type_traits>
// the number type
#include <CGAL/MP_Float.h>
// example that uses an exact number type
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>
#include <CGAL/enum.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Weighted_point_2.h>

NAMESPACE_BEGIN(mscut)
namespace core
{
	class PowerDiagramAdaptor
	{
	public:
		// typedefs for the traits and the algorithm
		typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef CGAL::Polygon_2<K>                                  Polygon_2;
		typedef CGAL::Triangle_2<K>									Triangle_2;
		typedef CGAL::Bbox_2										Bbox_2;
		typedef K::Point_2											Point_2;
		typedef K::Point_3											Point_3;
		typedef K::Vector_2											Vector_2;
		typedef K::Vector_3											Vector_3;
		typedef K::Iso_rectangle_2									Iso_rectangle_2;
		typedef K::Segment_2										Segment_2;
		typedef K::Segment_3										Segment_3;
		typedef K::Ray_2											Ray_2;
		typedef K::Line_2											Line_2;

		typedef CGAL::Regular_triangulation_2<K>					Regular_triangulation;
		typedef K::Weighted_point_2									Wp_2;

		/* Member type */
		using PowerDiagramPoint_2 = Point_2;
		using PowerDiagramLine_2 = Segment_2;

		/* A class to recover Voronoi diagram and shrink in a given boundary from stream. */
		struct CroppedVoronoiFromPower {
			std::vector<Segment_2> m_cropSegments;
			Iso_rectangle_2 m_bbox; // m_boundary的bbox
			Polygon_2 m_boundary;

			CroppedVoronoiFromPower() noexcept = default;
			CroppedVoronoiFromPower(const Polygon_2& _boundary) :m_boundary(_boundary) {
				m_bbox = m_boundary.bbox()/* + Iso_rectangle_2(-1, -1, 1, 1)*/; // TODO: 是否有必要将boundingbox扩大一点？
			}
			CroppedVoronoiFromPower(const CroppedVoronoiFromPower& other) {
				m_cropSegments = other.m_cropSegments;
				m_bbox = other.m_bbox;
				m_boundary = other.m_boundary;
			}

			bool isValidSeg(const Segment_2& seg)
			{
				// 一定要先判断是否为nan再判断是否是degenerate
				return (!std::isnan(seg.vertex(0).x()) && !std::isnan(seg.vertex(1).y()) &&
					!std::isnan(seg.vertex(0).x()) && !std::isnan(seg.vertex(1).y()) &&
					!seg.is_degenerate());
			}

			std::vector<Point_2> Intersect(const Segment_2& seg)
			{
				std::vector<Point_2> intersection_points;
				for (const Segment_2& e : m_boundary.edges()) {
					if (CGAL::do_intersect(seg, e))
					{
						const auto intersect_result = CGAL::intersection(seg, e);
						if (const Point_2* p = boost::get<Point_2>(&*intersect_result)) {
							intersection_points.push_back(*p);
						}
					}
				}
				return intersection_points;
			}

			void crop_and_extract_segment(const CGAL::Object& rsl) {
				const Segment_2* s = nullptr;
				// 先与bbox求交，因为传入的rsl可能是Ray_2类型
				if (CGAL::object_cast<Ray_2>(&rsl))
				{
					const Ray_2 ray_2 = *(CGAL::object_cast<Ray_2>(&rsl));
					CGAL::Object obj = CGAL::intersection(ray_2, m_bbox);
					s = CGAL::object_cast<Segment_2>(&obj);
				}
				else if (CGAL::object_cast<Segment_2>(&rsl))
				{
					s = CGAL::object_cast<Segment_2>(&rsl);
				}
				else return;

				// 再与boundary求交
				if (s && isValidSeg(*s))
				{
					const Segment_2 seg = *s;

					const CGAL::Bounded_side side1 = m_boundary.bounded_side(seg.source());
					const CGAL::Bounded_side side2 = m_boundary.bounded_side(seg.target());

					// 都不在外部
					if (side1 != CGAL::ON_UNBOUNDED_SIDE && side2 != CGAL::ON_UNBOUNDED_SIDE)
					{
						m_cropSegments.push_back(seg);
						return;
					}

					std::vector<Point_2> intersection_points = Intersect(seg);

					if (intersection_points.empty()) return;

					Segment_2 clipSeg;
					if (intersection_points.size() == 1)
					{
						// 都不在内部，说明交点是两个外部点相切而来的
						if (side1 != CGAL::ON_BOUNDED_SIDE && side2 != CGAL::ON_BOUNDED_SIDE) return;
						Point_2 innerPoint = seg.source();
						if (side2 == CGAL::ON_BOUNDED_SIDE) innerPoint = seg.target(); // 如果target对应的端点在内部
						clipSeg = Segment_2(innerPoint, intersection_points[0]);
					}
					else
					{
						clipSeg = Segment_2(intersection_points[0], intersection_points[1]);
					}
					if (!isValidSeg(clipSeg)) return;

					m_cropSegments.push_back(clipSeg);

					/*const Segment_2 seg = *s;
					m_cropSegments.push_back(seg);*/
				}
			}

			void reset() {
				m_cropSegments.clear(); m_cropSegments.shrink_to_fit();
			}
		};

	private:
		/* Data */
		CroppedVoronoiFromPower vor;

	public:
		/* Constructor and Destructor */
		PowerDiagramAdaptor() noexcept = default;

		PowerDiagramAdaptor(const PowerDiagramAdaptor& _other) {
			vor = _other.vor;
		}

	public:
		template<class Site, class BD>
		std::vector<Segment_2> computePDForBoundary(const std::vector<Site>& sites/*, const std::vector<bool>& flags*/, const BD& boundary)
		{
			static_assert(std::is_same_v<typename Site::_Point, Point_2>, "POINT_IN_POWER_DIAGRAM_MUST_BE_2-DIM");
			static_assert(std::is_same_v<BD, Polygon_2> || std::is_same_v<BD, Triangle_2>, "TYPE_OF_BOUNDARY_IN_POWER_DIAGRAM_MUST_BE_Polygon_2_OR_Triangle_2");

			std::vector<Segment_2> adCropSegments; // 保存最后的结果

			std::unordered_map<Point_2, int> siteToIdx;

			// iterate to extract Voronoi diagram edges around each vertex
			Regular_triangulation rt;
			for (int i = 0; i < sites.size(); ++i)
			{
				// 扰动项
				K::Vector_2 vec((rand() / (double)RAND_MAX - 0.5) * 0.001,
					(rand() / (double)RAND_MAX - 0.5) * 0.001);

				const Point_2 disturbPos = sites[i].pos + vec;
				Wp_2 site(disturbPos, sites[i].weight);
				//std::cout << "Weight of site #" << i << " = " << sites[i].weight << std::endl;
				siteToIdx.insert(std::make_pair(disturbPos, i));
				rt.insert(site);
			}
			rt.is_valid();

			vor = CroppedVoronoiFromPower(boundary);
			for (Regular_triangulation::Edge_iterator eit = rt.edges_begin(); eit != rt.edges_end(); ++eit) {
				//std::cout << "Vertex " << vit->site().point() << std::endl;
				// ec为与生成的Power Diagram的顶点vit相邻的所有边，边的另一端顶点可能是一个无限远的顶点，所以需要裁剪操作
				//Regular_triangulation::Edge_circulator ec = rt.incident_edges(vit);
				auto site_1 = eit->first->vertex(CGAL::Triangulation_cw_ccw_2::ccw(eit->second))->point().point();
				auto site_2 = eit->first->vertex(CGAL::Triangulation_cw_ccw_2::cw(eit->second))->point().point();

				const double f_val_1 = sites[siteToIdx[site_1]].f_val;
				const double f_val_2 = sites[siteToIdx[site_2]].f_val;

				/*bool isLargeZero_1 = (f_val_1 > 1e-9); // TODO; 等于0的情况
				bool isLargeZero_2 = (f_val_2 > 1e-9);
				bool isEqualZero_1 = (std::fabs(f_val_1) < 1e-9);
				bool isEqualZero_2 = (std::fabs(f_val_2) < 1e-9);
				// 除去都等于0的情况，只要一个等于0或者二者异号就会保留
				if ((isLargeZero_1 ^ isLargeZero_2) ||
					((isEqualZero_1 || isEqualZero_2) && (isEqualZero_1 ^ isEqualZero_2)))
				{
					CGAL::Object o = rt.dual(eit);
					vor.crop_and_extract_segment(o);
				}*/

				bool isLargeZero_1 = (f_val_1 > 0);
				bool isLargeZero_2 = (f_val_2 > 0);
				if (isLargeZero_1 ^ isLargeZero_2)
				{
					CGAL::Object o = rt.dual(eit);
					vor.crop_and_extract_segment(o);
				}

				//print the cropped Voronoi diagram edges as segments
				/*std::copy(vor.m_cropped_vd.begin(), vor.m_cropped_vd.end(),
					std::ostream_iterator<Segment_2>(std::cout, "\n"));
				std::cout << "=========\n";*/

				adCropSegments.insert(adCropSegments.end(), vor.m_cropSegments.begin(), vor.m_cropSegments.end());

				vor.reset();
			}

			return adCropSegments;
		}
	};

} // namespace core
NAMESPACE_END(mscut)