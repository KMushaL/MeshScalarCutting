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

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

NAMESPACE_BEGIN(mscut)
namespace core
{
	class ApolloniusGraphAdaptor
	{
	public:
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
		typedef K::Ray_2											Ray_2;
		typedef K::Line_2											Line_2;

		// typedefs for the traits and the algorithm
		typedef CGAL::Apollonius_graph_traits_2<K>					Traits;
		typedef CGAL::Apollonius_graph_2<Traits>					Apollonius_graph;

		/* Member type */
		using ApolloniusDiagramPoint_2 = Point_2;
		using ApolloniusDiagramLine_2 = Segment_2;

		/* A class to recover Voronoi diagram and shrink in a given boundary from stream. */
		struct CroppedVoronoiFromApollonius {
			std::vector<Segment_2> m_cropSegments;
			Iso_rectangle_2 m_bbox; // m_boundary的bbox
			Polygon_2 m_boundary;

			CroppedVoronoiFromApollonius() noexcept = default;
			CroppedVoronoiFromApollonius(const Polygon_2& _boundary) :m_boundary(_boundary) {
				m_bbox = m_boundary.bbox()/* + Bbox_2(-1, -1, 1, 1)*/; // TODO: 是否有必要将boundingbox扩大一点？
			}
			CroppedVoronoiFromApollonius(const CroppedVoronoiFromApollonius& other) {
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
							// ���result��Point_2���͵�ָ�룬��˵��ֻ��һ������
							intersection_points.push_back(*p);
						}
						//else if (const Segment_2* s = boost::get<Segment_2>(&*intersect_result)) {
						//	// ���result��Segment_2���͵�ָ�룬��˵����һ���߶���Ϊ����
						//}
					}
				}
				/*Segment_2 clipSeg(seg);
				if (intersection_points.size() == 1)
				{
					clipSeg = Segment_2(seg[0], intersection_points[0]);
				}
				else if (intersection_points.size() == 2)
				{
					clipSeg = Segment_2(intersection_points[0], intersection_points[1]);
				}*/
				return intersection_points;
			}

			template <class RSL>
			void crop_and_extract_segment(const RSL& rsl) {
				// 先与bbox求交，因为传入的rsl可能是Ray_2类型
				CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
				const Segment_2* s = CGAL::object_cast<Segment_2>(&obj);
				// 再与boundary求交
				if (s && isValidSeg(*s))
				{
					const Segment_2 seg = *s;

					CGAL::Bounded_side side1 = m_boundary.bounded_side(seg.source());
					CGAL::Bounded_side side2 = m_boundary.bounded_side(seg.target());

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
						if (side1 != CGAL::ON_BOUNDED_SIDE && side2 != CGAL::ON_BOUNDED_SIDE) return;
						Point_2 innerPoint = seg.source();
						if (side2 == CGAL::ON_BOUNDED_SIDE) innerPoint = seg.target();
						clipSeg = Segment_2(innerPoint, intersection_points[0]);
					}
					else
					{
						clipSeg = Segment_2(intersection_points[0], intersection_points[1]);
					}
					if (!isValidSeg(clipSeg)) return;

					m_cropSegments.push_back(clipSeg);
				}
			}

			void operator<<(const Ray_2& ray) { crop_and_extract_segment(ray); }
			void operator<<(const Line_2& line) { crop_and_extract_segment(line); }
			void operator<<(const Segment_2& seg) { crop_and_extract_segment(seg); }

			void reset() {
				m_cropSegments.clear(); m_cropSegments.shrink_to_fit();
			}
		};

	private:
		/* Data */
		CroppedVoronoiFromApollonius vor;

	public:
		/* Constructor and Destructor */
		ApolloniusGraphAdaptor() noexcept = default;

		ApolloniusGraphAdaptor(const ApolloniusGraphAdaptor& _other) {
			vor = _other.vor;
		}

	public:
		template<class Site, class BD>
		std::vector<Segment_2> computeAGForBoundary(const std::vector<Site>& sites/*, const std::vector<bool>& flags*/, const BD& boundary)
		{
			static_assert(std::is_same_v<typename Site::_Point, Point_2>, "POINT_IN_APOLLONIUS_GRAPH_MUST_BE_2-DIM");
			static_assert(std::is_same_v<BD, Polygon_2> || std::is_same_v<BD, Triangle_2>, "TYPE_OF_BOUNDARY_IN_APOLLONIUS_GRAPH_MUST_BE_Polygon_2_OR_Triangle_2");

			std::vector<Segment_2> adCropSegments; // 保存最后的结果

			std::unordered_map<Point_2, int> siteToIdx;

			// iterate to extract Voronoi diagram edges around each vertex
			Apollonius_graph ag;
			for (int i = 0; i < sites.size(); ++i)
			{
				// 扰动项
				K::Vector_2 vec((rand() / (double)RAND_MAX - 0.5) * 0.00001,
					(rand() / (double)RAND_MAX - 0.5) * 0.00001);

				const Point_2 disturbPos = sites[i].pos + vec;
				Apollonius_graph::Site_2 site(disturbPos, sites[i].weight);
				siteToIdx.insert(std::make_pair(disturbPos, i));
				ag.insert(site);
			}

			vor = CroppedVoronoiFromApollonius(boundary);
			//BoundaryApolloniusGraph bag(boundary);
			for (auto vit = ag.finite_vertices_begin(); vit != ag.finite_vertices_end(); ++vit) {
				//std::cout << "Vertex " << vit->site().point() << std::endl;
				// ec为与生成的Apollonius Diagram的顶点vit相邻的所有边，边的另一端顶点可能是一个无限远的顶点，所以需要裁剪操作
				Apollonius_graph::Edge_circulator ec = ag.incident_edges(vit), done(ec);
				if (ec != 0) {
					do {
						auto site_1 = ec->first->vertex(CGAL::Triangulation_cw_ccw_2::ccw(ec->second))->site().point();
						auto site_2 = ec->first->vertex(CGAL::Triangulation_cw_ccw_2::cw(ec->second))->site().point();

						bool isLargeZero_1 = (sites[siteToIdx[site_1]].f_val > 0); // TODO; 等于0的情况
						bool isLargeZero_2 = (sites[siteToIdx[site_2]].f_val > 0);
						if (isLargeZero_1 ^ isLargeZero_2) ag.draw_dual_edge(*ec, vor);
					} while (++ec != done);
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