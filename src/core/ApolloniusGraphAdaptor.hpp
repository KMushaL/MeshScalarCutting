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
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>
#include <CGAL/enum.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Triangle_2<K>									Triangle_2;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K::Vector_2 Vector_2;
typedef K::Vector_3 Vector_3;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

// typedefs for the traits and the algorithm
typedef CGAL::Apollonius_graph_traits_2<K>   Traits;
typedef CGAL::Apollonius_graph_2<Traits>     Apollonius_graph;

NAMESPACE_BEGIN(mscut)
namespace core
{
	class ApolloniusGraphAdaptor
	{
	public:
		/* Member type */
		using ApolloniusDiagramPoint_2 = Point_2;
		using ApolloniusDiagramLine_2 = Segment_2;

		//A class to recover Voronoi diagram from stream.
		struct CroppedVoronoiFromApollonius {
			std::list<Segment_2> m_cropped_vd;
			Iso_rectangle_2 m_bbox;

			CroppedVoronoiFromApollonius() noexcept = default;
			CroppedVoronoiFromApollonius(const Iso_rectangle_2& bbox) :m_bbox(bbox) {}
			CroppedVoronoiFromApollonius(const CroppedVoronoiFromApollonius& other) {
				m_cropped_vd = other.m_cropped_vd;
				m_bbox = other.m_bbox;
			}

			template <class RSL>
			void crop_and_extract_segment(const RSL& rsl) {
				CGAL::Object obj = CGAL::intersection(rsl, m_bbox);
				const Segment_2* s = CGAL::object_cast<Segment_2>(&obj);
				if (s) m_cropped_vd.push_back(*s);
			}

			void operator<<(const Ray_2& ray) { crop_and_extract_segment(ray); }
			void operator<<(const Line_2& line) { crop_and_extract_segment(line); }
			void operator<<(const Segment_2& seg) { crop_and_extract_segment(seg); }

			void reset() {
				m_cropped_vd.erase(m_cropped_vd.begin(), m_cropped_vd.end());
			}
		};

		/* shrink操作，将appollonius diagram中那些跑出边界(三角形)的边强制保留在边界内部 */
		struct BoundaryApolloniusGraph
		{
			Polygon_2 m_boundary;
			std::vector<Segment_2> clipSegments; // �ü�������б�

			BoundaryApolloniusGraph() noexcept = default;
			BoundaryApolloniusGraph(const Polygon_2& boundary) noexcept : m_boundary(boundary) {}
			BoundaryApolloniusGraph(const BoundaryApolloniusGraph& other) {
				m_boundary = other.m_boundary;
				clipSegments = other.clipSegments;
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

			void operator<<(const Segment_2& seg)
			{
				if (seg.is_degenerate()) return;
				CGAL::Bounded_side side1 = m_boundary.bounded_side(seg.source());
				CGAL::Bounded_side side2 = m_boundary.bounded_side(seg.target());

				// ����������ⲿֱ�����
				if (side1 != CGAL::ON_UNBOUNDED_SIDE && side2 != CGAL::ON_UNBOUNDED_SIDE)
				{
					clipSegments.push_back(seg);
					return;
				}
				std::vector<Point_2> intersection_points = Intersect(seg);
				if (intersection_points.empty()) return;

				Segment_2 clipSeg;
				if (intersection_points.size() == 1)
				{
					if (side1 != CGAL::ON_BOUNDED_SIDE && side2 != CGAL::ON_BOUNDED_SIDE) return; // �����������ڵ�ֱ��return
					Point_2 innerPoint = seg.source();
					if (side2 == CGAL::ON_BOUNDED_SIDE) innerPoint = seg.target();
					clipSeg = Segment_2(innerPoint, intersection_points[0]);
				}
				else
				{
					clipSeg = Segment_2(intersection_points[0], intersection_points[1]);
				}
				if (clipSeg.is_degenerate()) return;// �˻���

				clipSegments.push_back(clipSeg);
			}

			void operator<<(const K::Ray_2& ray)
			{
				std::cerr << "cannot be dealt\n";
			}

			void operator<<(const K::Line_2& line)
			{
				std::cerr << "cannot be dealt\n";
			}
		};

	private:
		/* Data */
		CroppedVoronoiFromApollonius vor;

	public:
		/* Constructor and Destructor */
		ApolloniusGraphAdaptor() {
			Iso_rectangle_2 bbox(-2000, -2000, 2000, 2000);
			vor.m_bbox = bbox;
		}

		ApolloniusGraphAdaptor(const ApolloniusGraphAdaptor& _other) {
			vor = _other.vor;
		}

	public:
		template<class Site, class BD>
		std::vector<Segment_2> computeAGForBoundary(const std::vector<Site>& sites, const BD& boundary)
		{
			static_assert(std::is_same_v<typename Site::_Point, Point_2>, "POINT_IN_APOLLONIUS_GRAPH_MUST_BE_2-DIM");
			static_assert(std::is_same_v<BD, Polygon_2> || std::is_same_v<BD, Triangle_2>, "TYPE_OF_BOUNDARY_IN_APOLLONIUS_GRAPH_MUST_BE_Polygon_2_OR_Triangle_2");

			// iterate to extract Voronoi diagram edges around each vertex
			Apollonius_graph ag;
			for (int i = 0; i < sites.size(); ++i)
			{
				Apollonius_graph::Site_2 site(sites[i].pos, sites[i].weight);
				ag.insert(site);
			}

			BoundaryApolloniusGraph bag(boundary);
			for (auto vit = ag.finite_vertices_begin(); vit != ag.finite_vertices_end(); ++vit) {
				//std::cout << "Vertex " << vit->site().point() << std::endl;
				Apollonius_graph::Edge_circulator ec = ag.incident_edges(vit), done(ec);
				if (ec != 0) {
					do {
						ag.draw_dual_edge(*ec, vor);
					} while (++ec != done);
				}
				//print the cropped Voronoi diagram edges as segments
				std::copy(vor.m_cropped_vd.begin(), vor.m_cropped_vd.end(),
					std::ostream_iterator<Segment_2>(std::cout, "\n"));
				std::cout << "=========\n";
				for (const auto& vor_seg : vor.m_cropped_vd)
				{
					bag << vor_seg;
				}

				vor.reset();
			}

			std::cout << 333 << std::endl;

			return bag.clipSegments;
		}
	};

} // namespace core
NAMESPACE_END(mscut)