#pragma once
#include <iostream>
#include <fstream>
#include <cassert>
#include <map>
#include <igl/read_triangle_mesh.h>


// the number type
#include <CGAL/MP_Float.h>


// example that uses an exact number type

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iterator>
#include <CGAL/enum.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;

// typedefs for the traits and the algorithm

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

typedef CGAL::Apollonius_graph_traits_2<K>   Traits;
typedef CGAL::Apollonius_graph_2<Traits>     Apollonius_graph;


//A class to recover Voronoi diagram from stream.
struct Cropped_voronoi_from_apollonius {
	std::list<Segment_2> m_cropped_vd;
	Iso_rectangle_2 m_bbox;

	Cropped_voronoi_from_apollonius(const Iso_rectangle_2& bbox) :m_bbox(bbox) {}

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

struct BoundaryApolloniusGraph
{
	Polygon_2 m_boundary;
	std::vector<Segment_2> clipSegments; // 裁减后的所有边

	BoundaryApolloniusGraph(const Polygon_2& boundary)
		: m_boundary(boundary) {}

	std::vector<Point_2> Intersect(const Segment_2& seg)
	{
		std::vector<Point_2> intersection_points;
		for (const Segment_2& e : m_boundary.edges()) {
			if (CGAL::do_intersect(seg, e))
			{
				const auto intersect_result = CGAL::intersection(seg, e);
				if (const Point_2* p = boost::get<Point_2>(&*intersect_result)) {
					// 如果result是Point_2类型的指针，则说明只有一个交点
					intersection_points.push_back(*p);
				}
				//else if (const Segment_2* s = boost::get<Segment_2>(&*intersect_result)) {
				//	// 如果result是Segment_2类型的指针，则说明有一条线段作为交点
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

		// 如果都不在外部直接输出
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
			if (side1 != CGAL::ON_BOUNDED_SIDE && side2 != CGAL::ON_BOUNDED_SIDE) return; // 两个都不在内的直接return
			Point_2 innerPoint = seg.source();
			if (side2 == CGAL::ON_BOUNDED_SIDE) innerPoint = seg.target();
			clipSeg = Segment_2(innerPoint, intersection_points[0]);
		}
		else
		{
			clipSeg = Segment_2(intersection_points[0], intersection_points[1]);
		}
		if (clipSeg.is_degenerate()) return;// 退化点

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

