#pragma once
#include "Config.hpp"
#include "utils/Log.hpp"
#include "utils/String.hpp"
#include "detail/BasicDataType.hpp"
#include <string>
#include <iomanip>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangle_2.h>

typedef  CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef  K::Point_2											  Point_2;
typedef  K::Segment_2										  Segment_2;
typedef  K::Point_3											  Point_3;
typedef  K::Segment_3										  Segment_3;
typedef  CGAL::Triangle_2<K>								  Triangle_2;

NAMESPACE_BEGIN(mscut)
using namespace detail;

namespace geometry::gvis
{
	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ofstream& output, const Point_2& v)
	{
		output << "v " << v.x() << " " << v.y() << " " << 0 << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ofstream& output, const Point_3& v)
	{
		output << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ofstream& output, const Point_2& v, const Vector3& rgb)
	{
		output << "v " << v.x() << " " << v.y() << " " << 0 << " " << rgb.x() << " " << rgb.y() << " " << rgb.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ofstream& output, const Point_3& v, const Vector3& rgb)
	{
		output << "v " << v.x() << " " << v.y() << " " << v.z() << " " << rgb.x() << " " << rgb.y() << " " << rgb.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex_to_xyz(std::ofstream& output, const Point_2& v)
	{
		output << v.x() << " " << v.y() << " " << 0 << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex_to_xyz(std::ofstream& output, const Point_3& v)
	{
		output << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}

	// Helper function to write line
	static void write_line(std::ofstream& output, int& vertIdx, const Segment_2& l)
	{
		if (vertIdx <= 0) vertIdx = 1;

		write_vertex(output, l.vertex(0));
		write_vertex(output, l.vertex(1));
		output << "l " << vertIdx << " " << vertIdx + 1 << std::endl;

		vertIdx += 2;
	}

	// Helper function to write line
	static void write_line(std::ofstream& output, int& vertIdx, const Segment_3& l)
	{
		if (vertIdx <= 0) vertIdx = 1;

		write_vertex(output, l.vertex(0));
		write_vertex(output, l.vertex(1));
		output << "l " << vertIdx << " " << vertIdx + 1 << std::endl;

		vertIdx += 2;
	}

	static void write_triangle_2_to_obj(std::ofstream& output, int& vertIdx, const Triangle_2& triangle_2)
	{
		if (vertIdx <= 0) vertIdx = 1;

		write_vertex(output, triangle_2.vertex(0));
		write_vertex(output, triangle_2.vertex(1));
		write_vertex(output, triangle_2.vertex(2));

		output << "f " << vertIdx << " " << vertIdx + 1 << " " << vertIdx + 2 << std::endl;

		vertIdx += 3;
	}

} // namespace geometry::gvis

NAMESPACE_END(mscut)