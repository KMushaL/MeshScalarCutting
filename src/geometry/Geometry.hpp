#pragma once
#include "Config.hpp"
#include "utils/Log.hpp"
#include "utils/String.hpp"
#include "detail/BasicDataType.hpp"
#include <string>
#include <iomanip>
#include <fstream>

NAMESPACE_BEGIN(mscut)
using namespace detail;

namespace geometry::gvis
{
	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ofstream& output, const Vector3& v)
	{
		output << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex(std::ofstream& output, const Vector3& v, const Vector3& rgb)
	{
		output << "v " << v.x() << " " << v.y() << " " << v.z() << " " << rgb.x() << " " << rgb.y() << " " << rgb.z() << std::endl;
	}

	// Helper function to write single vertex to OBJ file
	static void write_vertex_to_xyz(std::ofstream& output, const Vector3& v)
	{
		output << v.x() << " " << v.y() << " " << v.z() << std::endl;
	}

	// Helper function to write face
	static void write_face(std::ofstream& output, const Vector3i& f)
	{
		output << "f " << f.x() << " " << f.y() << " " << f.z() << std::endl;
	}

	// Helper function to write face
	static void write_face(std::ofstream& output, const Eigen::Vector4i& f)
	{
		output << "f " << f.x() << " " << f.y() << " " << f.z() << " " << f.w() << std::endl;
	}

	// Helper function to write line
	static void write_line(std::ofstream& output, const Eigen::Vector2i& l)
	{
		output << "l " << l.x() << " " << l.y() << std::endl;
	}

	// Helper function to write full cube (using relative vertex positions in the OBJ file - support for this should be widespread by now)
	inline void writeCube(const Vector3& nodeOrigin, const Vector3& unit, const size_t& faceBegIdx, std::ofstream& output)
	{
		//	   2-------1
		//	  /|      /|
		//	 / |     / |
		//	7--|----8  |
		//	|  4----|--3
		//	| /     | /
		//	5-------6
		// Create vertices
		Vector3 v1 = nodeOrigin + Vector3(0, unit.y(), unit.z());
		Vector3 v2 = nodeOrigin + Vector3(0, 0, unit.z());
		Vector3 v3 = nodeOrigin + Vector3(0, unit.y(), 0);
		Vector3 v4 = nodeOrigin;
		Vector3 v5 = nodeOrigin + Vector3(unit.x(), 0, 0);
		Vector3 v6 = nodeOrigin + Vector3(unit.x(), unit.y(), 0);
		Vector3 v7 = nodeOrigin + Vector3(unit.x(), 0, unit.z());
		Vector3 v8 = nodeOrigin + Vector3(unit.x(), unit.y(), unit.z());

		// write them in reverse order, so relative position is -i for v_i
		write_vertex(output, v1);
		write_vertex(output, v2);
		write_vertex(output, v3);
		write_vertex(output, v4);
		write_vertex(output, v5);
		write_vertex(output, v6);
		write_vertex(output, v7);
		write_vertex(output, v8);

		// create faces
#if defined(MESH_WRITE)
		// back
		write_face(output, Eigen::Vector3i(faceBegIdx + 1, faceBegIdx + 3, faceBegIdx + 4));
		write_face(output, Eigen::Vector3i(faceBegIdx + 1, faceBegIdx + 4, faceBegIdx + 2));
		// bottom								   
		write_face(output, Eigen::Vector3i(faceBegIdx + 4, faceBegIdx + 3, faceBegIdx + 6));
		write_face(output, Eigen::Vector3i(faceBegIdx + 4, faceBegIdx + 6, faceBegIdx + 5));
		// right								   	 
		write_face(output, Eigen::Vector3i(faceBegIdx + 3, faceBegIdx + 1, faceBegIdx + 8));
		write_face(output, Eigen::Vector3i(faceBegIdx + 3, faceBegIdx + 8, faceBegIdx + 6));
		// top									   	 
		write_face(output, Eigen::Vector3i(faceBegIdx + 1, faceBegIdx + 2, faceBegIdx + 7));
		write_face(output, Eigen::Vector3i(faceBegIdx + 1, faceBegIdx + 7, faceBegIdx + 8));
		// left									   	 
		write_face(output, Eigen::Vector3i(faceBegIdx + 2, faceBegIdx + 4, faceBegIdx + 5));
		write_face(output, Eigen::Vector3i(faceBegIdx + 2, faceBegIdx + 5, faceBegIdx + 7));
		// front								   	  
		write_face(output, Eigen::Vector3i(faceBegIdx + 5, faceBegIdx + 6, faceBegIdx + 8));
		write_face(output, Eigen::Vector3i(faceBegIdx + 5, faceBegIdx + 8, faceBegIdx + 7));
#elif defined(CUBE_WRITE)
		// back
		write_face(output, Eigen::Vector4i(faceBegIdx + 3, faceBegIdx + 4, faceBegIdx + 2, faceBegIdx + 1));
		// bottom								   
		write_face(output, Eigen::Vector4i(faceBegIdx + 6, faceBegIdx + 5, faceBegIdx + 4, faceBegIdx + 3));
		// right								   	 
		write_face(output, Eigen::Vector4i(faceBegIdx + 1, faceBegIdx + 8, faceBegIdx + 6, faceBegIdx + 3));
		// top									   	 
		write_face(output, Eigen::Vector4i(faceBegIdx + 1, faceBegIdx + 2, faceBegIdx + 7, faceBegIdx + 8));
		// left									   	 
		write_face(output, Eigen::Vector4i(faceBegIdx + 4, faceBegIdx + 5, faceBegIdx + 7, faceBegIdx + 2));
		// front								   	  
		write_face(output, Eigen::Vector4i(faceBegIdx + 8, faceBegIdx + 7, faceBegIdx + 5, faceBegIdx + 6));
#else
		write_line(output, Eigen::Vector2i(faceBegIdx + 1, faceBegIdx + 2));
		write_line(output, Eigen::Vector2i(faceBegIdx + 2, faceBegIdx + 7));
		write_line(output, Eigen::Vector2i(faceBegIdx + 7, faceBegIdx + 8));
		write_line(output, Eigen::Vector2i(faceBegIdx + 8, faceBegIdx + 1));

		write_line(output, Eigen::Vector2i(faceBegIdx + 3, faceBegIdx + 4));
		write_line(output, Eigen::Vector2i(faceBegIdx + 4, faceBegIdx + 5));
		write_line(output, Eigen::Vector2i(faceBegIdx + 5, faceBegIdx + 6));
		write_line(output, Eigen::Vector2i(faceBegIdx + 6, faceBegIdx + 3));

		write_line(output, Eigen::Vector2i(faceBegIdx + 3, faceBegIdx + 1));
		write_line(output, Eigen::Vector2i(faceBegIdx + 4, faceBegIdx + 2));
		write_line(output, Eigen::Vector2i(faceBegIdx + 5, faceBegIdx + 7));
		write_line(output, Eigen::Vector2i(faceBegIdx + 6, faceBegIdx + 8));
#endif

		//faceBegIdx += 8;
	}

	inline void writePointCloud(const std::vector<Vector3>& points, std::ofstream& output)
	{
		for (size_t i = 0; i < points.size(); ++i)
			write_vertex(output, points[i]);
	}

	inline void writePointCloud_xyz(const std::vector<Vector3>& points, std::ofstream& output)
	{
		for (size_t i = 0; i < points.size(); ++i)
			write_vertex_to_xyz(output, points[i]);
	}

	inline void writePointCloud(const std::vector<Vector3>& points, const std::vector<Vector3>& rgbs, std::ofstream& output)
	{
		for (size_t i = 0; i < points.size(); ++i)
			write_vertex(output, points[i], rgbs[i]);
	}

	inline void writePointCloud(const MatrixX& points, std::ofstream& output)
	{
		for (size_t i = 0; i < points.size(); ++i)
			write_vertex(output, points.row(i));
	}

	inline void writePointCloud_xyz(const MatrixX& points, std::ofstream& output)
	{
		for (size_t i = 0; i < points.rows(); ++i)
			write_vertex_to_xyz(output, points.row(i));
	}

	inline void writePointCloud(const MatrixX& points, const std::vector<Vector3>& rgbs, std::ofstream& output)
	{
		for (size_t i = 0; i < points.rows(); ++i)
			write_vertex(output, points.row(i), rgbs[i]);
	}

	inline void writePointCloud(const Vector3& point, const Vector3& rgb, std::ofstream& output)
	{
		write_vertex(output, point, rgb);
	}
} // namespace geometry::gvis

namespace geometry
{
	template <typename Real>
	struct Edge
	{
		using type = std::pair<Real, Real>;
	};

	// An Axis Aligned Box (AAB) of a certain Real - to be initialized with a boxOrigin and boxEnd
	template <typename Real>
	struct AABox {
		//using Real = typename Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
		Real boxOrigin;
		Real boxEnd;
		Real boxWidth;

		_CUDA_GENERAL_CALL_ AABox() : boxOrigin(Real()), boxEnd(Real()), boxWidth(Real()) {}
		_CUDA_GENERAL_CALL_ AABox(const Real& _boxOrigin, const Real& _boxEnd) : boxOrigin(_boxOrigin), boxEnd(_boxEnd), boxWidth(_boxEnd - _boxOrigin) {}

		_CUDA_GENERAL_CALL_ void scaleAndTranslate(const double& scale_factor, const Vector3& translation)
		{
			// Calculate centroid point
			const Real center = (boxOrigin + boxEnd) / 2.0;

			// Zoom and translation based on the centroid point
			const Real scaled_min_point = (boxOrigin - center) * scale_factor + center + translation;
			const Real scaled_max_point = (boxEnd - center) * scale_factor + center + translation;

			// Update coordinate
			boxOrigin = scaled_min_point;
			boxEnd = scaled_max_point;
		}

		_CUDA_GENERAL_CALL_ AABox<Real>(const AABox<Real>& _box)
		{
			boxOrigin = _box.boxOrigin;
			boxEnd = _box.boxEnd;
			boxWidth = _box.boxWidth;
		}

		_CUDA_GENERAL_CALL_ AABox<Real>& operator=(const AABox<Real>& _box)
		{
			boxOrigin = _box.boxOrigin;
			boxEnd = _box.boxEnd;
			boxWidth = _box.boxWidth;

			return *this;
		}

		bool isInBox(const Real& queryPoint)
		{
			int dim = boxOrigin.rows();
			for (int i = 0; i < dim; ++i)
				if (!(boxOrigin(i) <= queryPoint(i) && queryPoint(i) <= boxEnd(i))) return false;
			return true;
		}

		void output(const std::string& filename)
		{
			str_util::checkDir(filename);
			std::ofstream out(filename);
			if (!out) { LOG::qpError("I/O: File ", filename.c_str(), " could not be opened!"); return; }
			LOG::qpInfo("Output bounding-box to ", std::quoted(filename), " ...");

			gvis::writeCube(boxOrigin, boxWidth, 0, out);

			out.close();
		}
	};

	template <typename Real>
	struct Triangle
	{
		Real p1, p2, p3;
		Real normal;
		double area;
		double dir;

		_CUDA_GENERAL_CALL_ Triangle() {}

		_CUDA_GENERAL_CALL_ Triangle(const Real& _p1, const Real& _p2, const Real& _p3) :p1(_p1), p2(_p2), p3(_p3) {
			normal = (p2 - p1).cross(p3 - p1);
		}
	};

	template <typename Real = Scalar, class T = Vector<Real, 3>>
	struct Tet_
	{
		typedef T Point;
		typedef std::array<Point, 4> Tetrahedron;
		typedef std::array<int, 4>   TetNeighbor;
	};
} // namespace geometry

NAMESPACE_END(mscut)