#pragma once
#include "geometry/PolyMesh.hpp"
#include "core/ApolloniusGraphAdaptor.hpp"
#include "core/PowerDiagramAdaptor.hpp"
#include "core/CDTAdaptor.hpp"
#include <functional>
#include <unordered_set>

NAMESPACE_BEGIN(mscut)
namespace core
{
	using namespace detail;
	using namespace geometry;


	/* 为了你方便以后调试，先不用模板 */
	//template<bool isNormalized = false>
	class MSCuttingModel :public PolyMesh
	{
	public:
		/* Member type */
		using Point_2 = typename PowerDiagramAdaptor::Point_2;
		using Point_3 = typename PowerDiagramAdaptor::Point_3;
		using Segment_2 = typename PowerDiagramAdaptor::Segment_2;
		using Segment_3 = typename PowerDiagramAdaptor::Segment_3;
		using Polygon_2 = typename PowerDiagramAdaptor::Polygon_2;
		/*
			因为我们只有在调用Apollonius时用到了CGAL，为了方便起见，
			三维下的点仍用Eigen的Data type，二维下的点使用CGAL的Point_2
		*/
		using PowerDiagramPoint_2 = typename PowerDiagramAdaptor::PowerDiagramPoint_2;
		using PowerDiagramLine_2 = typename PowerDiagramAdaptor::PowerDiagramLine_2;
		using PowerDiagramPoint_3 = Eigen::Vector3d;
		using PowerDiagramLine_3 = std::pair<PowerDiagramPoint_3, PowerDiagramPoint_3>;
		using CDTriangle = typename CDTAdaptor::CDTriangle;

		using ScalarValFunc = std::function<Scalar(const PowerDiagramPoint_3&)>;
		using GradFunc = std::function<Eigen::Vector3d(const PowerDiagramPoint_3&)>;

		struct ScalarFunc {
			ScalarValFunc val;
			GradFunc grad;
		};

		enum class CELL_TYPE {
			INSIDE,
			OUTSIDE,
			EQUAL
		};

		/* 采样点定义 */
		struct SamplePoint {
			int mEdgeIdx = -1;					// 属于的边(注意是model edge，而不是half edge)
			Vector3 pos;							// 采样点坐标
			double val;

			SamplePoint() {}
			SamplePoint(const Vector3& _pos, int _mEdgeIdx, double _val) :pos(_pos), mEdgeIdx(_mEdgeIdx), val(_val) {}

			SamplePoint(const SamplePoint& other) { mEdgeIdx = other.mEdgeIdx; pos = other.pos; val = other.val; }
		};

		/* (二维)带权重点(站点)定义 */
		struct WeightPoint {
			using _Point = PowerDiagramPoint_2;
			_Point pos;
			double weight = 0;								// 先都设为0
			double f_val = DINF;

			WeightPoint() {}
			WeightPoint(const _Point& _pos, double _weight) :pos(_pos), weight(_weight) {}
		};

		/* 包含了采样点、站点(权重点)的面定义 */
		struct SampleFacet {
			using Basis = typename Eigen::Vector3d;

			Basis x_basis, y_basis;							// 局部坐标系
			Eigen::Vector3d origin;							// 局部坐标系的原点
			std::set<Point3> aroundSamplePointsSet;
			std::vector<SamplePoint> aroundSamplePoints;	// 包含的(边)采样点
			std::vector<WeightPoint> sites;					// 所有站点

			SampleFacet() {}
			SampleFacet(const std::vector<SamplePoint>& _aroundSamplePoints, const std::vector<WeightPoint>& _sites)
				:aroundSamplePoints(_aroundSamplePoints), sites(_sites) {}
		};

		/* 计算CDT时所传入的constrained cell定义 */
		struct CDT_Cell {
			using _Point = Point_3;
			using _Edge = Segment_3;

			std::vector<_Point> bdPoints;
			std::vector<_Edge> bdEdges;

			CDT_Cell() noexcept = default;
			CDT_Cell(const std::vector<_Edge>& _bdEdges) :bdEdges(_bdEdges) {}
		};

	private:
		/* Data */
		int m_numSamplesPerEdge;
		int numSamplesPerEdge;  // 每条边固定的采样数量 TODO：后期不可能固定

		std::vector<CELL_TYPE> meshVertsType;

		std::vector<Vector3> singulars;
		double singularEpsilon = 0.5;

		std::vector<SamplePoint> samplePoints; // 整个模型所有边上的采样点(包括边的端点)

		// TODO: 目前合理的条件是相邻点对的隐函数值异号
		std::vector<SamplePoint> validSamplePoints; // 整个模型所有边上那些合理的采样点(可能包括边的端点)

		std::vector<SampleFacet> sampleFacets; // 包含了边采样点的模型面

		ApolloniusGraphAdaptor agAdaptor;

		PowerDiagramAdaptor pdAdaptor;

		CDTAdaptor cdtAdaptor;

		ScalarFunc scalarFunc; // 标量函数

		using MeshComponent = geometry::MeshComponent<Scalar, Vector3>;
		std::vector<MeshComponent> insideMesh;
		std::vector<MeshComponent> outsideMesh;
		int insideMeshVertIdx = 0, outsideMeshVertIdx = 0; // 用于控制内外mesh的输出点索引，注意是从0开始

	public:
		/* Constructor and Destructor */
		MSCuttingModel() noexcept = default;

		MSCuttingModel(const std::string& filename, const std::string& norm_out_file, const ScalarFunc& _scalarFunc, int _numSamples) noexcept :
			PolyMesh(filename, norm_out_file), scalarFunc(_scalarFunc), numSamplesPerEdge(_numSamples) {
			samplePoints.reserve(numSamplesPerEdge * numMeshEdges);
			sampleFacets.resize(numMeshFaces, SampleFacet());

			meshVertsType.resize(numMeshVerts);

			agAdaptor = ApolloniusGraphAdaptor();
			pdAdaptor = PowerDiagramAdaptor();
		}

		MSCuttingModel(const std::string& filename,
			const ScalarFunc& _scalarFunc, const std::vector<Vector3>& _singulars,
			int _m_numSamples, int _n_numSamples) noexcept :
			PolyMesh(filename), scalarFunc(_scalarFunc), singulars(_singulars),
			m_numSamplesPerEdge(_m_numSamples), numSamplesPerEdge(_n_numSamples)
		{
			samplePoints.reserve(numSamplesPerEdge * numMeshEdges);
			sampleFacets.resize(numMeshFaces, SampleFacet());

			meshVertsType.resize(numMeshVerts);

			agAdaptor = ApolloniusGraphAdaptor();
			pdAdaptor = PowerDiagramAdaptor();
		}

		MSCuttingModel(const std::string& filename, int _numSamples) noexcept :
			PolyMesh(filename), numSamplesPerEdge(_numSamples) {
			samplePoints.reserve(_numSamples * numMeshEdges);
			sampleFacets.resize(numMeshFaces, SampleFacet());

			meshVertsType.resize(numMeshVerts);

			agAdaptor = ApolloniusGraphAdaptor();
			pdAdaptor = PowerDiagramAdaptor();
		}

		MSCuttingModel(const std::string& filename) noexcept : MSCuttingModel(filename, 0) {}

		~MSCuttingModel() noexcept = default;

	private:
		/* Details */
		/*
		   Compute local coordinate(while updating related data) and get the boundary
		   for sample facet before computing Apollonius Graph
		*/
		int updateFacetLocalCoord(int facetIdx, Polygon_2& facetBoundary);

		/* Transform coordinate of point from local to global after computing Apollonius Graph */
		Eigen::Vector3d getGlobalCoordInFacet(int facetIdx, const PowerDiagramPoint_2& point_2);

		int computePDForFacet(int facetIdx,
			int& globalOutVertIdx, // 用于控制等值线的输出索引值
			std::vector<PowerDiagramPoint_3>& apolloniusDiagramPoints,
			std::vector<std::pair<int, int>>& apolloniusDiagramLines,
			std::unordered_map<int, std::vector<int>>& edgeTable,
			std::map<PowerDiagramPoint_3, int>&);

		static constexpr double PROJ_EPSILON = 1e-6;
		std::vector<PowerDiagramPoint_3> postProcessFacetPoints(int faceIdx,
			const std::vector<PowerDiagramPoint_3>&,
			const std::map<PowerDiagramPoint_3, int>&,
			std::array<std::vector<PowerDiagramPoint_3>, 3>&,
			std::map<PowerDiagramPoint_3, int>&,
			std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<PowerDiagramPoint_3, int>>>>&);

		using CellTriangle = std::vector<CDTriangle>;
		std::vector<CellTriangle> computeCDTForFacet(int faceIdx,
			const std::vector<PowerDiagramPoint_3>&,
			const std::unordered_map<int, std::vector<int>>&,
			const std::array<std::vector<PowerDiagramPoint_3>, 3>&,
			const std::map<PowerDiagramPoint_3, int>&);

	private:
		/* Methods for Our Algorithm */
		int samplePointPerEdge();

		void computePowerDiagram(std::ofstream& out);

	private:
		/* Visualization */
		void outputSamplePoints(std::ofstream& out);

		void outputSamplePoints(std::ofstream& out_1, std::ofstream& out_2);

		void outputCutMesh(std::ofstream& out_1, std::ofstream& out_2);

	public:
		/* API for user */
		bool launch(const std::string& ad_vis_file,
			const std::string& insideMeshVisFile = "",
			const std::string& outsideMeshVisFile = ""); // 分别传递的等值线的保存位置以及切割后的内外mesh保存位置

	public:
		/* Test APIs for us */
		bool testSamplingPoints(std::ofstream& out_1, std::ofstream& out_2);

		bool testLocalGlobalTransform(int facetIdx);

		bool testComputeADForFacet(int facetIdx, std::ofstream& out);
	};

} // namespace core
NAMESPACE_END(mscut)