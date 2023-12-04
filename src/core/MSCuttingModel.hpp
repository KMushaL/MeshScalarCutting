#pragma once
#include "geometry/PolyMesh.hpp"
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

	public:
		/* struct type */
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
			Vector3 pos;						// 采样点坐标
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
		/* Constane Variables */
		constexpr static int WEIGHT_MAX_ITER = 20;			// 计算权重时的最大迭代次数

		constexpr static double SINGULAR_DIS_EPSILON = 0.0; // 与不可导点的距离值

		constexpr static double NUMERICAL_EPSILON = 1e-9;   // 控制后处理中的停止迭代误差(包括梯度接近于0、值接近于c等)

		constexpr static double PROJECT_EPSILON = 1e-6;     // 控制后处理中的投影误差

	private:
		/* Data */
		/* 标量函数 */
		ScalarFunc scalarFunc;

		/* 采样相关 */
		int numSamplesPerEdge;  // 每条边固定的采样数量 TODO：后期不可能固定

		std::vector<CELL_TYPE> meshVertsType;

		std::vector<SamplePoint> samplePoints; // 整个模型所有边上的采样点(包括边的端点)

		std::vector<SampleFacet> sampleFacets; // 包含了边采样点的模型面

		std::vector<Vector3> singulars;

		/* 后处理相关 */
		int POST_PROCESSING_MAX_ITER;

	private:
		/* Adaptors */
		PowerDiagramAdaptor pdAdaptor;

		CDTAdaptor cdtAdaptor;

	private:
		/* 用于Cutting的数据结构 */
		using MeshComponent = geometry::MeshComponent<Scalar, Vector3>;

		// 用于控制内外mesh的输出点索引，注意是从0开始
		int insideMeshVertIdx = 0;
		int outsideMeshVertIdx = 0;

		std::vector<MeshComponent> insideMesh;
		std::vector<MeshComponent> outsideMesh;

	public:
		/* Constructor and Destructor */
		MSCuttingModel() noexcept = default;

		/* 非MeshNorm版构造函数 */
		MSCuttingModel(const std::string& filename,
			int _numSamples,
			const ScalarFunc& _scalarFunc,
			const std::vector<Vector3>& _singulars) noexcept :
			PolyMesh(filename), numSamplesPerEdge(_numSamples), scalarFunc(_scalarFunc), singulars(_singulars)
		{
			samplePoints.reserve(_numSamples * numMeshEdges);
			sampleFacets.resize(numMeshFaces, SampleFacet());

			meshVertsType.resize(numMeshVerts);
			pdAdaptor = PowerDiagramAdaptor();
		}

		/* MeshNorm版构造函数 */
		MSCuttingModel(const std::string& filename,
			int _numSamples,
			const ScalarFunc& _scalarFunc,
			const std::vector<Vector3>& _singulars,
			const std::string& norm_out_file,
			double _scaleFactor = 1.0) noexcept :
			PolyMesh(filename, true, _scaleFactor, norm_out_file),
			numSamplesPerEdge(_numSamples), scalarFunc(_scalarFunc), singulars(_singulars)
		{
			samplePoints.reserve(numSamplesPerEdge * numMeshEdges);
			sampleFacets.resize(numMeshFaces, SampleFacet());

			meshVertsType.resize(numMeshVerts);
			pdAdaptor = PowerDiagramAdaptor();
		}

		~MSCuttingModel() noexcept = default;

	private:
		/* Details */

		/* 设置单个面的采样点二维局部坐标和权重 */
		int updateFacetLocalCoord(int facetIdx, Polygon_2& facetBoundary);

		/* Transform coordinate of point from local to global after computing isolines */
		Eigen::Vector3d getGlobalCoordInFacet(int facetIdx, const PowerDiagramPoint_2& point_2);

		/* 计算单个面的power diagram */
		int computePDForFacet(int facetIdx,
			int& globalOutVertIdx, // 用于控制等值线的输出索引值
			std::vector<PowerDiagramPoint_3>& apolloniusDiagramPoints,
			std::vector<std::pair<int, int>>& apolloniusDiagramLines,
			std::unordered_map<int, std::vector<int>>& edgeTable,
			std::map<PowerDiagramPoint_3, int>&);

		/* 对单个面的等值线进行后处理 */
		std::vector<PowerDiagramPoint_3> postProcessFacetPoints(int faceIdx,
			const std::vector<std::vector<PowerDiagramPoint_3>>&,
			const std::map<PowerDiagramPoint_3, int>&,
			std::array<std::vector<PowerDiagramPoint_3>, 3>&,
			std::map<PowerDiagramPoint_3, int>&,
			std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<PowerDiagramPoint_3, int>>>>&);

		/* 计算单个面的CDT */
		using CellTriangle = std::vector<CDTriangle>;
		std::vector<CellTriangle> computeCDTForFacet(int faceIdx,
			const std::vector<PowerDiagramPoint_3>&,
			const std::unordered_map<int, std::vector<int>>&,
			const std::array<std::vector<PowerDiagramPoint_3>, 3>&,
			const std::map<PowerDiagramPoint_3, int>&);

	private:
		/* 对mesh的所有边进行采样 */
		int samplePointPerEdge();

		/* 计算所有面的power diagram，并输出至out，以及返回mae */
		double computeIsoline(std::ofstream& out);

		double computeMAE();

	private:
		/* Visualization */
		void outputSamplePoints(std::ofstream& out);

		void outputCutMesh(std::ofstream& out_1, std::ofstream& out_2);

	public:
		/* API for user */
		/*
			第一个参数控制后处理中的迭代次数
			后面三个参数分别传递的等值线的保存位置以及切割后内外mesh的保存位置
		*/
		bool launch(int iter,
			std::ofstream& mae_out,
			const std::string& isolineVisFile,
			const std::string& insideMeshVisFile = "",
			const std::string& outsideMeshVisFile = "");

#if ENABLE_TEST
	public:
		/* Test APIs for us */
		bool testSamplingPoints(std::ofstream& out_1, std::ofstream& out_2);

		bool testLocalGlobalTransform(int facetIdx);

		bool testComputeADForFacet(int facetIdx, std::ofstream& out);
#endif


		/* [[DEPRECATED]] */
		/*
		int m_numSamplesPerEdge;
		[[deprecated]] MSCuttingModel(const std::string& filename,
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
		}*/

		/* [[DEPRECATED]] */
		/* TODO: 目前合理的条件是相邻点对的隐函数值异号
		std::vector<SamplePoint> validSamplePoints; // 整个模型所有边上那些合理的采样点(可能包括边的端点)

		void outputSamplePoints(std::ofstream& out_1, std::ofstream& out_2);
		*/
	};

} // namespace core
NAMESPACE_END(mscut)