﻿#pragma once
#include "geometry/PolyMesh.hpp"
#include "core/ApolloniusGraphAdaptor.hpp"
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
		/*
			因为我们只有在调用Apollonius时用到了CGAL，为了方便起见，
			三维下的点仍用Eigen的Data type，二维下的点使用CGAL的Point_2
		*/
		using ApolloniusDiagramPoint_2 = typename ApolloniusGraphAdaptor::ApolloniusDiagramPoint_2;
		using ApolloniusDiagramLine_2 = typename ApolloniusGraphAdaptor::ApolloniusDiagramLine_2;
		using ApolloniusDiagramPoint_3 = Eigen::Vector3d;
		using ApolloniusDiagramLine_3 = std::pair<ApolloniusDiagramPoint_3, ApolloniusDiagramPoint_3>;

		/* 采样点定义 */
		struct SamplePoint {
			int mEdgeIdx = -1;								// 属于的边(注意是model edge，而不是half edge)
			Eigen::Vector3d pos;							// 采样点坐标

			SamplePoint() {}
			SamplePoint(const Eigen::Vector3d& _pos, int _mEdgeIdx) :pos(_pos), mEdgeIdx(_mEdgeIdx) {}
		};

		/* (二维)带权重点(站点)定义 */
		struct WeightPoint {
			using _Point = ApolloniusDiagramPoint_2;
			_Point pos;
			double weight = 0;								// 先都设为0

			WeightPoint() {}
			WeightPoint(const _Point& _pos, double _weight) :pos(_pos), weight(_weight) {}
		};

		/* 包含了采样点、站点(权重点)的面定义 */
		struct SampleFacet {
			using Basis = typename Eigen::Vector3d;

			Basis x_basis, y_basis;							// 局部坐标系
			Eigen::Vector3d origin;							// 局部坐标系的原点
			std::vector<SamplePoint> aroundSamplePoints;	// 包含的(边)采样点
			std::vector<WeightPoint> sites;					// 所有站点

			SampleFacet() {}
			SampleFacet(const std::vector<SamplePoint>& _aroundSamplePoints, const std::vector<WeightPoint>& _sites)
				:aroundSamplePoints(_aroundSamplePoints), sites(_sites) {}
		};

	private:
		/* Data */
		int numSamplesPerEdge;  // 每条边固定的采样数量 TODO：后期不可能固定

		std::vector<SamplePoint> samplePoints; // 整个模型所有边上的采样点
		std::vector<SampleFacet> sampleFacets; // 包含了边采样点的模型面

		ApolloniusGraphAdaptor agAdaptor;

	public:
		/* Constructor and Destructor */
		MSCuttingModel() noexcept = default;

		MSCuttingModel(const std::string& filename, int _numSamples) noexcept :
			PolyMesh(filename), numSamplesPerEdge(_numSamples) {
			samplePoints.reserve(_numSamples * numMeshEdges);
			sampleFacets.resize(numMeshFaces, SampleFacet());

			agAdaptor = ApolloniusGraphAdaptor();
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
		Eigen::Vector3d getGlobalCoordInFacet(int facetIdx, const ApolloniusDiagramPoint_2& point_2);

		int computeApolloniusGraphForFacet(int facetIdx,
			int& globalOutVertIdx,
			std::vector<ApolloniusDiagramPoint_3>& apolloniusDiagramPoints,
			std::vector<std::pair<int, int>>& apolloniusDiagramLines);

	private:
		/* Methods for Our Algorithm */
		int samplePointPerEdge();

		void computeApolloniusGraph(std::ofstream& out);

	private:
		/* Visualization */
		void outputSamplePoints(std::ofstream& out);

	public:
		/* API for user */
		bool launch(const std::string& ad_vis_file); // 传递的ad_vis_file是最后计算结果保存的文件位置

	public:
		/* Test APIs for us */
		bool testSamplingPoints(std::ofstream& out);

		bool testLocalGlobalTransform(int facetIdx);

		bool testComputeADForFacet(int facetIdx, std::ofstream& out);
	};

} // namespace core
NAMESPACE_END(mscut)