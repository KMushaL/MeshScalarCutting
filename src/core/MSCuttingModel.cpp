#pragma once
#include "MSCuttingModel.hpp"
#include "detail/Handler.hpp"
#include <map>
#include <fstream>
#include <igl/AABB.h>

NAMESPACE_BEGIN(mscut)
namespace core
{
	/* Details */
	/**
	* @param facetIdx: 待转换面的索引
	* @param facetBoundary: 待转换面边界(在函数体内被更新)
	* @return: 站点数量
	*/
	int MSCuttingModel::updateFacetLocalCoord(int facetIdx, Polygon_2& facetBoundary)
	{
		using Basis = typename SampleFacet::Basis;

		SampleFacet& curSampleFacet = sampleFacets[facetIdx];
		const auto facetNormal = this->getFaceVec()[facetIdx]->normal;
		const auto& aroundEndVerts = this->getFaceVec()[facetIdx]->aroundVerts();
		const auto& aroundSamplePoints = curSampleFacet.aroundSamplePoints;

		std::vector<Point_2> facet_dim2;

		// 注意facet不一定是三角形，但是没关系，我们仍然可以将 v1 当作局部坐标系的原点
		// 并且仍可通过 v1--v2 和 v1--v3 以及法向量确定一个局部坐标系
		Eigen::Vector3d v1 = aroundEndVerts[0]->pos, v2 = aroundEndVerts[1]->pos, v3 = aroundEndVerts[2]->pos;

		// 原点坐标
		curSampleFacet.origin = v1;
		// x轴和y轴
		Basis x_basis = (v2 - v1).normalized();
		Basis y_basis = facetNormal.cross(x_basis).normalized();
		curSampleFacet.x_basis = x_basis;
		curSampleFacet.y_basis = y_basis;

		// 计算所有端点的局部坐标(TODO: 目前将端点也当作站点的做法是否合理？)
		for (int i = 0; i < aroundEndVerts.size(); ++i)
		{
			WeightPoint site;

			auto endVert_dim3 = aroundEndVerts[i]->pos;
			Point_2 endVert_dim2 = Point_2((endVert_dim3 - v1).dot(x_basis), (endVert_dim3 - v1).dot(y_basis));

			site.pos = endVert_dim2;

			meshVertsType[aroundEndVerts[i]->index] = (VERT_TYPE)(scalarFunc.val(endVert_dim3) > 0);

			/*site.weight = -scalarFunc.val(endVert_dim3) / scalarFunc.grad(endVert_dim3).norm();
			curSampleFacet.sites.emplace_back(site);*/

			// 更新边界信息
			facetBoundary.push_back(site.pos);
		}

		// 计算所有采样点的局部坐标
		for (int i = 0; i < aroundSamplePoints.size(); ++i)
		{
			WeightPoint site;

			auto sampleVert_dim3 = aroundSamplePoints[i].pos;
			Point_2 sampleVert_dim2 = Point_2((sampleVert_dim3 - v1).dot(x_basis), (sampleVert_dim3 - v1).dot(y_basis));

			site.pos = sampleVert_dim2;

			Scalar f_val = scalarFunc.val(sampleVert_dim3);
			Vector3 f_grad = scalarFunc.grad(sampleVert_dim3);

			site.f_val = f_val;
			/*site.weight = abs(-f_val / f_grad.norm()) * 0.9*/;
			double _weight = f_val / f_grad.norm();
			site.weight = _weight * _weight * 1.1;

			//std::cout << "pos: " << sampleVert_dim3.transpose() << ", weight = " << site.weight << std::endl;

			curSampleFacet.sites.emplace_back(site);
		}

		return curSampleFacet.sites.size();
	}

	/**
	* @param facetIdx: 待转换的局部坐标点所在的面的索引
	* @param point_2: 待转换的局部坐标点
	* @return: 转换后的全局坐标
	*/
	Eigen::Vector3d MSCuttingModel::getGlobalCoordInFacet(int facetIdx, const ApolloniusDiagramPoint_2& point_2)
	{
		// TODO: 检查point_2是否在facetIdx对应的平面上

		using Basis = typename SampleFacet::Basis;
		const SampleFacet& curSampleFacet = sampleFacets[facetIdx];

		// 获取原点和坐标系
		Eigen::Vector3d origin = curSampleFacet.origin;
		Basis x_basis = curSampleFacet.x_basis;
		Basis y_basis = curSampleFacet.y_basis;

		Eigen::Vector3d point_3;
		// 计算全局坐标
		point_3 = origin + point_2.x() * x_basis + point_2.y() * y_basis;
		// 避开小数值的扰动
		point_3 = (point_3.array().abs() > ZERO_EPSILON).select(point_3, 0);

		return point_3;
	}

	/**
	* brief: 为facetIdx这个面计算Apollonius Diagram
	* @param facetIdx: 面的索引
	* @return: Apollonius Diagram边的数量
	*/
	int MSCuttingModel::computeApolloniusGraphForFacet(int facetIdx,
		int& globalOutVertIdx,
		std::vector<ApolloniusDiagramPoint_3>& apolloniusDiagramPoints,
		std::vector<std::pair<int, int>>& apolloniusDiagramLines,
		std::unordered_map<int, std::vector<int>>& edgeTable)
	{
		const int originGlobalOutVertIdx = globalOutVertIdx;

		// 从local point到输出下标的映射，主要用于防止输出重复的点
		// 为什么不弄个全局的(也就是对所有面的)？为了后续稍微好并行(对每个面同时计算Apollonius Diagram)一点
		std::map<Point_2, int> pointToGlobalOutIdx;

		// 计算局部坐标和边界
		Polygon_2 facetBoundary;
		updateFacetLocalCoord(facetIdx, facetBoundary);

		// 计算Apollonius Diagram
		//auto adLineSegmentList = agAdaptor.computeAGForBoundary(sampleFacets[facetIdx].sites, facetBoundary);
		auto adLineSegmentList = pdAdaptor.computePDForBoundary(sampleFacets[facetIdx].sites, facetBoundary);

		// 遍历每条边，转回全局坐标
		for (const auto& lineSegemnt : adLineSegmentList)
		{
			const ApolloniusDiagramPoint_2 vert_1 = lineSegemnt.vertex(0);
			const ApolloniusDiagramPoint_2 vert_2 = lineSegemnt.vertex(1);

			int outVertIdx_1, outVertIdx_2; // 当前这两个点应输出的下标
			if (pointToGlobalOutIdx.count(vert_1)) outVertIdx_1 = pointToGlobalOutIdx.at(vert_1);
			else
			{
				ApolloniusDiagramPoint_3 apolloniusDiagramPoint = getGlobalCoordInFacet(facetIdx, vert_1);
				apolloniusDiagramPoints.emplace_back(apolloniusDiagramPoint);

				outVertIdx_1 = globalOutVertIdx++;
				pointToGlobalOutIdx.insert(std::make_pair(vert_1, outVertIdx_1));
			}

			if (pointToGlobalOutIdx.count(vert_2)) outVertIdx_2 = pointToGlobalOutIdx.at(vert_2);
			else
			{
				ApolloniusDiagramPoint_3 apolloniusDiagramPoint = getGlobalCoordInFacet(facetIdx, vert_2);
				apolloniusDiagramPoints.emplace_back(apolloniusDiagramPoint);

				outVertIdx_2 = globalOutVertIdx++;
				pointToGlobalOutIdx.insert(std::make_pair(vert_2, outVertIdx_2));
			}

			// 保存这条边的两个顶点的索引
			// TODO: 注意，我们这里假定每个点当且仅当只会和一个点相邻
			apolloniusDiagramLines.emplace_back(std::make_pair(outVertIdx_1, outVertIdx_2));

			// edgeTable保存的是局部索引
			edgeTable[outVertIdx_1 - originGlobalOutVertIdx].emplace_back(outVertIdx_2 - originGlobalOutVertIdx);
			edgeTable[outVertIdx_2 - originGlobalOutVertIdx].emplace_back(outVertIdx_1 - originGlobalOutVertIdx);
		}

		return apolloniusDiagramLines.size();
	}

	/* Post-processing */
	std::vector<MSCuttingModel::ApolloniusDiagramPoint_3>
		MSCuttingModel::postProcessFacetPoints(const int faceIdx,
			const std::vector<ApolloniusDiagramPoint_3>& curFacetADPoints,
			std::array<std::vector<ApolloniusDiagramPoint_3>, 3>& edgePoints,
			std::map<ApolloniusDiagramPoint_3, std::pair<int, int>>& edgePointToIdx)
	{
		int numADPoints = curFacetADPoints.size();
		std::vector<ApolloniusDiagramPoint_3> resPoints;

		// aabb initialization
		igl::AABB<MatrixX, 3> tri_aabb;
		MatrixX singleTriV(3, 3);
		const auto& triVerts = this->getFaceVec()[faceIdx]->aroundVerts();
		singleTriV << triVerts[0]->pos.transpose(),
			triVerts[1]->pos.transpose(),
			triVerts[2]->pos.transpose();
		MatrixXi singleTriF(1, 3);
		singleTriF << 0, 1, 2;
		tri_aabb.init(singleTriV, singleTriF);

		using VertIdx = std::pair<int, int>;
		const std::array<VertIdx, 3> triEdges = {
			VertIdx(0, 1),
			VertIdx(1, 2),
			VertIdx(2, 0)
		};
		const std::array<Vector3, 3> triEdgesDir = {
			triVerts[1]->pos - triVerts[0]->pos,
			triVerts[2]->pos - triVerts[1]->pos,
			triVerts[0]->pos - triVerts[2]->pos
		};

#pragma omp parallel for
		for (int i = 0; i < numADPoints; ++i) {
			const ApolloniusDiagramPoint_3 originalADPoint = curFacetADPoints[i];

			// 判断点是否在边上
			int onEdgeIdx = -1;
			for (int i = 0; i < 3; ++i) {
				const Vector3 v1 = triVerts[triEdges[i].first]->pos;

				if (onEdgeIdx == -1 && triEdgesDir[i].cross(originalADPoint - v1).norm() < 1e-9) { onEdgeIdx = i; break; }
			}

			VectorX tri_sqrD; VectorXi tri_I; MatrixX tri_C; MatrixX projQ(1, 3);

			// for computing projection point
			ApolloniusDiagramPoint_3 lastPoint = originalADPoint;
			ApolloniusDiagramPoint_3 projPoint;
			Scalar lastVal; Vector3 lastGrad; double alpha;

			int iter = 0;
			while (iter < 50) {
				++iter;

				lastVal = scalarFunc.val(lastPoint);
				lastGrad = scalarFunc.grad(lastPoint);

				alpha = -lastVal / (lastGrad.squaredNorm() + 1e-6);

				// project to triangle
				if (onEdgeIdx == -1)
				{
					projPoint = alpha * lastGrad + lastPoint;

					projQ.row(0) = projPoint;
					tri_aabb.squared_distance(singleTriV, singleTriF, projQ, tri_sqrD, tri_I, tri_C);
					projPoint = tri_C.row(0);
				}
				else
				{
					const Vector3 proj_grad = triEdgesDir[onEdgeIdx] * (alpha * lastGrad.dot(triEdgesDir[onEdgeIdx]) / (triEdgesDir[onEdgeIdx].norm()));
					projPoint = proj_grad + lastPoint;
				}

				if ((projPoint - lastPoint).norm() < PROJ_EPSILON) break;
				else lastPoint = projPoint;
			};

			//resPoints[i] = projPoint;

			// TODO: 将去重工作移到外面，提升并行效率
#pragma omp critical
			{
				// 目前对后处理后的所有点做了一个去重工作
				if (edgePointToIdx.find(projPoint) == edgePointToIdx.end())
				{
					edgePointToIdx[projPoint] = std::make_pair(onEdgeIdx, i); // 得到该边界点在哪条边上的信息以及在edgeTable中的索引
					resPoints.emplace_back(projPoint);
					if (onEdgeIdx != -1) edgePoints[onEdgeIdx].emplace_back(projPoint);
				}

				if (i <= 20) LOG::qpInfo("#i = ", i, ", went ", iter, " iterations.");
				else LOG::qpWarn("#i = ", i, ", went ", iter, " iterations.");
			}
		}

		// 将边上的点统一按逆时针排序
		struct CCWSort {
			ApolloniusDiagramPoint_3 edgeBaseVert;
			Vector3 edgeDir;
			CCWSort(const ApolloniusDiagramPoint_3& _edgeBaseVert, const Vector3& _edgeDir) :edgeBaseVert(_edgeBaseVert), edgeDir(_edgeDir) {}

			bool operator()(const ApolloniusDiagramPoint_3& lhs, const ApolloniusDiagramPoint_3& rhs)
			{
				return ((lhs - edgeBaseVert).dot(edgeDir) < (rhs - edgeBaseVert).dot(edgeDir));
			}
		};
		donut::Loop<int, 3>([&](int i) {
			if (!edgePoints[i].empty()) {
				CCWSort ccwSort(triVerts[i]->pos, triEdgesDir[i]);
				std::sort(edgePoints[i].begin(), edgePoints[i].end(), ccwSort);
			}
			});

		return resPoints;
	}

	std::vector<typename MSCuttingModel::CellTriangle>
		MSCuttingModel::computeCDTForFacet(const int faceIdx,
			const std::vector<ApolloniusDiagramPoint_3>& isoLineVerts,
			//const std::vector<std::pair<int, int>>& isoLineVertIdx,
			const std::unordered_map<int, std::vector<int>>& edgeTable,
			const std::array<std::vector<ApolloniusDiagramPoint_3>, 3>& isoEdgePoints,
			const std::map<ApolloniusDiagramPoint_3, std::pair<int, int>>& isoEdgePointToIdx)
	{
		// 1. 搜寻边界边并确定区域
		// 1.1 按逆时针将所有边上的点(包括面的端点)按逆时针排序
		//std::queue<int> edgeCCWPointsIdxQueue;

		using VertIdxPair = std::pair<int, int>;
		std::queue<VertIdxPair> bdSegQueue;
		std::set<VertIdxPair> isVisSeg;

		//std::vector<ApolloniusDiagramPoint_3> edgeCCWPoints;

		// Initialize cycle points
		struct EdgePointWithType {
			ApolloniusDiagramPoint_3 pos;
			VERT_TYPE type = VERT_TYPE::EQUAL;

			EdgePointWithType(const ApolloniusDiagramPoint_3& _pos) :pos(_pos) {}
			EdgePointWithType(const ApolloniusDiagramPoint_3& _pos, const VERT_TYPE& _type) :pos(_pos), type(_type) {}
		};
		std::vector<EdgePointWithType> edgePWTVec;

		// TODO: 目前没有对edgeCyclePoints做一个去重工作，也就是说有可能三角形端点和等值线上的边界点几乎是一个点
		const auto& triVerts = this->getFaceVec()[faceIdx]->aroundVerts();
		edgePWTVec.emplace_back(EdgePointWithType(triVerts[0]->pos, meshVertsType[triVerts[0]->index]));
		for (int i = 0; i < isoEdgePoints[0].size(); ++i) edgePWTVec.emplace_back(EdgePointWithType(isoEdgePoints[0][i]));
		edgePWTVec.emplace_back(EdgePointWithType(triVerts[1]->pos, meshVertsType[triVerts[1]->index]));
		for (int i = 0; i < isoEdgePoints[1].size(); ++i) edgePWTVec.emplace_back(EdgePointWithType(isoEdgePoints[1][i]));
		edgePWTVec.emplace_back(EdgePointWithType(triVerts[2]->pos, meshVertsType[triVerts[2]->index]));
		for (int i = 0; i < isoEdgePoints[2].size(); ++i) edgePWTVec.emplace_back(EdgePointWithType(isoEdgePoints[2][i]));

		int numEdgeCCWPoints = edgePWTVec.size();
		for (int i = 0; i < numEdgeCCWPoints; ++i)
		{
			bdSegQueue.push(std::make_pair(i, (i + 1) % numEdgeCCWPoints));
			//edgeCCWPointsIdxQueue.push(i);
		}

		std::vector<CellTriangle> facetAllCDTs;
		std::vector<typename CDT_Cell::_Point> _bdPoints;
		std::vector<typename CDT_Cell::_Edge> _bdEdges;
		while (!bdSegQueue.empty())
		{
			VertIdxPair fro_seg = bdSegQueue.front();
			bdSegQueue.pop();
			if (isVisSeg.count(fro_seg)) continue;
			isVisSeg.insert(fro_seg);

			// 1.2 Set boundary of current cell
			_bdPoints.clear(); _bdEdges.clear();

			VERT_TYPE cellType = VERT_TYPE::EQUAL;
			int beg_cellPointIdx = fro_seg.second;
			int cur_cellPointIdx = beg_cellPointIdx; // 始终存储的是边界点的索引

			int last_cellIsoPointIdx = -1;
			do {
				// 如果边界点是某段等值线上的点
				ApolloniusDiagramPoint_3 curPoint = edgePWTVec[cur_cellPointIdx].pos;
				if (isoEdgePointToIdx.find(curPoint) != isoEdgePointToIdx.end() &&
					isoEdgePointToIdx.at(curPoint).first != last_cellIsoPointIdx)
				{
					if (cellType == VERT_TYPE::EQUAL && edgePWTVec[cur_cellPointIdx].type != VERT_TYPE::EQUAL)
					{
						cellType = (VERT_TYPE)!int(edgePWTVec[cur_cellPointIdx].type);
					}

					// 让它走完一整段等值线
					// 获得在edgeTable中的索引
					int beg_isoVertTableIdx = isoEdgePointToIdx.at(curPoint).second;
					int cur_isoVertTableIdx = beg_isoVertTableIdx;
					do
					{
						const ApolloniusDiagramPoint_3 iso_vert_1 = isoLineVerts[cur_isoVertTableIdx];
						cur_isoVertTableIdx = edgeTable.at(cur_isoVertTableIdx)[0];
						const ApolloniusDiagramPoint_3 iso_vert_2 = isoLineVerts[cur_isoVertTableIdx];
						_bdEdges.emplace_back(Segment_3(
							Point_3(iso_vert_1.x(), iso_vert_1.y(), iso_vert_1.z()),
							Point_3(iso_vert_2.x(), iso_vert_2.y(), iso_vert_2.z())
						));

						curPoint = iso_vert_2;
					} while (isoEdgePointToIdx.find(curPoint) == isoEdgePointToIdx.end()); // 直到等值线的某个点是边上的边界点为止
					const int end_isoVertIdx = isoEdgePointToIdx.at(curPoint).first; // 等值线的末端点索引

					if (cellType != VERT_TYPE::EQUAL)
					{
						edgePWTVec[cur_cellPointIdx].type = cellType; // cur_cellPointIdx是等值线的起始端点索引
						edgePWTVec[end_isoVertIdx].type = cellType;
					}

					cur_cellPointIdx = end_isoVertIdx;
					last_cellIsoPointIdx = cur_cellPointIdx;
				}
				else
				{
					// 如果边界点是三角形的端点，则该cell的type由该端点确定
					if (cellType == VERT_TYPE::EQUAL) cellType = edgePWTVec[cur_cellPointIdx].type;

					// 走向下一条边界边，并将边界边的vis设为true
					int next_cellPointIdx = (cur_cellPointIdx + 1) % numEdgeCCWPoints;
					const auto next_seg = std::make_pair(cur_cellPointIdx, next_cellPointIdx);
					if (!isVisSeg.count(next_seg)) isVisSeg.insert(next_seg);

					const ApolloniusDiagramPoint_3 bd_edge_vert_1 = edgePWTVec[cur_cellPointIdx].pos;
					const ApolloniusDiagramPoint_3 bd_edge_vert_2 = edgePWTVec[next_cellPointIdx].pos;
					_bdEdges.emplace_back(Segment_3(
						Point_3(bd_edge_vert_1.x(), bd_edge_vert_1.y(), bd_edge_vert_1.z()),
						Point_3(bd_edge_vert_2.x(), bd_edge_vert_2.y(), bd_edge_vert_2.z())
					));

					cur_cellPointIdx = next_cellPointIdx;
				}
			} while (beg_cellPointIdx != cur_cellPointIdx);

			//// 2. 对每个区域计算CDT
			//CDT_Cell cell(_bdEdges);
			//const auto& cdtForCurCell = cdtAdaptor.computeCDForBoundary(cell);

			//// 3. push结果
			//facetAllCDTs.emplace_back(cdtForCurCell);
		}

		return facetAllCDTs;
	}

	/* Methods for Our Algorithm */
	/**
	* @brief: 对mesh每条边进行固定数量的采样
	* @return: 采样点数量
	*/
	int MSCuttingModel::samplePointPerEdge()
	{
		if (numSamplesPerEdge <= 0)
		{
			LOG::qpWarn("The number of sample points is smaller than 0, the process of sampling is terminated.");
			return 0;
		}

		using Point3 = typename HEVert::Point3; // 其实就是Eigen::Vector3
		using EdgeDir = typename Vector3;

		// 获得所有model edge
		const auto& modelEdges = this->getMEdgeVec();
		int numSplits = numSamplesPerEdge + 1;

		// 对每条边计算numSamplesPerEdge个采样点
		for (const auto& mEdge : modelEdges)
		{
			// 两个端点坐标
			const Point3 vert_1 = mEdge->firstVertex()->pos;
			const Point3 vert_2 = mEdge->secondVertex()->pos;
			// 关联的两个(或一个)面
			int facet_1 = -1, facet_2 = -1;
			if (!mEdge->halfEdge()->isBoundary())
				facet_1 = mEdge->halfEdge()->boundFace()->index;
			if (!mEdge->halfEdge()->oppoHEdge()->isBoundary())
				facet_2 = mEdge->halfEdge()->oppoHEdge()->boundFace()->index;

			const EdgeDir edgeDir = vert_2 - vert_1;

			Point3 preSamplePointPos = vert_1;
			Scalar preSamplePointVal = scalarFunc.val(vert_1);
			SamplePoint preSamplePoint(preSamplePointPos, mEdge->index, preSamplePointVal);

			int preInsertedIdx = -1;
			//bool preIsLargeZero = (preSamplePointVal > 0);
			for (int k = -1; k < numSamplesPerEdge + 1; ++k) // +1 是为了涵盖另一个端点vert_2
			{
				// 计算采样点位置
				Point3 curSamplePointPos = vert_1 + (k + 1) * edgeDir / numSplits; // 先乘后除，或许可以减小点误差
				// 计算当前采样点的值
				Scalar curSamplePointVal = scalarFunc.val(curSamplePointPos);

				SamplePoint curSamplePoint(curSamplePointPos, mEdge->index, curSamplePointVal);

				//bool curIsLargeZero = (curSamplePointVal > 0); // TODO: = 0 的情况该怎么处理？判断是不是切点
				//if (curIsLargeZero ^ preIsLargeZero)
				{
					/*if (preInsertedIdx != k - 1 || k == 0) validSamplePoints.emplace_back(preSamplePoint);
					validSamplePoints.emplace_back(curSamplePoint);*/

					/*std::cout << "preSamplePoint = " << preSamplePointPos.transpose() << " val = " << preSamplePointVal << "\n";
					std::cout << "curSamplePoint = " << curSamplePointPos.transpose() << " val = " << curSamplePointVal << "\n";
					LOG::qpSplit();*/

					// push结果
					if (facet_1 != -1)
					{
						///*if (preInsertedIdx != k - 1 || k == 0) */sampleFacets[facet_1].aroundSamplePoints.emplace_back(preSamplePoint);
						sampleFacets[facet_1].aroundSamplePoints.emplace_back(curSamplePoint);

					}
					if (facet_2 != -1)
					{
						///*if (preInsertedIdx != k - 1 || k == 0) */sampleFacets[facet_2].aroundSamplePoints.emplace_back(preSamplePoint);
						sampleFacets[facet_2].aroundSamplePoints.emplace_back(curSamplePoint);
					}

					preInsertedIdx = k;
				}

				/*preSamplePointPos = curSamplePointPos;
				preIsLargeZero = curIsLargeZero;
				preSamplePoint = curSamplePoint;*/

				samplePoints.emplace_back(curSamplePoint);
			}
		}

		//assert(samplePoints.size() == (numSamplesPerEdge * modelEdges.size()));

		return samplePoints.size();
	}

	/**
	* @brief: 计算Apollonius Graph(对所有采样面)
	* @param out: 输出流
	*/
	void MSCuttingModel::computeApolloniusGraph(std::ofstream& out)
	{
		int globalOutVertIdx = 1; // .obj的顶点下标从1开始
		for (int i = 0; i < numMeshFaces; ++i)
		{
			// 保存当前面的ApolloniusDiagram顶点的信息
			std::vector<ApolloniusDiagramPoint_3> curFacetADPoints;

			// 保存当前面的ApolloniusDiagram边的信息(通过顶点的输出索引进行存储)
			std::vector<std::pair<int, int>> curFacetADLines;
			std::unordered_map<int, std::vector<int>> edgeTable; // 保存该面power diagram点的邻边表
			computeApolloniusGraphForFacet(i, globalOutVertIdx, curFacetADPoints, curFacetADLines, edgeTable);

			std::array<std::vector<ApolloniusDiagramPoint_3>, 3> edgePoints; // 保存在边上的(后处理过的)点
			std::vector<ApolloniusDiagramPoint_3> projCurFacetADPoints = postProcessFacetPoints(i, curFacetADPoints, edgePoints);

			// 输出Apollonius Diagram中的所有全局坐标点
			for (const auto& adPoint : projCurFacetADPoints)
			{
				out << "v " << adPoint.transpose() << std::endl;
			}

			// 输出Apollonius Diagram中的所有边
			for (const auto& vertIdxPair : curFacetADLines)
			{
				out << "l " << vertIdxPair.first << " " << vertIdxPair.second << std::endl;
			}
		}
	}

	/* Visualization */
	/**
	* @brief: Visualization of sample points on each edge
	* @param out: 可视化结果的输出流
	*/
	void MSCuttingModel::outputSamplePoints(std::ofstream& out_1, std::ofstream& out_2)
	{
		if (out_1)
		{
			for (const auto& samplePoint : samplePoints)
			{
				out_1 << "v " << samplePoint.pos.transpose() << std::endl;
			}
		}

		if (out_2)
		{
			for (const auto& validSamplePoint : validSamplePoints)
			{
				out_2 << "v " << validSamplePoint.pos.transpose() << std::endl;
			}
		}
	}

	/* API for user */
	/**
	* @brief: 暴露给用户调用的接口，用于运行整体算法
	* @param filename: 输出结果所保存的文件位置
	* @return: 算法运行成功/错误
	*/
	bool MSCuttingModel::launch(const std::string& ad_vis_file)
	{
		str_util::checkDir(ad_vis_file);
		std::ofstream ad_vis_out(ad_vis_file);
		if (!ad_vis_out) { LOG::qpError("I/O: File ", ad_vis_file.c_str(), " could not be opened!"); return false; }

		samplePointPerEdge();

		const std::string sample_vis_file = str_util::concatFilePath(VIS_DIR, modelName, (std::string)"sample_points.obj");
		str_util::checkDir(sample_vis_file);
		std::ofstream sample_vis_out(sample_vis_file);
		if (!sample_vis_out) { LOG::qpError("I/O: File ", sample_vis_file.c_str(), " could not be opened!"); return false; }

		const std::string valid_sample_vis_file = str_util::concatFilePath(VIS_DIR, modelName, (std::string)"valid_sample_points.obj");
		str_util::checkDir(valid_sample_vis_file);
		std::ofstream valid_sample_vis_out(valid_sample_vis_file);
		if (!valid_sample_vis_out) { LOG::qpError("I/O: File ", sample_vis_file.c_str(), " could not be opened!"); return false; }

		LOG::qpInfo("Output Sample Points to ", std::quoted(ad_vis_file), " ...");
		outputSamplePoints(sample_vis_out, valid_sample_vis_out);
		sample_vis_out.close();

		LOG::qpInfo("Output Apollonius Diagram to ", std::quoted(ad_vis_file), " ...");
		computeApolloniusGraph(ad_vis_out);
		ad_vis_out.close();

		return true;
	}

	/* Test APIs for us */
	bool MSCuttingModel::testSamplingPoints(std::ofstream& out_1, std::ofstream& out_2)
	{
		samplePointPerEdge();

		outputSamplePoints(out_1, out_2);

		return true;
	}

	bool MSCuttingModel::testLocalGlobalTransform(int facetIdx)
	{
		samplePointPerEdge();

		const auto& curSampleFacet = sampleFacets[facetIdx];
		const auto& aroundEndVerts = this->getFaceVec()[facetIdx]->aroundVerts();
		const auto& aroundSamplePoints = curSampleFacet.aroundSamplePoints;
		// 输出所有原始采样点的坐标
		LOG::qpTest("Original points:");
		for (const auto& point_3 : aroundEndVerts)
		{
			LOG::qpNormal("v: ", point_3->pos.transpose());
		}
		for (const auto& point_3 : aroundSamplePoints)
		{
			LOG::qpNormal("v: ", point_3.pos.transpose());
		}

		LOG::qpSplit();

		Polygon_2 facetBoundary;
		updateFacetLocalCoord(facetIdx, facetBoundary);

		// 输出转换到局部坐标系下的点坐标
		LOG::qpTest("Local coordinate of original points:");
		for (const auto& point_2 : curSampleFacet.sites)
		{
			LOG::qpNormal("local v: ", point_2.pos);
		}

		LOG::qpSplit();

		// 再转换回去
		LOG::qpTest("Transforming these local points back to global:");
		for (int i = 0; i < curSampleFacet.sites.size(); ++i)
		{
			const auto point_2 = curSampleFacet.sites[i];
			const auto point_3 = getGlobalCoordInFacet(facetIdx, point_2.pos);
			LOG::qpNormal("global v: ", point_3.transpose());

			//OFFSET_ENSURE(point_3.isApprox(curSampleFacet.aroundSamplePoints[i].pos, 1e-10));
		}

		return true;
	}

	bool MSCuttingModel::testComputeADForFacet(int facetIdx, std::ofstream& out)
	{
		samplePointPerEdge();

		// 保存当前面的ApolloniusDiagram顶点的信息
		std::vector<ApolloniusDiagramPoint_3> facetADPoints;

		// 保存当前面的ApolloniusDiagram边的信息(通过顶点的输出索引进行存储)
		std::vector<std::pair<int, int>> facetADLines;

		int globalOutVertIdx = 1; // .obj的顶点下标从1开始

		std::unordered_map<int, std::vector<int>> edgeTable; // 保存该面power diagram点的邻边表
		int numADLines = computeApolloniusGraphForFacet(facetIdx, globalOutVertIdx, facetADPoints, facetADLines, edgeTable);

		LOG::qpTest("The number of edges of Apollonius Diagram = ", numADLines);

		// 输出facet
		const std::string tri_out_file = R"(F:\VisualStudioProgram\MeshScalarCutting\vis\cube\tri.obj)";
		std::ofstream tri_out(tri_out_file);
		for (const auto& triPoint : this->getFaceVec()[facetIdx]->aroundVerts())
		{
			tri_out << "v " << triPoint->pos.transpose() << "\n";
		}
		tri_out << "f 1 2 3\n";

		// 输出站点
		const std::string site_out_file = R"(F:\VisualStudioProgram\MeshScalarCutting\vis\cube\sites.obj)";
		std::ofstream site_out(site_out_file);
		for (const auto& site : sampleFacets[facetIdx].sites)
		{
			site_out << "v " << getGlobalCoordInFacet(facetIdx, site.pos).transpose() << "\n";
		}
		tri_out.close();
		site_out.close();

		// 输出Apollonius Diagram中的所有全局坐标点
		for (const auto& adPoint : facetADPoints)
		{
			out << "v " << adPoint.transpose() << std::endl;
		}

		// 输出Apollonius Diagram中的所有边
		for (const auto& vertIdxPair : facetADLines)
		{
			out << "l " << vertIdxPair.first << " " << vertIdxPair.second << std::endl;
		}

		return true;
	}

} // namespace core
NAMESPACE_END(mscut)