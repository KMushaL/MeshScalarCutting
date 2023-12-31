﻿#pragma once
#include "MSCuttingModel.hpp"
#include "detail/Handler.hpp"
#include <map>
#include <fstream>
#include <igl/AABB.h>
#include <utils/Donut.hpp>

NAMESPACE_BEGIN(mscut)
namespace core
{
	///////////////////////
	// Algorithm Details //
	///////////////////////
	/**
	* @brief: 对mesh每条边进行固定数量的采样
	* @return: 采样点数量
	*/
	int MSCuttingModel::samplePointPerEdge()
	{
		if (numSamplesPerEdge < 1)
		{
			numSamplesPerEdge = 1;
			LOG::qpWarn("The number of sample points for each facet is smaller than zero.");
		}

		using Point3 = typename HEVert::Point3; // 其实就是Eigen::Vector3
		using EdgeDir = typename Vector3;

		// 获得所有model edge
		const auto& modelEdges = this->getMEdgeVec();
		int numSplits = numSamplesPerEdge - 1/* + 1*/;

		std::set<Point3> samplePointSet;

		// 对每条边计算numSamplesPerEdge个采样点
		for (const auto& mEdge : modelEdges)
		{
			// 两个端点坐标
			const Point3 vert_1 = mEdge->firstVertex()->pos;
			const Point3 vert_2 = mEdge->secondVertex()->pos;

			// 关联的两个(或一个)面
			int facet_1 = -1, facet_2 = -1;
			if (!mEdge->halfEdge()->isBoundary()) facet_1 = mEdge->halfEdge()->boundFace()->index;
			if (!mEdge->halfEdge()->oppoHEdge()->isBoundary()) facet_2 = mEdge->halfEdge()->oppoHEdge()->boundFace()->index;

			const EdgeDir edgeDir = vert_2 - vert_1;
			for (int k = -1; k < numSamplesPerEdge - 1; ++k) // +1 是为了涵盖另一个端点vert_2
			{
				// 计算采样点位置
				Point3 curSamplePointPos = vert_1 + (k + 1) * edgeDir / numSplits; // 先乘后除，或许可以减小点误差
				//if (std::fabs(curSamplePointPos.x() - 1) < 1e-6 || std::fabs(curSamplePointPos.y()) < 1e-6) continue;
				int s = 0;
				for (s; s < singulars.size(); ++s)
					if ((curSamplePointPos - singulars[s]).norm() <= SINGULAR_DIS_EPSILON) break;
				if (s != singulars.size()) continue;

				// 计算当前采样点的值
				Vector3 curSamplePointGrad = scalarFunc.grad(curSamplePointPos);
				if (std::isinf(curSamplePointGrad(0)) || std::isinf(curSamplePointGrad(1)) || std::isinf(curSamplePointGrad(2)) ||
					std::isnan(curSamplePointGrad(0)) || std::isnan(curSamplePointGrad(1)) || std::isnan(curSamplePointGrad(2)) ||
					curSamplePointGrad.norm() < 1e-9) continue;
				Scalar curSamplePointVal = scalarFunc.val(curSamplePointPos);
				SamplePoint curSamplePoint(curSamplePointPos, mEdge->index, curSamplePointVal);

				double alpha = std::fabs(curSamplePointVal / (curSamplePointGrad.squaredNorm()));

				if (facet_1 != -1 && !sampleFacets[facet_1].aroundSamplePointsSet.count(curSamplePointPos))
				{
					sampleFacets[facet_1].aroundSamplePointsSet.insert(curSamplePointPos);
					sampleFacets[facet_1].aroundSamplePoints.emplace_back(curSamplePoint);
				}
				if (facet_2 != -1 && !sampleFacets[facet_2].aroundSamplePointsSet.count(curSamplePointPos))
				{
					sampleFacets[facet_2].aroundSamplePointsSet.insert(curSamplePointPos);
					sampleFacets[facet_2].aroundSamplePoints.emplace_back(curSamplePoint);
				}

				if (!samplePointSet.count(curSamplePointPos))
				{
					samplePointSet.insert(curSamplePointPos);
					samplePoints.emplace_back(curSamplePoint);
				}
			}
		}

		return samplePoints.size();
	}

	/**
	* @brief: 对单个面计算局部的二维坐标以及权重
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

			meshVertsType[aroundEndVerts[i]->index] = (CELL_TYPE)(scalarFunc.val(endVert_dim3) > 0); // TODO: 等于0时的数值误差处理

			/*site.weight = -scalarFunc.val(endVert_dim3) / scalarFunc.grad(endVert_dim3).norm();
			curSampleFacet.sites.emplace_back(site);*/

			// 更新边界信息
			facetBoundary.push_back(site.pos);
		}

		auto lineSearch = [&](const Vector3& x_p, const Vector3& grad_p, const Vector3& drt, const double& f_val,
			Vector3& x, Vector3& grad_x, double& f_x)
			{
				double step = 1e-1 * 1.0 / (grad_p.norm()); // 初始步长

				const double ftol = 1e-4;
				const double wolfe = 0.9;

				// Projection of gradient on the search direction
				const Scalar f_x_init = f_val;
				const Vector3 grad_x_init = grad_p;
				Scalar dg_init = grad_p.dot(drt); // < 0
				Scalar test_decr = ftol * dg_init; // 保证test_decr < 0

				// Upper and lower end of the current line search range
				Scalar step_lo = 0,
					step_hi = std::numeric_limits<Scalar>::infinity();
				Scalar min_step = Scalar(1e-20);
				Scalar max_step = Scalar(1e+20);

				int max_linesearch = 5;
				const Vector3 xp = x;
				for (int iter = 0; iter < max_linesearch; ++iter)
				{
					//std::cout << "iter = " << iter << std::endl;
					// x_{k+1} = x_k + step * d_k
					x = xp + step * drt;
					// Evaluate this candidate
					f_x = 0.5 * scalarFunc.val(x) * scalarFunc.val(x);
					grad_x = f_x * scalarFunc.grad(x);

					if (std::isinf(grad_x(0)) || std::isinf(grad_x(1)) || std::isinf(grad_x(2)) ||
						std::isnan(grad_x(0)) || std::isnan(grad_x(1)) || std::isnan(grad_x(2)))
					{
						if (iter == max_linesearch - 1)
						{
							x = xp;
							f_x = f_x_init;
							grad_x = grad_x_init;
						}
						continue;
					}
					/*double t_f_x = f_x;
					Vector3 t_grad_x = grad_x;
					if (f_val < 0)
					{
						t_f_x = -t_f_x;
						t_grad_x = -t_grad_x;
					}*/

					if (f_x > f_x_init + step * test_decr || (f_x != f_x))
					{
						step_hi = step;
					}
					else
					{
						// Armijo condition is met
						//break;

						Scalar dg = grad_x.dot(drt);

						if (dg < wolfe * dg_init)
						{
							step_lo = step;
						}
						else
						{
							// Regular Wolfe condition is met
							//break;

							if (dg > -wolfe * dg_init)
							{
								step_hi = step;
							}
							else
							{
								// Strong Wolfe condition is met
								break;
							}
						}
					}

					assert(step_lo < step_hi);

					if (step < min_step) break;
					//throw std::runtime_error("the line search step became smaller than the minimum value allowed");

					if (step > max_step) break;
					//throw std::runtime_error("the line search step became larger than the maximum value allowed");

				// continue search in mid of current search range
					step = std::isinf(step_hi) ? 2 * step : step_lo / 2 + step_hi / 2;
				}

				return step;
			};

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

			if (f_grad.isApprox(Vector3(0, 0, 0), 1e-9)) site.weight = 0;
			else
			{
				/*if (f_val == 0) {
					site.weight = .0;
					info_out << "v " << sampleVert_dim3(0) << " " << sampleVert_dim3(1) << " " << site.weight << "\n";
					continue;
				}*/

				Vector3 x = sampleVert_dim3;
				Vector3 last_x = x;

				double g_x = 0.5 * f_val * f_val;
				double last_g_x = g_x;

				Vector3 grad_x = f_val * f_grad;
				Vector3 last_grad_x = grad_x;

				Vector3 drt = -grad_x;
				int iter = 1;
				while (iter <= WEIGHT_MAX_ITER)
				{
					double step = lineSearch(last_x, last_grad_x, drt, last_g_x, x, grad_x, g_x);

					grad_x = scalarFunc.val(x) * scalarFunc.grad(x);
					if ((last_grad_x - grad_x).norm() < 1e-6) break;

					last_x = x;
					last_grad_x = grad_x;
					last_g_x = g_x;

					drt = -last_grad_x;

					++iter;
				}
				site.weight = (x - sampleVert_dim3).squaredNorm();
			}

			//std::cout << "pos: " << sampleVert_dim3.transpose() << ", weight = " << site.weight << std::endl;
			//info_out << "v " << sampleVert_dim3(0) << " " << sampleVert_dim3(1) << " " << site.weight << "\n";

			curSampleFacet.sites.emplace_back(site);
		}

		return curSampleFacet.sites.size();
	}

	/**
	* @param facetIdx: 待转换的局部坐标点所在的面的索引
	* @param point_2: 待转换的局部坐标点
	* @return: 转换后的全局坐标
	*/
	Eigen::Vector3d MSCuttingModel::getGlobalCoordInFacet(int facetIdx, const PowerDiagramPoint_2& point_2)
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
	* @brief: 为facetIdx这个面计算Power Diagram
	* @param facetIdx: 面的索引
	* @return: Power Diagram边的数量
	*/
	int MSCuttingModel::computePDForFacet(int facetIdx,
		int& globalOutVertIdx,
		std::vector<PowerDiagramPoint_3>& pdPoints,
		std::vector<std::pair<int, int>>& pdLines,
		std::unordered_map<int, std::vector<int>>& edgeTable,
		std::map<PowerDiagramPoint_3, int>& pdPointToOutIdx)
	{
		const int originGlobalOutVertIdx = globalOutVertIdx;

		// 从local point到输出下标的映射，主要用于防止输出重复的点
		// 为什么不弄个全局的(也就是对所有面的)？为了后续稍微好并行(对每个面同时计算Apollonius Diagram)一点
		std::map<Point_2, int> pointToGlobalOutIdx;

		// 计算局部坐标和边界
		Polygon_2 facetBoundary;
		updateFacetLocalCoord(facetIdx, facetBoundary);

		// 计算Power Diagram
		auto adLineSegmentList = pdAdaptor.computePDForBoundary(sampleFacets[facetIdx].sites, facetBoundary);

		if (adLineSegmentList.empty())
		{
			const auto& aroundVerts = this->getFaceVec()[facetIdx]->aroundVerts();
			std::vector<Vector3> verts; std::vector<int> vertIdx;
			CELL_TYPE type = meshVertsType[aroundVerts[0]->index];
			int meshVertIdx = (type == CELL_TYPE::INSIDE) ? insideMeshVertIdx : outsideMeshVertIdx;
			for (int i = 0; i < 3; ++i)
			{
				verts.emplace_back(aroundVerts[i]->pos);
				vertIdx.emplace_back(meshVertIdx + i);
			}

			if (type == CELL_TYPE::INSIDE) insideMeshVertIdx += verts.size(), insideMesh.emplace_back(verts, vertIdx);
			else if (type == CELL_TYPE::OUTSIDE) outsideMeshVertIdx += verts.size(), outsideMesh.emplace_back(verts, vertIdx);

			return 0;
		}

		// 遍历每条边，转回全局坐标
		for (const auto& lineSegemnt : adLineSegmentList)
		{
			const PowerDiagramPoint_2 vert_1 = lineSegemnt.vertex(0);
			const PowerDiagramPoint_2 vert_2 = lineSegemnt.vertex(1);

			int outVertIdx_1, outVertIdx_2; // 当前这两个点应输出的下标
			if (pointToGlobalOutIdx.count(vert_1)) outVertIdx_1 = pointToGlobalOutIdx.at(vert_1);
			else
			{
				PowerDiagramPoint_3 pdPoint = getGlobalCoordInFacet(facetIdx, vert_1);
				pdPoints.emplace_back(pdPoint);

				outVertIdx_1 = globalOutVertIdx++;
				pointToGlobalOutIdx.insert(std::make_pair(vert_1, outVertIdx_1));
				pdPointToOutIdx.insert(std::make_pair(pdPoint, outVertIdx_1));
			}

			if (pointToGlobalOutIdx.count(vert_2)) outVertIdx_2 = pointToGlobalOutIdx.at(vert_2);
			else
			{
				PowerDiagramPoint_3 pdPoint = getGlobalCoordInFacet(facetIdx, vert_2);
				pdPoints.emplace_back(pdPoint);

				outVertIdx_2 = globalOutVertIdx++;
				pointToGlobalOutIdx.insert(std::make_pair(vert_2, outVertIdx_2));
				pdPointToOutIdx.insert(std::make_pair(pdPoint, outVertIdx_2));
			}

			// 保存这条边的两个顶点的索引
			// TODO: 注意，我们这里假定每个点当且仅当只会和一个点相邻
			pdLines.emplace_back(std::make_pair(outVertIdx_1, outVertIdx_2));

			// edgeTable保存的是局部索引
			edgeTable[outVertIdx_1 - originGlobalOutVertIdx].emplace_back(outVertIdx_2 - originGlobalOutVertIdx);
			edgeTable[outVertIdx_2 - originGlobalOutVertIdx].emplace_back(outVertIdx_1 - originGlobalOutVertIdx);
		}

		return pdLines.size();
	}

	/* Post-processing */
	std::vector<MSCuttingModel::PowerDiagramPoint_3>
		MSCuttingModel::postProcessFacetPoints(const int faceIdx,
			const std::vector<std::vector<PowerDiagramPoint_3>>& facetPDPoints,
			const std::map<PowerDiagramPoint_3, int>& pdPointToOutIdx,
			std::array<std::vector<PowerDiagramPoint_3>, 3>& edgePoints,
			std::map<PowerDiagramPoint_3, int>& edgePointToIdx,
			std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<PowerDiagramPoint_3, int>>>>& edgeToPointsInNeighFace)
	{
		int numPDPoints = facetPDPoints[faceIdx].size();
		std::vector<PowerDiagramPoint_3> resPoints(numPDPoints);

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
		const auto hEdge = this->getFaceVec()[faceIdx]->halfEdge();
		std::array<int, 3> faceEdgeIdx;
		auto _hEdge = hEdge;
		std::array<bool, 3> isConsistencyOnEdge = { false, false, false };
		for (int i = 0; i < 3; ++i) {
			// 邻居面没有PD经过，则不需要使用投影梯度法
			if (_hEdge->oppoHEdge()->isBoundary() ||
				facetPDPoints[_hEdge->oppoHEdge()->boundFace()->index].empty())
				isConsistencyOnEdge[i] = true;

			faceEdgeIdx[i] = _hEdge->mEdge()->index;
			_hEdge = _hEdge->nextHEdge();
		}

		const std::array<VertIdx, 3> triEdges = {
			VertIdx(0, 1),
			VertIdx(1, 2),
			VertIdx(2, 0)
		};
		const auto mEdgeVec = this->getMEdgeVec();
		const std::array<Vector3, 3> triEdgesDir = {
			(triVerts[1]->pos - triVerts[0]->pos).normalized(),
			(triVerts[2]->pos - triVerts[1]->pos).normalized(),
			(triVerts[0]->pos - triVerts[2]->pos).normalized()
		};

#pragma omp parallel for
		for (int i = 0; i < numPDPoints; ++i) {
			const PowerDiagramPoint_3 originalPDPoint = facetPDPoints[faceIdx][i];

			// 判断点是否在边上
			int onEdgeIdx = -1;
			int globalOutIdx = -1;
			for (int i = 0; i < 3; ++i) {
				const Vector3 v1 = triVerts[triEdges[i].first]->pos;

				// TODO: Exact Predicates
				if (!isConsistencyOnEdge[i] &&
					onEdgeIdx == -1 && triEdgesDir[i].cross((originalPDPoint - v1).normalized()).norm() < 1e-9) {
					onEdgeIdx = i;
					globalOutIdx = pdPointToOutIdx.at(originalPDPoint);
					break;
				}
			}

			VectorX tri_sqrD; VectorXi tri_I; MatrixX tri_C; MatrixX projQ(1, 3);

			// for computing projection point
			PowerDiagramPoint_3 lastPoint = originalPDPoint;
			PowerDiagramPoint_3 projPoint = lastPoint;
			Scalar lastVal; Vector3 lastGrad; double alpha;

			int iter = 0;
			while (iter < POST_PROCESSING_MAX_ITER) {
				lastVal = scalarFunc.val(lastPoint);
				if (std::fabs(lastVal) < NUMERICAL_EPSILON)
				{
					projPoint = lastPoint;
					break;
				}

				lastGrad = scalarFunc.grad(lastPoint);

				if (std::isinf(lastGrad.norm()) || std::isnan(lastGrad.norm()) ||
					std::fabs(lastGrad.norm()) < NUMERICAL_EPSILON)
				{
					projPoint = lastPoint;
					break;
				}

				++iter;
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
					double max_move;
					Vector3 proj_grad;
					double t = (alpha * lastGrad).dot(triEdgesDir[onEdgeIdx]);
					if (t >= 0)
					{
						max_move = std::min(t, (lastPoint - (triVerts[triEdges[onEdgeIdx].second])->pos).norm());
						proj_grad = triEdgesDir[onEdgeIdx] * max_move;
					}
					else
					{
						max_move = std::min(-t, (lastPoint - (triVerts[triEdges[onEdgeIdx].first])->pos).norm());
						proj_grad = -triEdgesDir[onEdgeIdx] * max_move;
					}
					projPoint = proj_grad + lastPoint;
				}

				if ((projPoint - lastPoint).norm() < PROJECT_EPSILON) break;
				else lastPoint = projPoint;
			}

			resPoints[i] = projPoint;

			// TODO: 将去重工作移到外面，提升并行效率
//#pragma omp critical
			{
				// 目前对后处理后的所有点做了一个去重工作
				//if (edgePointToIdx.find(projPoint) == edgePointToIdx.end())
				{
					//resPoints[i] = projPoint;
					if (onEdgeIdx != -1)
					{
						// 用于连接相邻面的等值线
						edgeToPointsInNeighFace[faceEdgeIdx[onEdgeIdx]][faceIdx].emplace_back(projPoint, globalOutIdx);

						//// 用于构建CDT，目前用不到先注释掉
						//edgePoints[onEdgeIdx].emplace_back(projPoint);
						//edgePointToIdx[projPoint] = i; // 得到该边界点在哪条边上的信息以及在edgeTable中的索引
					}
				}

				/*if (i <= 50) LOG::qpInfo("#i = ", i, ", went ", iter, " iterations.");
				else LOG::qpWarn("#i = ", i, ", went ", iter, " iterations.");*/
			}
		}

		//// 将边上的点统一按逆时针排序
		//struct CCWSort {
		//	PowerDiagramPoint_3 edgeBaseVert;
		//	Vector3 edgeDir;
		//	CCWSort(const PowerDiagramPoint_3& _edgeBaseVert, const Vector3& _edgeDir) :edgeBaseVert(_edgeBaseVert), edgeDir(_edgeDir) {}
		//
		//	bool operator()(const PowerDiagramPoint_3& lhs, const PowerDiagramPoint_3& rhs)
		//	{
		//		return ((lhs - edgeBaseVert).dot(edgeDir) < (rhs - edgeBaseVert).dot(edgeDir));
		//	}
		//};
		//donut::Loop<int, 3>([&](int i) {
		//	if (!edgePoints[i].empty()) {
		//		CCWSort ccwSort(triVerts[i]->pos, triEdgesDir[i]);
		//		std::sort(edgePoints[i].begin(), edgePoints[i].end(), ccwSort);
		//	}
		//	});

		return resPoints;
	}

	std::vector<typename MSCuttingModel::CellTriangle>
		MSCuttingModel::computeCDTForFacet(const int faceIdx,
			const std::vector<PowerDiagramPoint_3>& isoLineVerts,
			const std::unordered_map<int, std::vector<int>>& edgeTable,
			const std::array<std::vector<PowerDiagramPoint_3>, 3>& isoEdgePoints,
			const std::map<PowerDiagramPoint_3, int>& isoEdgePointToIdx)
	{

		// 1. 搜寻边界边并确定区域
		// 1.1 按逆时针将所有边上的点(包括面的端点)按逆时针排序
		/*std::queue<int> edgeCCWPointsIdxQueue;
		std::vector<PowerDiagramPoint_3> edgeCCWPoints;*/
		using VertIdxPair = std::pair<int, int>;
		std::queue<VertIdxPair> bdSegQueue;
		std::set<VertIdxPair> isVisSeg;

		// Initialize cycle(counter clock-wise) points
		struct EdgePointWithType {
			PowerDiagramPoint_3 pos;
			CELL_TYPE type = CELL_TYPE::EQUAL;

			EdgePointWithType(const PowerDiagramPoint_3& _pos) :pos(_pos) {}
			EdgePointWithType(const PowerDiagramPoint_3& _pos, const CELL_TYPE& _type) :pos(_pos), type(_type) {}
		};
		std::vector<EdgePointWithType> edgePWTVec;

		// TODO: 目前没有对edgeCyclePoints做一个去重工作，也就是说即使三角形端点和等值线上的边界点几乎是一个点但目前也认为是两个点
		int prefixSum = 0; std::map<PowerDiagramPoint_3, int> isoEdgePointToPWTIdx; // 存储等值线的边界点到edgePWTVec中的索引
		const auto& triVerts = this->getFaceVec()[faceIdx]->aroundVerts();
		auto initCCWEdgePointData = [&](int triVertIdx) {
			edgePWTVec.emplace_back(triVerts[triVertIdx]->pos, meshVertsType[triVerts[triVertIdx]->index]);
			int numIsoEdgePoints = isoEdgePoints[triVertIdx].size(); // 在第triVertIdx条边上的等值线点数量
			for (int i = 0; i < numIsoEdgePoints; ++i)
			{
				const auto& isoEdgePoint = isoEdgePoints[triVertIdx][i];
				edgePWTVec.emplace_back(isoEdgePoint);
				isoEdgePointToPWTIdx[isoEdgePoint] = i + prefixSum + 1;
			}
			prefixSum += (numIsoEdgePoints + 1);
			};
		donut::Loop<int, 3>([&](int i) {
			initCCWEdgePointData(i);
			});

		int numEdgeCCWPoints = edgePWTVec.size();
		for (int i = 0; i < numEdgeCCWPoints; ++i)
		{
			bdSegQueue.push(std::make_pair(i, (i + 1) % numEdgeCCWPoints));
			//edgeCCWPointsIdxQueue.push(i);
		}

		// 保存结果的vectors
		std::vector<CellTriangle> facetAllCDTs;
		std::vector<typename CDT_Cell::_Point> constrainedPoints;
		std::vector<typename CDT_Cell::_Edge> constrainedEdges;

		// 添加constrained edge的helper function
		auto addConstrainedEdge = [&constrainedEdges](const PowerDiagramPoint_3& vert_1,
			const PowerDiagramPoint_3& vert_2) {
				constrainedEdges.emplace_back(
					Point_3(vert_1.x(), vert_1.y(), vert_1.z()),
					Point_3(vert_2.x(), vert_2.y(), vert_2.z()));
			};
		// 完整走一段等值线的helper function
		auto walkIsoLine = [&](int beg_isoVertTableIdx) {
			PowerDiagramPoint_3 end_isoVert;
			int cur_isoVertTableIdx = beg_isoVertTableIdx;
			std::unordered_set<int> isVisIsoVert;
			do
			{
				isVisIsoVert.insert(cur_isoVertTableIdx);

				const PowerDiagramPoint_3 iso_vert_1 = isoLineVerts[cur_isoVertTableIdx];
				for (const auto& idx : edgeTable.at(cur_isoVertTableIdx))
					if (!isVisIsoVert.count(idx)) cur_isoVertTableIdx = idx;  // TODO: 我们希望等值线上的每个点只有一个邻接点
				const PowerDiagramPoint_3 iso_vert_2 = isoLineVerts[cur_isoVertTableIdx];
				addConstrainedEdge(iso_vert_1, iso_vert_2); // push结果

				end_isoVert = iso_vert_2;
			} while (isoEdgePointToIdx.find(end_isoVert) == isoEdgePointToIdx.end()); // 直到等值线的某个点是边上的边界点为止
			return end_isoVert;
			};


		while (!bdSegQueue.empty())
		{
			VertIdxPair fro_bdSeg = bdSegQueue.front();
			bdSegQueue.pop();
			if (isVisSeg.count(fro_bdSeg)) continue;
			isVisSeg.insert(fro_bdSeg);

			// 1.2 Set boundary of current cell
			/*constrainedPoints.clear(); */constrainedEdges.clear();

			CELL_TYPE cellType = CELL_TYPE::EQUAL;
			std::unordered_set<int> cur_cellBDIsoVertIdxSet; // 存储当前cell所经历过的在边界上的等值线点，用于更新等值线的type

			int beg_cellBDPointIdx = fro_bdSeg.second;
			int cur_cellBDPointIdx = beg_cellBDPointIdx; // 始终存储的是边界点的索引，该索引为edgePWTVec这一vector中的索引

			//int last_cellIsoPointIdx = -1; // 防止重复走一条等值线
			do {
				PowerDiagramPoint_3 cur_cellBDPoint = edgePWTVec[cur_cellBDPointIdx].pos; // 存储当前cell正在处理的边界点

				// 如果边界点是某段等值线上的点
				if (isoEdgePointToIdx.count(cur_cellBDPoint)/* &&
					isoEdgePointToIdx.at(curBDPoint).second != last_cellIsoPointIdx*/)
				{
					// 设置当前cell所属的type
					if (cellType == CELL_TYPE::EQUAL && edgePWTVec[cur_cellBDPointIdx].type != CELL_TYPE::EQUAL)
						cellType = (CELL_TYPE)!int(edgePWTVec[cur_cellBDPointIdx].type);

					const int beg_isoVertIdx = cur_cellBDPointIdx; // 等值线的起始端点在edgePWTVec中的索引
					cur_cellBDPoint = walkIsoLine(isoEdgePointToIdx.at(cur_cellBDPoint)); // 获得在edgeTable中的索引并走完一条等值线，此时cur_cellBDPoint为等值线的末端点
					const int end_isoVertIdx = isoEdgePointToPWTIdx.at(cur_cellBDPoint); // 等值线的末端点在edgePWTVec中的索引

					cur_cellBDIsoVertIdxSet.insert(beg_isoVertIdx);
					cur_cellBDIsoVertIdxSet.insert(end_isoVertIdx);

					cur_cellBDPointIdx = end_isoVertIdx;
					//last_cellIsoPointIdx = cur_cellBDPointIdx;
				}
				else if (cellType == CELL_TYPE::EQUAL)
					cellType = edgePWTVec[cur_cellBDPointIdx].type; // 如果边界点是三角形的端点，则该cell的type由该端点确定

				// 走向下一条边界边，并将边界边的vis设为true
				int next_cellPointIdx = (cur_cellBDPointIdx + 1) % numEdgeCCWPoints;
				const auto next_seg = std::make_pair(cur_cellBDPointIdx, next_cellPointIdx);
				if (!isVisSeg.count(next_seg)) isVisSeg.insert(next_seg);

				// push结果
				const PowerDiagramPoint_3 bd_edge_vert_1 = edgePWTVec[cur_cellBDPointIdx].pos;
				const PowerDiagramPoint_3 bd_edge_vert_2 = edgePWTVec[next_cellPointIdx].pos;
				addConstrainedEdge(bd_edge_vert_1, bd_edge_vert_2);

				cur_cellBDPointIdx = next_cellPointIdx;
			} while (cur_cellBDPointIdx != beg_cellBDPointIdx);

			// 更新经过的所有等值线端点的type
			if (cellType == CELL_TYPE::EQUAL) std::cerr << "Error! Cell type is not be updated!\n";
			std::for_each(cur_cellBDIsoVertIdxSet.begin(), cur_cellBDIsoVertIdxSet.end(), [&](int idx) {
				edgePWTVec[idx].type = cellType;
				});

			std::vector<Vector3> verts; std::vector<int> vertIdx;
			int meshVertIdx = (cellType == CELL_TYPE::INSIDE) ? insideMeshVertIdx : outsideMeshVertIdx;

			/*if (faceIdx == 2173 && cellType == CELL_TYPE::INSIDE)
			{
				const auto aroundVerts = this->getFaceVec()[faceIdx]->aroundVerts();
				std::cout << "端点#0: " << aroundVerts[0]->pos.transpose() << std::endl;
				std::cout << "端点#1: " << aroundVerts[1]->pos.transpose() << std::endl;
				std::cout << "端点#2: " << aroundVerts[2]->pos.transpose() << std::endl;
			}*/

			for (int i = 0; i < constrainedEdges.size(); ++i)
			{
				const auto& vert = constrainedEdges[i].vertex(0);
				verts.emplace_back(vert.x(), vert.y(), vert.z());
				vertIdx.emplace_back(meshVertIdx + i);

				/*if (faceIdx == 2173 && cellType == CELL_TYPE::INSIDE)
				{
					std::cout << "vert: " << vert << "\nvertIdx: meshVertIdx + i" << std::endl;
					system("pause");
				}*/
			}
			/*for (auto it = constrainedEdges.begin(); it != constrainedEdges.end(); ++it)
			{
				const auto& vert = it->vertex(0);
				verts.emplace_back(vert.x(), vert.y(), vert.z());
				vertIdx.emplace_back(meshVertIdx + it - constrainedEdges.begin());
			}*/

			if (cellType == CELL_TYPE::INSIDE) insideMeshVertIdx += verts.size(), insideMesh.emplace_back(verts, vertIdx);
			else if (cellType == CELL_TYPE::OUTSIDE) outsideMeshVertIdx += verts.size(), outsideMesh.emplace_back(verts, vertIdx);

			//// 2. 对每个区域计算CDT
			//CDT_Cell cell(_bdEdges);
			//const auto& cdtForCurCell = cdtAdaptor.computeCDForBoundary(cell);

			//// 3. push结果
			//facetAllCDTs.emplace_back(cdtForCurCell);
		}

		return facetAllCDTs;
	}

	///////////////////////
	//  Algorithm  Core  //
	///////////////////////
	/**
	* @brief: 对采样点计算Isoline
	* @param out: 输出流，保存Isoline结果
	*/
	double MSCuttingModel::computeIsoline(std::ofstream& out)
	{
		int globalOutVertIdx = 1; // .obj的顶点下标从1开始
		std::map<PowerDiagramPoint_3, int> pdPointToOutIdx;
		std::unordered_map<int, std::unordered_map<int, std::vector<std::pair<PowerDiagramPoint_3, int>>>> edgeToPointsInNeighFace; // 保存边到其关联的两个face的所有存在该条边上点的映射

		// 保存当前面的PowerDiagram顶点的信息
		std::vector<std::vector<PowerDiagramPoint_3>> facetPDPoints(numMeshFaces);

		// 保存当前面的PowerDiagram边的信息(通过顶点的输出索引进行存储)
		std::vector<std::vector<std::pair<int, int>> > facetPDLines(numMeshFaces);
		std::vector<std::unordered_map<int, std::vector<int>> > edgeTable(numMeshFaces); // 保存该面power diagram点的邻边表

		/* Step1. 对所有采样面计算Power Diagram并筛选 */
		for (int i = 0; i < numMeshFaces; ++i)
		{
			computePDForFacet(i, globalOutVertIdx, facetPDPoints[i], facetPDLines[i], edgeTable[i], pdPointToOutIdx);
			//LOG::qpInfo("-- Compute coarse power diagram is finished!\n");
		}

		/* Step2. 对筛选后的粗糙等值线结果进行后处理 */
		double mae = .0;
		int numPDPoints = .0;
		for (int i = 0; i < numMeshFaces; ++i)
		{
			std::array<std::vector<PowerDiagramPoint_3>, 3> isoEdgePoints; // 保存在边上的(后处理过的)点
			std::map<PowerDiagramPoint_3, int> isoEdgePointToIdx; // 保存边上的(后处理过的)点到edgeTable索引的映射
			const std::vector<PowerDiagramPoint_3>& projCurFacetPDPoints =
				postProcessFacetPoints(i, facetPDPoints, pdPointToOutIdx, isoEdgePoints, isoEdgePointToIdx, edgeToPointsInNeighFace);
			//LOG::qpInfo("-- Post-processing is finished!\n");

			// 计算CDT
			//const std::vector<CellTriangle>& cdt_res = computeCDTForFacet(i, projCurFacetADPoints, edgeTable, isoEdgePoints, isoEdgePointToIdx);
			//LOG::qpInfo("-- Cutting is finished!\n");

			// 输出Power Diagram中的所有全局坐标点
			if (POST_PROCESSING_MAX_ITER <= 0)
			{
				for (const auto& pdPoint : facetPDPoints[i])
				{
					out << "v " << pdPoint.transpose() << std::endl;
					mae += std::fabs(scalarFunc.val(pdPoint));
				}
				numPDPoints += facetPDPoints[i].size();
			}
			else
			{
				for (const auto& pdPoint : projCurFacetPDPoints)
				{
					out << "v " << pdPoint.transpose() << std::endl;
					mae += std::fabs(scalarFunc.val(pdPoint));
				}
				numPDPoints += projCurFacetPDPoints.size();
			}

			// 输出Power Diagram中的所有边
			for (const auto& vertIdxPair : facetPDLines[i])
			{
				out << "l " << vertIdxPair.first << " " << vertIdxPair.second << std::endl;
			}
		}
		mae /= numPDPoints;

		/* Step3. 将相邻两个面的等值线连接起来 */
		// 将边上的点统一按逆时针排序
		struct DirSort {
			PowerDiagramPoint_3 edgeBaseVert;
			Vector3 edgeDir;
			DirSort(const PowerDiagramPoint_3& _edgeBaseVert, const Vector3& _edgeDir) :edgeBaseVert(_edgeBaseVert), edgeDir(_edgeDir) {}

			bool operator()(const std::pair<PowerDiagramPoint_3, int>& lhs, const std::pair<PowerDiagramPoint_3, int>& rhs)
			{
				return ((lhs.first - edgeBaseVert).dot(edgeDir) < (rhs.first - edgeBaseVert).dot(edgeDir));
			}
		};
		const auto& mEdgeVec = this->getMEdgeVec();
		std::vector<std::pair<std::pair<PowerDiagramPoint_3, int>, std::pair<PowerDiagramPoint_3, int>>> groupEdgePoint;
		for (auto& edge_face_point : edgeToPointsInNeighFace)
		{
			if (edge_face_point.second.size() < 2) continue;
			if (edge_face_point.second.size() != 2) { LOG::qpError("Non-manifold edge!"); continue; }

			const int mEdgeIdx = edge_face_point.first;
			const auto mEdgeVerts = mEdgeVec[mEdgeIdx]->verts();
			const Vector3 mEdgeDir = (mEdgeVerts.second->pos - mEdgeVerts.first->pos).normalized();
			const auto baseVert = mEdgeVerts.first->pos;

			DirSort dirSort(baseVert, mEdgeDir);

			int idx = 0;
			std::array<int, 2> edgeNeighFace;
			for (auto& face_point : edge_face_point.second)
			{
				edgeNeighFace[idx++] = face_point.first;
				std::sort(face_point.second.begin(), face_point.second.end(), dirSort);
			}

			for (int j = 0; j < edge_face_point.second.at(edgeNeighFace[0]).size(); ++j)
			{
				auto p1 = edge_face_point.second.at(edgeNeighFace[0])[j];
				auto p2 = edge_face_point.second.at(edgeNeighFace[1])[j];
				if ((p1.first - baseVert).dot(mEdgeDir) > (p2.first - baseVert).dot(mEdgeDir)) std::swap(p1, p2);
				//groupEdgePoint.emplace_back(p1, p2);
				out << "l " << p1.second << " " << p2.second << std::endl;
			}

			/*if (groupEdgePoint.size() >= 2)
				for (int i = 0; i < groupEdgePoint.size() - 1; ++i)
				{
					if (scalarFunc.grad(groupEdgePoint[i].second.first).normalized().
						isApprox(scalarFunc.grad(groupEdgePoint[i + 1].first.first).normalized()))
						out << "l " << groupEdgePoint[i].second.second << " " << groupEdgePoint[i + 1].first.second << std::endl;
				}*/
		}
		return mae;
	}

	///////////////////////
	//    Visualiztion   //
	///////////////////////
	/**
	* @brief: Visualization of sample points on each edge
	* @param out: 用于输出所有采样点的输出流
	*/
	void MSCuttingModel::outputSamplePoints(std::ofstream& out)
	{
		if (out)
		{
			for (const auto& samplePoint : samplePoints)
			{
				out << "v " << samplePoint.pos.transpose() << std::endl;
			}
		}
	}

	/**
	* @brief: Visualization of sample points on each edge
	* @param out_1: 用于输出inside mesh
	* @param out_2: 用于输出outside mesh
	*/
	void MSCuttingModel::outputCutMesh(std::ofstream& out_1, std::ofstream& out_2)
	{
		for (const auto& comp : insideMesh)
		{
			std::for_each(comp.pointVec.begin(), comp.pointVec.end(), [&out_1](Vector3 point) {
				out_1 << "v " << point.transpose() << "\n";
				});
			out_1 << "f ";
			std::for_each(comp.pointIdxVec.begin(), comp.pointIdxVec.end(), [&out_1](int pointIdx) {
				out_1 << pointIdx + 1 << " ";
				});
			out_1 << "\n";
		}

		for (const auto& comp : outsideMesh)
		{
			std::for_each(comp.pointVec.begin(), comp.pointVec.end(), [&out_2](Vector3 point) {
				out_2 << "v " << point.transpose() << "\n";
				});
			out_2 << "f ";
			std::for_each(comp.pointIdxVec.begin(), comp.pointIdxVec.end(), [&out_2](int pointIdx) {
				out_2 << pointIdx + 1 << " ";
				});
			out_2 << "\n";
		}
	}

	///////////////////////
	//    API for user   //
	///////////////////////
	/**
	* @brief: 暴露给用户调用的接口，用于运行整体算法
	* @param iter:				 设置后处理的迭代次数
	* @param isolineVisFile:     等值线输出结果所保存的文件位置
	* @param insideMeshVisFile:  内mesh输出结果所保存的文件位置
	* @param outsideMeshVisFile: 外mesh输出结果所保存的文件位置
	* @return: 算法运行成功/错误
	*/
	bool MSCuttingModel::launch(int iter,
		std::ofstream& mae_out,
		const std::string& isolineVisFile,
		const std::string& insideMeshVisFile,
		const std::string& outsideMeshVisFile)
	{
		/* 检查输出流 */
		str_util::checkDir(isolineVisFile);
		std::ofstream pd_vis_out(isolineVisFile);
		if (!pd_vis_out) { LOG::qpError("I/O: File ", isolineVisFile.c_str(), " could not be opened!"); return false; }

		/* 设置采样点输出路径 */
		const std::string sample_vis_file = str_util::concatFilePath(VIS_DIR, modelName, std::to_string(numSamplesPerEdge), "sample_points.obj");
		str_util::checkDir(sample_vis_file);
		std::ofstream sample_vis_out(sample_vis_file);
		if (!sample_vis_out) { LOG::qpError("I/O: File ", sample_vis_file.c_str(), " could not be opened!"); return false; }

		/* 采样 */
		samplePointPerEdge();
		/* 输出采样点 */
		LOG::qpInfo("Output Sample Points to ", std::quoted(sample_vis_file), " ...");
		outputSamplePoints(sample_vis_out);
		sample_vis_out.close();

		/* 计算等值线 */
		double mae = computeIsoline(pd_vis_out);
		mae_out << "iter: " << iter << " " << mae << "\n";
		/* 输出采样点 */
		LOG::qpInfo("Output Isoline to ", std::quoted(isolineVisFile), " ...");
		pd_vis_out.close();

		/* 计算Cutting结果 */
		bool cut = (!insideMeshVisFile.empty() && !outsideMeshVisFile.empty()) ? true : false;
		if (!cut) LOG::qpInfo("Cutting process is disabled.");
		else {
			str_util::checkDir(insideMeshVisFile);
			std::ofstream inside_vis_out(insideMeshVisFile);
			if (!inside_vis_out) { LOG::qpError("I/O: File ", insideMeshVisFile.c_str(), " could not be opened!"); return false; }

			str_util::checkDir(outsideMeshVisFile);
			std::ofstream outside_vis_out(outsideMeshVisFile);
			if (!outside_vis_out) { LOG::qpError("I/O: File ", outsideMeshVisFile.c_str(), " could not be opened!"); return false; }

			outputCutMesh(inside_vis_out, outside_vis_out);
			inside_vis_out.close(); inside_vis_out.close();
		}

		return true;
	}

#if ENABLE_TEST
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
		}

		return true;
	}

	bool MSCuttingModel::testComputeADForFacet(int facetIdx, std::ofstream& out)
	{
		samplePointPerEdge();

		// 保存当前面的PowerDiagram顶点的信息
		std::vector<PowerDiagramPoint_3> facetADPoints;

		// 保存当前面的PowerDiagram边的信息(通过顶点的输出索引进行存储)
		std::vector<std::pair<int, int>> facetADLines;

		int globalOutVertIdx = 1; // .obj的顶点下标从1开始

		std::unordered_map<int, std::vector<int>> edgeTable; // 保存该面power diagram点的邻边表
		std::map<PowerDiagramPoint_3, int> pdPointToOutIdx;
		int numADLines = computePDForFacet(facetIdx, globalOutVertIdx, facetADPoints, facetADLines, edgeTable, pdPointToOutIdx);

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
#endif

} // namespace core
NAMESPACE_END(mscut)