//
// Created by lei on 23-9-13.
//
#include "PolyMesh.hpp"
#include "utils/Log.hpp"
#include "utils/Donut.hpp"
#include <omp.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>

NAMESPACE_BEGIN(mscut)
namespace geometry
{
	/* Memory allocate and deallocate */
	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::reserveMemory(size_t numVerts, size_t numFaces) {
		vertVec.reserve(numVerts);
		//            vertPool.allocate(numVerts);

		faceVec.reserve(numFaces);
		//            facePool.allocate(numFaces);

					// at least larger than number of vertices for triangle mesh
		mEdgeVec.reserve(3 * numVerts);
		//            mEdgePool.allocate(3 * numVerts);

		hEdgeVec.reserve(6 * numVerts);
		//            hEdgePool.allocate(6 * numVerts);
	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::clear() {
		for (const auto& vert : vertVec)
			//                vertPool.destroy(vert);
			vertPool.deallocate(vert);
		vertVec.clear();
		vertVec.shrink_to_fit();

		for (const auto& face : faceVec)
			facePool.deallocate(face);
		faceVec.clear();
		faceVec.shrink_to_fit();

		for (const auto& mEdge : mEdgeVec)
			mEdgePool.deallocate(mEdge);
		mEdgeVec.clear();
		mEdgeVec.shrink_to_fit();

		for (const auto& hEdge : hEdgeVec)
			hEdgePool.deallocate(hEdge);
		hEdgeVec.clear();
		hEdgeVec.shrink_to_fit();
	}

	//template<bool isNormalized>
	HEVert* PolyMesh/*<isNormalized>*/::newVertex() {
		HEVert* _vert = vertPool.allocate();
		vertPool.construct(_vert, HEVert());
		//            new(_vert) HEVert();
		_vert->index = vertVec.size();
		vertVec.emplace_back(_vert);
		return _vert;
	}

	//template<bool isNormalized>
	HEFace* PolyMesh/*<isNormalized>*/::newFace() {
		HEFace* _face = facePool.allocate();
		facePool.construct(_face, HEFace());
		//            new(_face) HEFace();
		_face->index = faceVec.size();
		faceVec.emplace_back(_face);
		return _face;
	}

	//template<bool isNormalized>
	MEdge* PolyMesh/*<isNormalized>*/::newEdge() {
		MEdge* _mEdge = mEdgePool.allocate();
		mEdgePool.construct(_mEdge, MEdge());
		//            new(_mEdge) MEdge();
		_mEdge->index = mEdgeVec.size();
		mEdgeVec.emplace_back(_mEdge);
		return _mEdge;
	}

	//template<bool isNormalized>
	MEdge* PolyMesh/*<isNormalized>*/::newEdge(HEVert* v1, HEVert* v2) {
		MEdge* _mEdge = mEdgePool.allocate();
		mEdgePool.construct(_mEdge, MEdge(v1, v2));
		//            new(_mEdge) MEdge(v1, v2);
		_mEdge->index = mEdgeVec.size();
		mEdgeVec.emplace_back(_mEdge);
		return _mEdge;
	}

	//template<bool isNormalized>
	HEdge* PolyMesh/*<isNormalized>*/::newHalfEdge() {
		HEdge* _hEdge = hEdgePool.allocate();
		hEdgePool.construct(_hEdge, HEdge());
		//            new(_hEdge) HEdge();
		_hEdge->index = hEdgeVec.size();
		hEdgeVec.emplace_back(_hEdge);
		return _hEdge;
	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteVertexInMem(HEVert* vert) {
		assert(vert->index < vertVec.size());
		assert(vert == vertVec[vert->index]);

		int idx = vert->index;
		vertPool.destroy(vert);
		vertPool.deallocate(vert);

		// Place the last element at the position of 'vert' in 'vertVec'
		if (idx != vertVec.size() - 1) {
			auto* last = vertVec.back();
			last->index = idx; // update the index
			vertVec[idx] = last; // replace 'vert'
		}
		// Delete the last element in 'vertVec'
		vertVec.pop_back();
	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteFaceInMem(HEFace* face) {
		assert(face->index < faceVec.size());
		assert(face == faceVec[face->index]);

		int idx = face->index;
		facePool.destroy(face);
		facePool.deallocate(face);

		if (idx != faceVec.size() - 1) {
			auto* last = faceVec.back();
			last->index = idx;
			faceVec[idx] = last;
		}
		faceVec.pop_back();
	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteEdgeInMem(MEdge* mEdge) {
		assert(mEdge->index < mEdgeVec.size());
		assert(mEdge == mEdgeVec[mEdge->index]);

		int idx = mEdge->index;
		mEdgePool.destroy(mEdge);
		mEdgePool.deallocate(mEdge);

		if (idx != mEdgeVec.size() - 1) {
			auto* last = mEdgeVec.back();
			last->index = idx;
			mEdgeVec[idx] = last;
		}
		mEdgeVec.pop_back();
	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteHalfEdgeInMem(HEdge* hEdge) {
		assert(hEdge->index < hEdgeVec.size());
		assert(hEdge == hEdgeVec[hEdge->index]);

		int idx = hEdge->index;
		hEdgePool.destroy(hEdge);
		hEdgePool.deallocate(hEdge);

		if (idx != hEdgeVec.size() - 1) {
			auto* last = hEdgeVec.back();
			last->index = idx;
			hEdgeVec[idx] = last;
		}
		hEdgeVec.pop_back();
	}

	/* Data manipulator */
	//template<bool isNormalized>
	HEVert* PolyMesh/*<isNormalized>*/::addVertex(const Point3& point) {
		HEVert* _vert = newVertex();
		_vert->pos = point;
		return _vert;
	}

	//template<bool isNormalized>
	MEdge* PolyMesh/*<isNormalized>*/::addEdge(HEVert* v1, HEVert* v2) {
		if (v1 == nullptr || v2 == nullptr) return nullptr;

		MEdge* mEdge = edgeBetween(v1, v2);
		// already added
		if (mEdge != nullptr) return mEdge;

		mEdge = newEdge(v1, v2);
		HEdge* he_1 = newHalfEdge();
		HEdge* he_2 = newHalfEdge();

		mEdge->setHEdge(he_1);

		he_1->setNextEdge(he_2);
		he_1->setPrevEdge(he_2);
		he_1->setOppoEdge(he_2);
		he_1->setFroVertex(v1);
		he_1->setMEdge(mEdge);

		he_2->setNextEdge(he_1);
		he_2->setPrevEdge(he_1);
		he_2->setOppoEdge(he_1);
		he_2->setFroVertex(v2);
		he_2->setMEdge(mEdge);

		return mEdge;
	}

	//template<bool isNormalized>
	HEFace* PolyMesh/*<isNormalized>*/::addFace(std::vector<HEVert*>& aroundVerts) {
		if (aroundVerts.empty()) return nullptr;

		size_t numVerts = aroundVerts.size();

		// The value is true when the i'th edge is newly constructed
		std::vector<bool> isNewEdge(numVerts, false);
		// The value is true when the i'th vertex need adjust its half-edge
		// Each vertex should store one of its boundary half-edge if exist
		std::vector<bool> isNeedAdjust(numVerts, false);
		// Stores all inner-side half-edges of the new new_face
		std::vector<HEdge*> halfEdgesVec;

		// Firstly, create a new new_face
		HEFace* new_face = newFace();

		// Secondly, add all vertices of the new_face
		// he_i represents inner-side half-edge of the new new_face,
		// whereas he_j represents outer-side half-edge
		HEdge* he_i = nullptr, * oppo_he_i = nullptr;
		for (size_t i = 0; i < numVerts; ++i) {
			size_t j = (i + 1) % numVerts;
			auto vert_i = aroundVerts[i];
			auto vert_j = aroundVerts[j];

			// Check whether the edge is inserted before
			MEdge* edge_ij = edgeBetween(vert_i, vert_j);
			if (edge_ij == nullptr) {
				edge_ij = addEdge(vert_i, vert_j);
				isNewEdge[i] = true;
			}
			if (edge_ij->firstVertex() == vert_i) {
				he_i = edge_ij->halfEdge();
				oppo_he_i = he_i->oppoHEdge();
			}
			else {
				he_i = edge_ij->halfEdge()->oppoHEdge();
				oppo_he_i = he_i->oppoHEdge();
			}

			if (new_face->halfEdge() == nullptr)
				new_face->setHEdge(he_i);

			he_i->setFace(new_face);
			halfEdgesVec.emplace_back(he_i);
		}

		// Thirdly, for two adjacent edges of the new_face that have both been inserted
		// we need to adjust the global half-edge connection relationship.
		// (because the newly added new_face causes a 'complete(interconnected) patch',
		// but before that, there may exist patches that are not adjacent to it(isolated patches),
		// which are needed to integrate the connection of half-edge with the 'complete(interconnected) patch'.)
		HEdge* inner_he_i = nullptr, * inner_he_i_next = nullptr;
		HEdge* outer_he_i = nullptr, * outer_he_i_next = nullptr;
		for (size_t i = 0; i < numVerts; ++i) {
			size_t j = (i + 1) % numVerts;
			if (!isNewEdge[i] && !isNewEdge[j]) {
				inner_he_i = halfEdgesVec[i];
				inner_he_i_next = halfEdgesVec[j];

				if (inner_he_i->nextHEdge() != inner_he_i_next) { // exist isolated patches
					outer_he_i = inner_he_i->oppoHEdge();
					outer_he_i_next = inner_he_i_next->oppoHEdge();

					// Find the boundary half-edges of the complete patch
					// The isolated patches separate the two half-edges on two sides
					// i.e., the two boundary half-edges sandwich the isolated patches
					HEdge* boundary_he = nullptr, * boundary_he_next = nullptr; // stores boundary half-edge of the complete patch
					boundary_he = outer_he_i_next;
					do {
						boundary_he = boundary_he->nextHEdge()->oppoHEdge();
					} while (!boundary_he->isBoundary());
					boundary_he_next = boundary_he->nextHEdge();

					// Get the beginning and ending half-edges of the isolated patches
					HEdge* isoPatchBeg = inner_he_i->nextHEdge();
					HEdge* isoPatchEnd = inner_he_i_next->prevHEdge();

					// Adjust the relationship
					boundary_he->setNextEdge(isoPatchBeg);
					isoPatchBeg->setPrevEdge(boundary_he);

					isoPatchEnd->setNextEdge(boundary_he_next);
					boundary_he_next->setPrevEdge(isoPatchEnd);

					inner_he_i->setNextEdge(inner_he_i_next);
					inner_he_i_next->setPrevEdge(inner_he_i);
				}
			}
		}

		// Fourthly, process the connection relationship between two adjacent half-edges, one(or each) of which is a newly added half-edge
		// and also update the variable 'isNeedAdjust' of the edges processed in the third step
		for (size_t i = 0; i < numVerts; ++i) {
			size_t j = (i + 1) % numVerts;

			HEVert* vert_ij = aroundVerts[j]; // The common vertex of two adjacent inner-side half-edges(i and j) of the newly new_face
			inner_he_i = halfEdgesVec[i];
			inner_he_i_next = halfEdgesVec[j];

			unsigned short idx = 0;
			if (isNewEdge[i]) idx |= 1;
			if (isNewEdge[j]) idx |= 2;

			if (idx) {
				outer_he_i = inner_he_i->oppoHEdge();
				outer_he_i_next = inner_he_i_next->oppoHEdge();

				HEdge* old_half_edge_prev = nullptr;
				HEdge* old_half_edge_next = nullptr;
				switch (idx) {
				case 1: // MEdge of the 'inner_he_i' is new but MEdge of the 'inner_he_i_next' is old
					// Currently, 'inner_he_i_next' is the outer-side half-edge of the old edge
					old_half_edge_prev = inner_he_i_next->prevHEdge();
					old_half_edge_prev->setNextEdge(outer_he_i);
					outer_he_i->setPrevEdge(old_half_edge_prev);
					// Update the half-edge of the common vertex to guarantee a boundary half-edge
					vert_ij->setHEdge(outer_he_i);
					break;
				case 2: // 'inner_he_i' is old but 'inner_he_i_next' is new
					old_half_edge_next = inner_he_i->nextHEdge();
					old_half_edge_next->setPrevEdge(outer_he_i_next);
					outer_he_i_next->setNextEdge(old_half_edge_next);
					vert_ij->setHEdge(old_half_edge_next);
					break;
				case 3: // both are new
					if (vert_ij->halfEdge() == nullptr) {
						vert_ij->setHEdge(outer_he_i);
						// reverse of inner-side
						outer_he_i_next->setNextEdge(outer_he_i);
						outer_he_i->setPrevEdge(outer_he_i_next);
					}
					else { // patches shared a common vertex(vert_ij) in this case
						// we need to adjust the connection related to these patches
						HEdge* half_edge_vert_ij = vert_ij->halfEdge();
						HEdge* prev_half_edge_vert_ij = half_edge_vert_ij->prevHEdge();

						prev_half_edge_vert_ij->setNextEdge(outer_he_i);
						outer_he_i->setPrevEdge(prev_half_edge_vert_ij);

						outer_he_i_next->setNextEdge(half_edge_vert_ij);
						half_edge_vert_ij->setPrevEdge(outer_he_i_next);
					}
					break;
				}
			}
			else {
				isNeedAdjust[j] = true;
			}

			inner_he_i->setNextEdge(inner_he_i_next);
			inner_he_i_next->setPrevEdge(inner_he_i);
		}

		// Finally, perform an adjust vertex operation on the edge if 'isNeedAdjust' is true
		for (size_t i = 0; i < numVerts; ++i) {
			if (isNeedAdjust[i]) aroundVerts[i]->adjustBoundaryHEdge();
		}

		// new_face->setAroundVerts(aroundVerts);
		new_face->update();
		return new_face;
	}

	//template<bool isNormalized>
	HEFace* PolyMesh/*<isNormalized>*/::addFace(const std::vector<HEdge*>& halfEdges) {
		if (halfEdges.empty()) return nullptr;

	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteVertex(HEVert* vert) {
		if (vert == nullptr) return;

		if (isIsolated<HEVert>(vert))
			deleteVertexInMem(vert);
		else {
			const auto oneRingEdges = vertOneRingEdges(vert);
			for (auto edge : oneRingEdges)
				deleteEdge(edge);
			deleteVertexInMem(vert);
		}
	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteFace(HEFace* face) {

	}

	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::deleteEdge(MEdge* edge) {

	}

	/*template<bool isNormalized>*/
	void PolyMesh/*<isNormalized>*/::constructMesh() {
		for (size_t i = 0; i < numMeshVerts; ++i) {
			addVertex(vertMat.row(i));
		}

		for (size_t i = 0; i < numMeshFaces; ++i) {
			std::vector<HEVert*> aroundVerts;
			for (size_t j = 0; j < faceMat.row(i).size(); ++j) {
				size_t vertIdx = faceMat(i, j);
				if (vertIdx == -1) break;
				aroundVerts.emplace_back(vertVec[vertIdx]);
			}
			donut::removeDupicates(aroundVerts);
			addFace(aroundVerts);
		}

		numMeshEdges = numEdges();
	}

	/* Methods for normalizing the mesh */
	//template<bool isNormalized>
	//template<typename>
	Matrix<Scalar, 4, 4> PolyMesh/*<isNormalized>*/::calcTransformMatrix()
	{
		Vector3 boxMin = vertMat.colwise().minCoeff();
		Vector3 boxMax = vertMat.colwise().maxCoeff();

		// Get the target solveRes (along the largest dimension)
		Scalar scale = boxMax[0] - boxMin[0];
		Scalar minScale = scale;
		for (int d = 1; d < 3; d++)
		{
			scale = std::max<Scalar>(scale, boxMax[d] - boxMin[d]);
			minScale = std::min<Scalar>(scale, boxMax[d] - boxMin[d]);
		}
		scale *= scaleFactor;
		Vector3 center = 0.5 * boxMax + 0.5 * boxMin;

		Matrix<Scalar, 4, 4> zoomMatrix = Matrix<Scalar, 4, 4>::Identity();
		Matrix<Scalar, 4, 4> transMatrix = Matrix<Scalar, 4, 4>::Identity();

		//// transform至[0,1]
		//for (int i = 0; i < 3; i++)
		//	center[i] -= scale / 2;
		//for (int i = 0; i < 3; i++)
		//{
		//	zoomMatrix(i, i) = 1. / scale;
		//	transMatrix(3, i) = -center[i];
		//}

		// transform至[-1,1]
		for (int i = 0; i < 3; i++)
		{
			zoomMatrix(i, i) = 2. / scale;
			transMatrix(3, i) = -center[i];
		}
		return zoomMatrix * transMatrix;
	}

	//template<bool isNormalized>
	//template<typename>
	void PolyMesh/*<isNormalized>*/::meshNorm(const std::string& out_file)
	{
		auto transMat = calcTransformMatrix();
#pragma omp parallel for
		for (int i = 0; i < numMeshVerts; ++i)
		{
			vertMat.row(i) += transMat.block(3, 0, 1, 3);
			vertMat.row(i) = vertMat.row(i) * transMat.block(0, 0, 3, 3);
		}

		if (!out_file.empty()) writeToOBJ(out_file);
	}

	void PolyMesh::zoomModel()
	{
		auto transMat = calcTransformMatrix();
#pragma omp parallel for
		for (int i = 0; i < vertMat.rows(); ++i) {
			vertMat.row(i) = vertMat.row(i) * transMat.block(0, 0, 3, 3);
		}
	}

	//template<bool isNormalized>
	//template<typename>
	void PolyMesh/*<isNormalized>*/::meshInverseNorm()
	{
		Matrix<Scalar, 4, 4> transMat = calcTransformMatrix();
		Matrix<Scalar, 3, 3> inverseTrans = transMat.block(0, 0, 3, 3).inverse();
#pragma omp parallel for
		for (int i = 0; i < numMeshVerts; ++i)
		{
			vertMat.row(i) = vertMat.row(i) * inverseTrans;
			vertMat.row(i) -= transMat.block(3, 0, 1, 3);
		}
	}

	/* Get vertex's neighbors */
	//template<bool isNormalized>
	std::vector<HEVert*> PolyMesh/*<isNormalized>*/::vertOneRingVerts(HEVert* vert) const {
		std::vector<HEVert*> oneRingVerts;
		if (vert == nullptr || vert->isIsolated()) return oneRingVerts;

		HEdge* he_begin = vert->halfEdge();
		HEdge* _he = he_begin;
		do {
			oneRingVerts.emplace_back(_he->nextVertex());
			_he = _he->rotateNext();
		} while (_he != he_begin);

		return oneRingVerts;
	}

	//template<bool isNormalized>
	std::vector<MEdge*> PolyMesh/*<isNormalized>*/::vertOneRingEdges(HEVert* vert) const {
		std::vector<MEdge*> oneRingEdges;
		if (vert == nullptr || vert->isIsolated()) return oneRingEdges;

		HEdge* he_begin = vert->halfEdge();
		HEdge* _he = he_begin;
		do {
			oneRingEdges.emplace_back(_he->mEdge());
			_he = _he->rotateNext();
		} while (_he != he_begin);

		return oneRingEdges;
	}

	//template<bool isNormalized>
	std::vector<HEdge*> PolyMesh/*<isNormalized>*/::vertOneRingHEdges(HEVert* vert) const {
		std::vector<HEdge*> oneRingHEdges;
		if (vert == nullptr || vert->isIsolated()) return oneRingHEdges;

		HEdge* he_begin = vert->halfEdge();
		HEdge* _he = he_begin;
		do {
			oneRingHEdges.emplace_back(_he);
			if (_he->oppoHEdge() != nullptr)
				oneRingHEdges.emplace_back(_he->oppoHEdge());
			_he = _he->rotateNext();
		} while (_he != he_begin);

		return oneRingHEdges;
	}

	//template<bool isNormalized>
	std::vector<HEFace*> PolyMesh/*<isNormalized>*/::vertOneRingFaces(HEVert* vert) const {
		std::vector<HEFace*> oneRingFaces;
		if (vert == nullptr || vert->isIsolated()) return oneRingFaces;

		HEdge* he_begin = vert->halfEdge();
		HEdge* _he = he_begin;
		do {
			if (!_he->isBoundary())
				oneRingFaces.emplace_back(_he->boundFace());
			_he = _he->rotateNext();
		} while (_he != he_begin);

		return oneRingFaces;
	}

	/* Get face's vertices, edges, half-edges and adjacent faces */
	//template<bool isNormalized>
	std::vector<HEVert*> PolyMesh/*<isNormalized>*/::faceVerts(HEFace* face) const {
		return face->aroundVerts();
	}

	//template<bool isNormalized>
	std::vector<MEdge*> PolyMesh/*<isNormalized>*/::faceMEdges(HEFace* face) const {
		std::vector<MEdge*> f_mEdges;
		if (face == nullptr) return f_mEdges;

		HEdge* he_begin = face->halfEdge();
		HEdge* _he = he_begin;
		do {
			f_mEdges.emplace_back(_he->mEdge());
			_he = _he->nextHEdge();
		} while (_he != he_begin);

		return f_mEdges;
	}

	//template<bool isNormalized>
	std::vector<HEdge*> PolyMesh/*<isNormalized>*/::faceHEdges(HEFace* face) const {
		std::vector<HEdge*> f_hEdges;
		if (face == nullptr) return f_hEdges;

		HEdge* he_begin = face->halfEdge();
		HEdge* _he = he_begin;
		do {
			f_hEdges.emplace_back(_he);
			if (_he->oppoHEdge() != nullptr)
				f_hEdges.emplace_back(_he->oppoHEdge());
			_he = _he->nextHEdge();
		} while (_he != he_begin);

		return f_hEdges;
	}

	//template<bool isNormalized>
	std::vector<HEFace*> PolyMesh/*<isNormalized>*/::faceAdjacentFaces(HEFace* face) const {
		std::vector<HEFace*> f_adjFaces;
		if (face == nullptr || face->isIsolated()) return f_adjFaces;

		HEdge* he_begin = face->halfEdge();
		HEdge* _he = he_begin;
		do {
			if (_he->oppoHEdge()->boundFace() != nullptr)
				f_adjFaces.emplace_back(_he->oppoHEdge()->boundFace());
			_he = _he->nextHEdge();
		} while (_he != he_begin);

		return f_adjFaces;
	}

	/* Get edge's adjacent faces */
	//template<bool isNormalized>
	std::vector<HEFace*> PolyMesh/*<isNormalized>*/::edgeAdjacentFaces(MEdge* edge) const {
		std::vector<HEFace*> e_adjFaces;
		if (edge == nullptr || edge->isIsolated()) return e_adjFaces;

		e_adjFaces.emplace_back(edge->halfEdge()->boundFace());
		e_adjFaces.emplace_back(edge->halfEdge()->oppoHEdge()->boundFace());

		return e_adjFaces;
	}

	/* Get edge associated half-edge by id */
	//template<bool isNormalized>
	HEdge* PolyMesh/*<isNormalized>*/::edgeHalfEdge(MEdge* mEdge, short id) const {
		switch (id) {
		case 0:
			return mEdge->halfEdge();
		case 1:
			return mEdge->halfEdge()->oppoHEdge();
		default:
			return nullptr;
		}
	}

	/* Get edge by its two ending vertices */
	//template<bool isNormalized>
	MEdge* PolyMesh/*<isNormalized>*/::edgeBetween(HEVert* v1, HEVert* v2) const {
		// if (v1 == nullptr || v2 == nullptr) return nullptr;
		// if (v1->halfEdge() == nullptr) return nullptr;
		// HEdge *he_begin = v1->halfEdge();
		// HEdge *_he = he_begin;
		// do {
		//     if (_he->nextVertex() == v2) return _he->mEdge();
		//     _he = _he->rotateNext();
		// } while (_he != he_begin);
		// return nullptr;

		HEdge* hEdge = halfEdgeBetween(v1, v2);
		return hEdge ? hEdge->mEdge() : nullptr;
	}

	/* Get half-edge by its two ending vertices(direction is important) */
	//template<bool isNormalized>
	HEdge* PolyMesh/*<isNormalized>*/::halfEdgeBetween(HEVert* v1, HEVert* v2) const {
		if (v1 == nullptr || v2 == nullptr) return nullptr;
		if (v1->halfEdge() == nullptr && v2->halfEdge() == nullptr) return nullptr;

		HEdge* he_begin;
		HEdge* _he;
		if (v1->halfEdge() != nullptr) {
			he_begin = v1->halfEdge();
			_he = he_begin;

			do {
				if (_he->nextVertex() == v2) return _he;
				_he = _he->rotateNext();
			} while (_he != he_begin);
		}
		else {
			he_begin = v2->halfEdge();
			_he = he_begin;

			do {
				if (_he->nextVertex() == v1) return _he;
				_he = _he->rotateNext();
			} while (_he != he_begin);
		}

		return nullptr;
	}

	/* Get half-edge from which a vertex in face */
	//template<bool isNormalized>
	HEdge* PolyMesh/*<isNormalized>*/::faceHalfEdgeFromVert(HEFace* face, HEVert* v) const {
		const auto hEdgesOfFace = faceHEdges(face);
		for (const auto& hEdge : hEdgesOfFace) {
			if (hEdge->froVertex() == v) return hEdge;
		}
		return nullptr;
	}

	/* Check whether two vertices are connected */
	//template<bool isNormalized>
	bool PolyMesh/*<isNormalized>*/::isConnected(HEVert* v1, HEVert* v2) const {
		return (halfEdgeBetween(v1, v2) != nullptr);
	}

	/* Check whether two vertices are connected */
	//template<bool isNormalized>
	bool PolyMesh/*<isNormalized>*/::isFaceContainsVertex(HEFace* face, HEVert* v) const {
		return std::find(face->aroundVerts().begin(), face->aroundVerts().end(), v) !=
			std::end(face->aroundVerts());
	}

	/* Basic topology check */
	//template<bool isNormalized>
	bool PolyMesh/*<isNormalized>*/::isTriMesh() const {
		return std::all_of(faceVec.cbegin(), faceVec.cend(), [](HEFace* face) {
			return face->numVerts == 3;
			}
		);
	}

	/* Low level apis */
	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::printInfo() const {
		LOG::qpInfo("Model: ", modelName);
		LOG::qpInfo("-- The number of vertices = ", numVertices());
		LOG::qpInfo("-- The number of faces = ", numFaces());
		LOG::qpInfo("-- The number of edges = ", numEdges());
		LOG::qpInfo("-- The number of half-edges = ", numHEdges());
	}

	/* I/O */
	/* Input */
	//template<bool isNormalized>
	void PolyMesh/*<isNormalized>*/::readMesh(const std::string& in_file) {
		// TODO: replace with a reader for any polygonal mesh
		igl::read_triangle_mesh(in_file, vertMat, faceMat);

		modelName = str_util::getFileName(in_file);

		numMeshVerts = vertMat.rows();
		numMeshFaces = faceMat.rows();
	}

	void PolyMesh::writeToOBJ(const std::string& out_file)
	{
		str_util::checkDir(out_file);

		igl::writeOBJ(out_file, vertMat, faceMat);
	}

} // namespace geometry
NAMESPACE_END(mscut)