//
// Created by lei on 23-9-13.
//
#include "BaseMesh.hpp"
#include "detail/Math.hpp"
#include <iostream>

NAMESPACE_BEGIN(mscut)
namespace geometry
{
	/**
	 * HEVert
	 */
	void HEVert::updateNormal(const std::function<Vector3(const std::vector<HEFace*>&)>& weightFunc) {
		std::vector<HEFace*> faces = getOneRingFaces();
		normal = weightFunc(faces);
	}

	void HEVert::update() {
		auto constantWeightFunc = [](const std::vector<HEFace*>& faces) -> Vector3 {
			Vector3 vertNormal = { 0, 0, 0 };
			std::for_each(faces.begin(), faces.end(), [&vertNormal](HEFace* face) {
				vertNormal += face->normal;
				});
			return vertNormal;
		};
		updateNormal(constantWeightFunc);
	}

	std::vector<HEFace*> HEVert::getOneRingFaces() const {
		std::vector<HEFace*> oneRingFaces;
		if (this->isIsolated()) return oneRingFaces;

		HEdge* _he = he;
		do {
			if (!_he->isBoundary())
				oneRingFaces.emplace_back(_he->boundFace());
			_he = _he->rotateNext();
		} while (_he != he);

		return oneRingFaces;
	}

	bool HEVert::isBoundary() const {
		if (this->isIsolated()) return true;

		HEdge* _he = he;
		do {
			if (_he->isBoundary()) return true;
			_he = _he->rotateNext();
		} while (_he != he);

		return false;
	}

	bool HEVert::isIsolated() const {
		return he == nullptr;
	}

	void HEVert::adjustBoundaryHEdge() {
		if (isIsolated()) return;

		HEdge* _he = he;
		do {
			if (_he->isBoundary()) {
				setHEdge(_he);
				break;
			}
			_he = _he->rotateNext();
		} while (_he != he);
	}

	/**
	 * HEFace
	 */
	void HEFace::updateVerts() {
		HEdge* t_he = he;
		do {
			++numVerts;
			// std::cout << "vertex: " << t_he->froVertex()->pos.transpose() << std::endl;
			verts.emplace_back(t_he->froVertex());
			t_he = t_he->nextHEdge();
		} while (t_he != he);
	}

	void HEFace::updateNormal() {
		const Vector3 dir_1 = he->dir();
		const Vector3 dir_2 = he->prevHEdge()->oppoDir();
		normal = dir_1.cross(dir_2);
		normal.normalize();
		// printf("normal: (%lf, %lf, %lf)\n", normal.x(), normal.y(), normal.z());
	}

	void HEFace::updateConvex() {
		if (numVerts < 3) {
			isConvex = ConvexType::INVALID;
			return;
		}
		HEdge* t_he = he;
		bool firstFlag = (t_he->dir().cross(t_he->nextHEdge()->dir())).dot(normal) > 0;
		bool flag;
		do {
			t_he = t_he->nextHEdge();
			flag = (t_he->dir().cross(t_he->nextHEdge()->dir())).dot(normal) > 0;
		} while (flag == firstFlag || t_he->nextHEdge() != he);

		isConvex = (flag == firstFlag) ? ConvexType::CONVEX : ConvexType::NO_CONVEX;
	}

	bool HEFace::isIsolated() const {
		if (he == nullptr) return true;

		HEdge* _he = he;
		HEdge* _he_oppo;
		do {
			_he_oppo = _he->oppoHEdge();
			if (!_he_oppo->isBoundary()) return false;
			_he = _he->nextHEdge();
		} while (_he != he);

		return true;
	}

	HEFace::Point3 HEFace::getFaceCenter() const {
		Point3 center = { 0, 0, 0 };
		std::for_each(verts.begin(), verts.end(), [&center](HEVert* heVert) {
			center += heVert->pos;
			});
		return center / numVerts;
	}

	int HEFace::isInside(const Point3& point) const {
		if (isConvex == ConvexType::CONVEX) {
			Scalar val = he->dir().cross(point - he->froVertex()->pos).dot(normal);
			return dcmp(val, GEO_EPSILON);
		}
		else if (isConvex == ConvexType::NO_CONVEX) {
			// TODO: ray method
		}
		else
			return -2;
	}

	std::vector<Scalar> HEFace::getBaryCoord() {
		return std::vector<Scalar>();
	}

	/**
	 * MEdge
	 */
	bool MEdge::isBoundary() const {
		if (he == nullptr) return true;

		HEdge* he_oppo = he->oppoHEdge();
		if (he_oppo == nullptr) return true;

		return (he->isBoundary() || he_oppo->isBoundary());
	}

	bool MEdge::isIsolated() const {
		if (he == nullptr) return true;

		HEdge* he_oppo = he->oppoHEdge();
		if (he_oppo == nullptr) return true;

		return (he->isBoundary() && he_oppo->isBoundary());
	}

	void MEdge::updateVerts() {
		vert_pair.first = he->froVertex();
		vert_pair.second = he->nextVertex();
	}

} // namespace geometry
NAMESPACE_END(mscut)
