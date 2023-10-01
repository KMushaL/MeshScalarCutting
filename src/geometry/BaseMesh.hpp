//
// Created by lei on 23-9-12.
//
#pragma once
#include "Config.hpp"
#include "detail/BasicDataType.hpp"
#include <array>
#include <utility>
#include <functional>

NAMESPACE_BEGIN(mscut)
namespace geometry
{
	using namespace detail;

	class HEVert;

	class HEdge;

	class HEFace;

	class MEdge;

	/**
	 * Class of 'vertex' for polygonal mesh in 'HalfEdge' structure
	 * index begins with 0
	 */
	class HEVert {
	public:
		using Point3 = Vector3;

	public:
		int index;
		Point3 pos;
		Vector3 normal;

		std::array<float, 4> color;
		std::array<float, 3> textCoord;

	private:
		HEdge* he; // half-edge beginning from the vertex

	public:
		HEVert() : index(-1), pos({ 0, 0, 0 }), normal({ 0, 0, 0 }),
			color({ 0, 0, 0, 0 }), textCoord({ 0, 0, 0 }), he(nullptr) {}

		explicit HEVert(Point3 _pos) : index(-1), pos(std::move(_pos)), normal({ 0, 0, 0 }),
			color({ 0, 0, 0, 0 }), textCoord({ 0, 0, 0 }), he(nullptr) {}

		HEVert(Scalar x, Scalar y, Scalar z) : index(-1), pos({ x, y, z }), normal({ 0, 0, 0 }),
			color({ 0, 0, 0, 0 }), textCoord({ 0, 0, 0 }), he(nullptr) {}

		virtual ~HEVert() noexcept = default;

	public:
		/* Getter */
		HEdge* halfEdge() noexcept { return he; }

		const HEdge* halfEdge() const noexcept { return he; }

		Scalar x() const noexcept { return pos.x(); }

		Scalar y() const noexcept { return pos.y(); }

		Scalar z() const noexcept { return pos.z(); }

		Scalar nx() const noexcept { return normal.x(); }

		Scalar ny() const noexcept { return normal.y(); }

		Scalar nz() const noexcept { return normal.z(); }

		/* Setter */
		void setHEdge(HEdge* const& _he) noexcept { he = _he; }

		void setColor(float r, float g, float b) { color = { r, g, b, 0 }; }

		void setColor(float r, float g, float b, float a) { color = { r, g, b, a }; }

		void setTexture(float u, float v) { textCoord = { u, v, 0 }; }

		void setTexture(float u, float v, float w) { textCoord = { u, v, w }; }

	private:
		void updateNormal(const std::function<Vector3(const std::vector<HEFace*>&)>&);

		std::vector<HEFace*> getOneRingFaces() const;

	public:
		void update();

		bool isBoundary() const;

		bool isIsolated() const;

		/* Ensure that each vertex stores a boundary half-edge
		 * (used for construction and to prevent errors during traversal operations) */
		void adjustBoundaryHEdge();
	};

	/**
	 * Class of 'face' for polygonal mesh in 'HalfEdge' structure
	 */
	class HEFace {
		using Point3 = typename HEVert::Point3;
	private:
		enum class ConvexType {
			CONVEX,
			NO_CONVEX,
			INVALID
		};

	public:
		int index;

		int numVerts; // number of vertices around the face

		Vector3 normal;

		ConvexType isConvex;

		std::array<float, 4> color;

	private:
		HEdge* he;

		std::vector<HEVert*> verts;

	public:
		HEFace() : index(-1), numVerts(0), normal({ 0, 0, 0 }), isConvex(ConvexType::INVALID), color({ 0, 0, 0, 0 }),
			he(nullptr) {}

		explicit HEFace(HEdge* _he)
			: index(-1), numVerts(0), normal({ 0, 0, 0 }), isConvex(ConvexType::INVALID), color({ 0, 0, 0, 0 }),
			he(_he) {}

		virtual ~HEFace() noexcept = default;

	public:
		/* Getter */
		HEdge* halfEdge() noexcept { return he; }

		const HEdge* halfEdge() const noexcept { return he; }

		std::vector<HEVert*> const& aroundVerts() noexcept { return verts; }

		/* Setter */
		void setFroVertex(std::vector<HEVert*> const& _verts) noexcept { verts = _verts; }

		void setHEdge(HEdge* const& _he) noexcept { he = _he; }

		void setAroundVerts(const std::vector<HEVert*>& _verts) noexcept { verts = _verts; }

		void setColor(float r, float g, float b) { color = { r, g, b, 0 }; }

		void setColor(float r, float g, float b, float a) { color = { r, g, b, a }; }

	private:
		void updateVerts();

		void updateNormal();

		void updateConvex();

	public:
		bool isIsolated() const;

		Point3 getFaceCenter() const;

		int isInside(const Point3& point) const;

		std::vector<Scalar> getBaryCoord();

		void update() {
			updateVerts();
			updateNormal();
			// updateConvex();
		}
	};

	/**
	 * class of 'edge' for polygonal mesh
	 * the index is different from 'index' in HEdge
	 */
	class MEdge {
	public:
		using VertPair = std::pair<HEVert*, HEVert*>;

	public:
		int index;

	private:
		VertPair vert_pair; // two vertices of the edge

		HEdge* he; // one of half-edges associated with the edge

	public:
		MEdge() : index(-1), vert_pair(std::make_pair(nullptr, nullptr)), he(nullptr) {}

		explicit MEdge(VertPair _verts) : index(-1), vert_pair(_verts), he(nullptr) {}

		MEdge(HEVert* _v1, HEVert* _v2) : MEdge(std::make_pair(_v1, _v2)) {}

		MEdge(VertPair _verts, HEdge* _he) : index(-1), vert_pair(_verts),
			he(_he) {}

		MEdge(HEVert* _v1, HEVert* _v2, HEdge* _he) : MEdge(std::make_pair(_v1, _v2), _he) {}

		virtual ~MEdge() noexcept = default;

	public:
		/* Getter */
		VertPair verts() noexcept { return vert_pair; }

		const VertPair verts() const noexcept { return vert_pair; }

		HEVert* const firstVertex() noexcept { return vert_pair.first; }

		const HEVert* firstVertex() const noexcept { return vert_pair.first; }

		HEVert* const secondVertex() noexcept { return vert_pair.second; }

		const HEVert* secondVertex() const noexcept { return vert_pair.second; }

		HEdge* const halfEdge() noexcept { return he; }

		const HEdge* halfEdge() const noexcept { return he; }

		/* Setter */
		void setVertex(VertPair const& _vert_pair) noexcept { vert_pair = _vert_pair; }

		void setFirstVertex(HEVert* const& _vert) noexcept { vert_pair.first = _vert; }

		void setSecondVertex(HEVert* const& _vert) noexcept { vert_pair.second = _vert; }

		void setHEdge(HEdge* const& _he) noexcept { he = _he; }

	public:
		bool isBoundary() const;

		bool isIsolated() const;

		void updateVerts();

		typename HEVert::Point3 center() {
			updateVerts();
			return (vert_pair.first->pos + vert_pair.second->pos) * 0.5;
		}

		Scalar length() {
			updateVerts();
			return (vert_pair.first->pos - vert_pair.second->pos).norm();
		}
	};

	/**
	 * Class of 'edge' for triangle mesh in 'HalfEdge' structure
	 */
	class HEdge {
	public:
		int index;

		std::array<float, 3> textCoord;

	private:
		HEVert* fro_vert; // vertex at the beginning of the half-edge

		MEdge* m_edge;   // associated edge

		HEdge* e_prev;   // previous half-edge
		HEdge* e_next;   // next half-edge around the face
		HEdge* e_oppo;   // opposite half-edge of the half-edge

		HEFace* face;    // face the half-edge borders

	public:
		HEdge() : index(-1), fro_vert(nullptr), m_edge(nullptr),
			e_prev(nullptr), e_next(nullptr), e_oppo(nullptr), face(nullptr), textCoord({ 0, 0, 0 }) {}

		HEdge(HEVert* _fro_vert, MEdge* _m_edge, HEdge* _e_prev, HEdge* _e_next, HEdge* _e_oppo,
			HEFace* _face)
			: index(-1), fro_vert(_fro_vert), m_edge(_m_edge),
			e_prev(_e_prev), e_next(_e_next), e_oppo(_e_oppo), face(_face), textCoord({ 0, 0, 0 }) {}

		virtual ~HEdge() noexcept = default;

	public:
		/* Getter */
		HEVert* froVertex() noexcept { return fro_vert; }

		const HEVert* froVertex() const noexcept { return fro_vert; }

		HEVert* nextVertex() noexcept { return e_next->fro_vert; }

		const HEVert* nextVertex() const noexcept { return e_next->fro_vert; }

		MEdge* mEdge() noexcept { return m_edge; }

		const MEdge* mEdge() const noexcept { return m_edge; }

		HEdge* prevHEdge() noexcept { return e_prev; }

		const HEdge* prevHEdge() const noexcept { return e_prev; }

		HEdge* nextHEdge() noexcept { return e_next; }

		const HEdge* nextHEdge() const noexcept { return e_next; }

		HEdge* oppoHEdge() noexcept { return e_oppo; }

		const HEdge* oppoHEdge() const noexcept { return e_oppo; }

		HEFace* boundFace() noexcept { return face; }

		const HEFace* boundFace() const noexcept { return face; }

		// Next half-edge which shared common vertex
		HEdge* rotateNext() noexcept { return oppoHEdge()->e_next; }

		const HEdge* rotateNext() const noexcept { return oppoHEdge()->e_next; }

		// Previous half-edge which shared common vertex
		HEdge* rotatePrev() noexcept { return prevHEdge()->e_oppo; }

		const HEdge* rotatePrev() const noexcept { return prevHEdge()->e_oppo; }

		/* Setter */
		void setFroVertex(HEVert* const& _vert) noexcept { fro_vert = _vert; }

		//            void setFroVertex(HEVert *_vert) noexcept { fro_vert = std::move(_vert); }
		void setMEdge(MEdge* const& _m_edge) noexcept { m_edge = _m_edge; }

		void setPrevEdge(HEdge* const& _prev) noexcept { e_prev = _prev; }

		void setNextEdge(HEdge* const& _next) noexcept { e_next = _next; }

		void setOppoEdge(HEdge* const& _oppo) noexcept { e_oppo = _oppo; }

		void setFace(HEFace* const& _face) noexcept { face = _face; }

		void setTexture(float u, float v) { textCoord = { u, v, 0 }; }

		void setTexture(float u, float v, float w) { textCoord = { u, v, w }; }

	public:
		bool isBoundary() const { return face == nullptr; }

		Vector3 dir() const {
			return froVertex()->pos - nextVertex()->pos;
		}

		Vector3 oppoDir() const {
			return -dir();
		}

		Vector3 tangent() const {
			return (froVertex()->pos - nextVertex()->pos).normalized();
		}
	};

} // namespace geometry
NAMESPACE_END(mscut)