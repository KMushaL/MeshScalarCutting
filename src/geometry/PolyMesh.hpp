#pragma once
#include "Config.hpp"
#include "Geometry.hpp"
#include "BaseMesh.hpp"
#include <string>
#include <vector>
#include <memory/MemoryPool.h>

NAMESPACE_BEGIN(mscut)
namespace geometry
{
	/* 为了你方便以后调试，先不用模板 */
	//template<bool isNormalized = false>
	class PolyMesh {
	public:
		/* Member types */
		using Point3 = HEVert::Point3;

	private:
		/* Basic data in matrix style(from libigl) */
		MatrixX vertMat;
		MatrixXi faceMat;
		MatrixX vertNormalMat;
		MatrixX faceNormalMat;

	private:
		/* Pointer(point to address of data in MemoryPool) is stored in vector */
		std::vector<HEVert*> vertVec;
		std::vector<HEFace*> faceVec;
		std::vector<MEdge*> mEdgeVec;
		std::vector<HEdge*> hEdgeVec;

		/* MemoryPool is where the data is actually stored */
		MemoryPool<HEVert> vertPool;
		MemoryPool<HEFace> facePool;
		MemoryPool<MEdge> mEdgePool;
		MemoryPool<HEdge> hEdgePool;

	protected:
		int numMeshVerts;
		int numMeshEdges;
		int numMeshFaces;

	protected:
		/* Axis-aligned bounding box */
		double boxScaleFactor;
		double diagonalLengthOfBBox;
		AABox<Vector3> modelBoundingBox;

	public:
		std::string modelName;
		std::string textureName;

	public:
		/* Constructor and Destructor */
		PolyMesh() noexcept = default;

		PolyMesh(const std::string& file, const bool is2Unifrom, const double _scaleFactor) noexcept :boxScaleFactor(_scaleFactor) {
			readMesh(file);

			//if constexpr (isNormalized) meshNorm();

			if (is2Unifrom) setBoundingBox<true>();
			else setBoundingBox<false>();

			modelBoundingBox.output(str_util::concatFilePath(VIS_DIR, modelName, (std::string)"modelBoundingBox.obj"));
			LOG::qpInfo("Diagonal length of the bounding-box = ", diagonalLengthOfBBox);

			constructMesh();
		}

		PolyMesh(const std::string& file, const bool is2Unifrom, const double _scaleFactor,
			const std::string& norm_out_file) noexcept :boxScaleFactor(_scaleFactor) {
			readMesh(file);

			meshNorm(norm_out_file);

			//if constexpr (isNormalized) meshNorm();

			if (is2Unifrom) setBoundingBox<true>();
			else setBoundingBox<false>();

			modelBoundingBox.output(str_util::concatFilePath(VIS_DIR, modelName, (std::string)"modelBoundingBox.obj"));
			LOG::qpInfo("Diagonal length of the bounding-box = ", diagonalLengthOfBBox);

			constructMesh();
		}

		PolyMesh(const std::string& file, const bool is2Unifrom) noexcept :PolyMesh(file, is2Unifrom, 1) {}

		explicit PolyMesh(const std::string& file) noexcept :PolyMesh(file, false, 1) {}

		explicit PolyMesh(const std::string& file, const std::string& norm_out_file) noexcept :PolyMesh(file, false, 1, norm_out_file) {}

		virtual ~PolyMesh() noexcept = default;

	public:
		/* Getter */
		size_t numVertices() const noexcept { return vertVec.size(); }

		size_t numFaces() const noexcept { return faceVec.size(); }

		size_t numEdges() const noexcept { return mEdgeVec.size(); }

		size_t numHEdges() const noexcept { return hEdgeVec.size(); }

		MatrixX& getVertMat() noexcept { return vertMat; }

		const MatrixX& getVertMat() const noexcept { return vertMat; }

		MatrixXi& getFaceMat() noexcept { return faceMat; }

		const MatrixXi& getFaceMat() const noexcept { return faceMat; }

		std::vector<HEVert*>& getVertVec() noexcept { return vertVec; }

		const std::vector<HEVert*>& getVertVec() const noexcept { return vertVec; }

		std::vector<HEFace*>& getFaceVec() noexcept { return faceVec; }

		const std::vector<HEFace*>& getFaceVec() const noexcept { return faceVec; }

		std::vector<MEdge*>& getMEdgeVec() noexcept { return mEdgeVec; }

		const std::vector<MEdge*>& getMEdgeVec() const noexcept { return mEdgeVec; }

		std::vector<HEdge*>& getHEdgeVec() noexcept { return hEdgeVec; }

		const std::vector<HEdge*>& getHEdgeVec() const noexcept { return hEdgeVec; }

		HEVert* vert(size_t id) noexcept { return (id < numVertices() ? vertVec[id] : nullptr); }

		const HEVert* vert(size_t id) const noexcept { return (id < numVertices() ? vertVec[id] : nullptr); }

		MEdge* edge(size_t id) noexcept { return (id < numEdges() ? mEdgeVec[id] : nullptr); }

		const MEdge* edge(size_t id) const noexcept { return (id < numEdges() ? mEdgeVec[id] : nullptr); }

		HEdge* halfEdge(size_t id) noexcept { return (id < numHEdges() ? hEdgeVec[id] : nullptr); }

		const HEdge* halfEdge(size_t id) const noexcept { return (id < numHEdges() ? hEdgeVec[id] : nullptr); }

	private:
		/* Memory allocate and deallocate */
		void reserveMemory(size_t numVerts, size_t numFaces);

		void clear();

		HEVert* newVertex();

		void deleteVertexInMem(HEVert* vert);

		HEFace* newFace();

		void deleteFaceInMem(HEFace* face);

		MEdge* newEdge();

		void deleteEdgeInMem(MEdge* mEdge);

		MEdge* newEdge(HEVert* v1, HEVert* v2);

		HEdge* newHalfEdge();

		void deleteHalfEdgeInMem(HEdge* hEdge);

	private:
		/* Data manipulator */
		HEVert* addVertex(const Point3& point);

		// Add an edge and two associated half-edges
		MEdge* addEdge(HEVert* v1, HEVert* v2);

		// Add a face by given vertices, it will not check added face is duplicate
		HEFace* addFace(std::vector<HEVert*>& aroundVerts);

		// Add a face by given half-edges, which are needed to be connected
		HEFace* addFace(const std::vector<HEdge*>& halfEdges);

		void deleteVertex(HEVert* vert);

		void deleteFace(HEFace* face);

		// Delete edge and associated half-edges
		void deleteEdge(MEdge* edge);

		void deleteMultiVerts(const std::vector<HEVert*>& verts);

		void deleteMultiFaces(const std::vector<HEFace*>& faces);

		void deleteMultiEdges(const std::vector<MEdge*>& mEdges);

		// Construct the whole half-edge structure from input mesh
		void constructMesh();

	public:
		/* Data manipulator(api for user) */
		size_t deleteIsolatedVerts();

		// Delete isolated edges and associated half-edges
		size_t deleteIsolatedEdges();

	private:
		/* Set (uniform)bounding-box */
		template<bool is2Uniform>
		void setBoundingBox();

	private:
		/* Methods for normalizing the mesh */
		//template<typename = std::enable_if_t<isNormalized, Matrix<Scalar, 4, 4>>>
		Matrix<Scalar, 4, 4> calcTransformMatrix();

		//template<typename = std::enable_if_t<isNormalized>>
		void meshInverseNorm();

	public:
		//template<typename = std::enable_if_t<isNormalized>>
		void meshNorm();

		void meshNorm(const std::string& out_file);

	protected:
		/* Methods for transforming the mesh */
		void zoomModel();

		void transformModel(const float& _scaleFactor);

	public:
		/* Data properties check */
		// only for 'HEVert, MEdge, HEdge'
		template<class U>
		bool isBoundary(U* data) const {
			return (data == nullptr || data->isBoundary());
		}

		// only for 'HEVert, HEFace, MEdge'
		template<class U>
		bool isIsolated(U* data) const {
			return (data == nullptr || data->isIsolated());
		}

	public:
		/* Basic(Traversal) operations for data */
		/* Get vertex's neighbors */
		std::vector<HEVert*> vertOneRingVerts(HEVert* vert) const;

		std::vector<MEdge*> vertOneRingEdges(HEVert* vert) const;

		std::vector<HEdge*> vertOneRingHEdges(HEVert* vert) const;

		std::vector<HEFace*> vertOneRingFaces(HEVert* vert) const;

		/* Get face's vertices, edges, half-edges and adjacent faces */
		std::vector<HEVert*> faceVerts(HEFace* face) const;

		std::vector<MEdge*> faceMEdges(HEFace* face) const;

		std::vector<HEdge*> faceHEdges(HEFace* face) const;

		std::vector<HEFace*> faceAdjacentFaces(HEFace* face) const;

		/* Get edge's adjacent faces */
		std::vector<HEFace*> edgeAdjacentFaces(MEdge* edge) const;

		/* Get edge associated half-edge by id */
		HEdge* edgeHalfEdge(MEdge* mEdge, short id) const;

		/* Get edge by its two ending vertices */
		MEdge* edgeBetween(HEVert* v1, HEVert* v2) const;

		/* Get half-edge by its two ending vertices(direction is important) */
		HEdge* halfEdgeBetween(HEVert* v1, HEVert* v2) const;

		/* Get half-edge from which a vertex in face */
		HEdge* faceHalfEdgeFromVert(HEFace* face, HEVert* v) const;

		/* Check whether two vertices are connected */
		bool isConnected(HEVert* v1, HEVert* v2) const;

		/* Check whether two vertices are connected */
		bool isFaceContainsVertex(HEFace* face, HEVert* v) const;

	public:
		/* Basic topology check */
		bool isTriMesh() const;

		bool isManifold() const;

	public:
		/* APIs for triangle mesh */

	public:
		/* Low level apis */
		void printInfo() const;

	public:
		/* I/O */
		/* Input */
		void readMesh(const std::string& in_file);

		/* Output */
		void writeToOBJ(const std::string& out_file);

		void writeToOFF(const std::string& out_file);

		void writeToPLY(const std::string& out_file);
	};

} // namespace geometry
NAMESPACE_END(mscut)

#include "PolyMesh.inl"