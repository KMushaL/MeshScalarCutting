#pragma once
#include "Config.hpp"
#include "geometry/PolyMesh.hpp"

NAMESPACE_BEGIN(mscut)
namespace unit_test
{
	using namespace detail;
	using namespace geometry;

	inline bool testMeshNorm(const std::string& mesh_file)
	{
		PolyMesh polyMesh(mesh_file);
		polyMesh.meshNorm();

		const std::string norm_mesh_file = str_util::concatFilePath(MODEL_DIR, polyMesh.modelName + (std::string)"_norm.obj");
		polyMesh.writeToOBJ(norm_mesh_file);
	}

} // namespace unit_test
NAMESPACE_END(mscut)