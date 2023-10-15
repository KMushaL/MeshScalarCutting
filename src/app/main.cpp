#include <iostream>
#include "utils/String.hpp"
#include "core/test/CoreTest.hpp"
#include "geometry/test/GeoTest.hpp"

using namespace mscut;

using ScalarFunc = core::MSCuttingModel::ScalarFunc;

Scalar val(const Eigen::Vector3d& p)
{
	return p.x() * p.x() + p.y() * p.y() + p.z() * p.z() - 1;
}

Eigen::Vector3d grad(const Eigen::Vector3d& p)
{
	return 2 * p;
}

int main()
{
	const std::string mesh_file(MODEL_DIR DELIMITER R"(bunny.obj)");
	/*const static int numSamplePoints = 0;

	ScalarFunc scalarFunc = { val, grad };
	core::MSCuttingModel mscModel(mesh_file, scalarFunc, numSamplePoints);
	mscModel.meshNorm();

	const std::string ad_vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"apollonius_diagram.obj");
	mscModel.launch(ad_vis_file);*/

	unit_test::testMeshNorm(mesh_file);
	/*unit_test::testSamplingPoint(mesh_file, 2);

	unit_test::testLocalGlobalTransoform(mesh_file, numSamplePoints);

	unit_test::testComputeAD(mesh_file, numSamplePoints);*/

	return 0;
}