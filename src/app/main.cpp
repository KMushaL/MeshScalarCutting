#include <iostream>
#include "utils/String.hpp"
#include "core/test/CoreTest.hpp"
#include "geometry/test/GeoTest.hpp"

using namespace mscut;

using ScalarFunc = core::MSCuttingModel::ScalarFunc;

Scalar val(const Eigen::Vector3d& p)
{
	return (p.x() + 1) * (p.x() + 1) + (p.y()) * (p.y()) + (p.z()) * (p.z()) - 2.25;
}

Eigen::Vector3d grad(const Eigen::Vector3d& p)
{
	return 2 * Eigen::Vector3d((p.x() + 1), p.y(), p.z());
}

int main()
{
	const std::string mesh_file(MODEL_DIR DELIMITER R"(bunny.obj)");
	const std::string modelName = str_util::getFileName("\\", mesh_file);
	const static int numSamplePoints = 5;

	const std::string norm_mesh_file(str_util::concatFilePath(VIS_DIR, modelName, modelName + (std::string)"_norm.obj"));
	ScalarFunc scalarFunc = { val, grad };
	core::MSCuttingModel mscModel(mesh_file, norm_mesh_file, scalarFunc, numSamplePoints);

	/*const std::string ad_vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"apollonius_diagram.obj");
	mscModel.launch(ad_vis_file);*/
	const std::string pd_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, (std::string)"power_diagram_" + std::to_string(numSamplePoints) + (std::string)".obj"
	);
	mscModel.launch(pd_vis_file);

	/*
	unit_test::testMeshNorm(mesh_file);

	ScalarFunc scalarFunc = { val, grad };
	unit_test::testSamplingPoint(mesh_file, scalarFunc, 2);

	unit_test::testLocalGlobalTransoform(mesh_file, numSamplePoints);

	unit_test::testComputeAD(mesh_file, numSamplePoints);*/

	return 0;
}