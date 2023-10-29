#include <iostream>
#include <CLI/CLI.hpp>
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

int main(int argc, char** argv)
{
	// parse input arguments
	CLI::App app("Mesh Scalar Cutting");
	//app.add_flag("-h,--help", "Print configuration and exit.");

	std::string modelArg = "bunny.obj";
	app.add_option("-f,--file", modelArg,
		"Input model's name with extension.");

	int numSamplePoints = 20;
	app.add_option("-n,--number", numSamplePoints,
		"Specify the number of sample points at each edge of input model.");

	try {
		argv = app.ensure_utf8(argv);
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError& e) {
		// 输出帮助信息
		return app.exit(e);
	}
	//std::cout << "modelArg = " << modelArg << std::endl;
	//printf("numSamplePoints = %d\n", numSamplePoints);

	const std::string mesh_file(MODEL_DIR + R"(\)" + modelArg);
	const std::string modelName = str_util::getFileName(mesh_file);

	const std::string norm_mesh_file(str_util::concatFilePath(VIS_DIR, modelName, modelName + "_norm.obj"));
	ScalarFunc scalarFunc = { val, grad };
	core::MSCuttingModel mscModel(mesh_file, norm_mesh_file, scalarFunc, numSamplePoints);

	/*const std::string ad_vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, "apollonius_diagram.obj");
	mscModel.launch(ad_vis_file);*/
	const std::string pd_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, std::to_string(numSamplePoints), "isoline.obj"
	);
	const std::string inside_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, std::to_string(numSamplePoints), "insideMesh.obj"
	);
	const std::string outside_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, std::to_string(numSamplePoints), "outsideMesh.obj"
	);
	mscModel.launch(pd_vis_file, inside_vis_file, outside_vis_file);

	/*unit_test::testMeshNorm(mesh_file);

	ScalarFunc scalarFunc = { val, grad };
	unit_test::testSamplingPoint(mesh_file, scalarFunc, 2);

	unit_test::testLocalGlobalTransoform(mesh_file, numSamplePoints);

	unit_test::testComputeAD(mesh_file, numSamplePoints);
	unit_test::testComputeCDT();*/

	return 0;
}