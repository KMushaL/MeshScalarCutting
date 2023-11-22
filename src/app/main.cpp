#include <iostream>
#include <CLI/CLI.hpp>
#include "utils/String.hpp"
#include "core/test/CoreTest.hpp"
#include "geometry/test/GeoTest.hpp"

using namespace mscut;

using ScalarFunc = core::MSCuttingModel::ScalarFunc;

// level-set定义: 函数式和梯度
Scalar val(const Eigen::Vector3d& p)
{
	return p.x() * p.x() + p.y() * p.y() - 1;
	//return p.y() * p.y() - p.x() * p.x() + p.x() * p.x() * p.x();
	//return p.x() - p.y();
	//return (p.x() + 1) * (p.x() + 1) + (p.y()) * (p.y()) + (p.z()) * (p.z()) - 2.25;
	//return (p.x()) * (p.x()) + (p.y()) * (p.y()) + (p.z()) * (p.z()) - 1;
	//return (2.0 / 3) * std::pow(p.x() * p.x(), 1.0 / 3) + (2.0 / 3) * std::pow(p.y() * p.y(), 1.0 / 3) - 1;
}
Eigen::Vector3d grad(const Eigen::Vector3d& p)
{
	return  2 * Eigen::Vector3d(p.x(), p.y(), 0);
	//return Eigen::Vector3d(3 * p.x() * p.x() - 2 * p.x(), 2 * p.y(), 0);
	//return Eigen::Vector3d(1, -1, 0);
	//return 2 * Eigen::Vector3d((p.x() + 1), p.y(), p.z());
	//return 2 * Eigen::Vector3d(p.x(), p.y(), p.z());
	//return (4.0 / 9) * Eigen::Vector3d(std::pow(std::fabs(p.x()), -1.0 / 3) * (p.x() < 0 ? -1 : 1), std::pow(std::fabs(p.y()), -1.0 / 3) * (p.y() < 0 ? -1 : 1), 0);
}

int main(int argc, char** argv)
{
	// parse input arguments
	CLI::App app("Mesh Scalar Cutting");
	//app.add_flag("-h,--help", "Print configuration and exit.");

	//std::string modelArg = R"(test2d.obj)";
	std::string modelArg = R"(square3.obj)";
	app.add_option("-f,--file", modelArg,
		"Input model's name with extension.")/*->required()*/;

	int m_numSamplePoints = 8;
	int n_numSamplePoints = 4;
	app.add_option("-m", m_numSamplePoints,
		"Specify the number of sample points at each edge of input model.");
	app.add_option("-n", n_numSamplePoints,
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

	ScalarFunc scalarFunc = { val, grad };

	// 传入奇点
	std::vector<Vector3> singulars;

	singulars.emplace_back(Vector3(-1.83713, 0, 0));
	singulars.emplace_back(Vector3(1.83713, 0, 0));
	singulars.emplace_back(Vector3(0, -1.83713, 0));
	singulars.emplace_back(Vector3(0, 1.83713, 0));
	singulars.emplace_back(Vector3(0, 0, 0));
	core::MSCuttingModel mscModel(mesh_file, scalarFunc, singulars, m_numSamplePoints, n_numSamplePoints);
	/*const std::string norm_mesh_file(str_util::concatFilePath(VIS_DIR, modelName, modelName + "_norm.obj"));
	core::MSCuttingModel mscModel(mesh_file, norm_mesh_file, scalarFunc, numSamplePoints);*/

	/*const std::string ad_vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, "apollonius_diagram.obj");
	mscModel.launch(ad_vis_file);*/
	const std::string pd_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, std::to_string(n_numSamplePoints), "isoline.obj"
	);
	mscModel.launch(pd_vis_file);

	/*const std::string inside_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, std::to_string(numSamplePoints), "insideMesh.obj"
	);
	const std::string outside_vis_file = str_util::concatFilePath(
		VIS_DIR, mscModel.modelName, std::to_string(numSamplePoints), "outsideMesh.obj"
	);
	mscModel.launch(pd_vis_file, inside_vis_file, outside_vis_file);*/

	/*unit_test::testMeshNorm(mesh_file);

	ScalarFunc scalarFunc = { val, grad };
	unit_test::testSamplingPoint(mesh_file, scalarFunc, 2);

	unit_test::testLocalGlobalTransoform(mesh_file, numSamplePoints);

	unit_test::testComputeAD(mesh_file, numSamplePoints);
	unit_test::testComputeCDT();*/

	return 0;
}