#include <iostream>
#include <CLI/CLI.hpp>
#include "utils/String.hpp"
#include "core/MSCuttingModel.hpp"

using namespace std;
using namespace mscut;

using ScalarFunc = core::MSCuttingModel::ScalarFunc;

const std::string funcName = "torus";
// level-set定义: 函数式和梯度
static Scalar val(const Eigen::Vector3d& p)
{
	//return p.x() * p.x() + p.y() * p.y() - 1;
	//return p.y() * p.y() - p.x() * p.x() + p.x() * p.x() * p.x();
	//return p.x() - p.y();
	//return (p.x() + 1) * (p.x() + 1) + (p.y()) * (p.y()) + (p.z()) * (p.z()) - 2.25;
	//return (p.x()) * (p.x()) + (p.y()) * (p.y()) + (p.z()) * (p.z()) - 1;
	/*double x = p.x() - 1;
	return (2.0 / 3) * std::pow(x * x, 1.0 / 3) + (2.0 / 3) * std::pow(p.y() * p.y(), 1.0 / 3) - 1;*/

	double x = p.x();
	double y = p.y();
	double z = p.z();
	//return std::pow(x * x + z * z, 0.5) - std::pow(y * y, 0.5);							//yuanzhuixian1
	//return std::pow(x * x + y * y, 0.5) - std::pow(z * z, 0.5);							//yuanzhuixian2
	//return std::pow((std::pow(8 * x * x + 8 * z * z, 0.5) - 2), 2) + 8 * y * y - 2;		//torus1
	//return std::pow((std::pow(8 * x * x + 8 * y * y, 0.5) - 2), 2) + 8 * z * z - 2;		//torus2
	//return std::pow((x / 2), 4) + std::pow((y / 2), 4) + std::pow((z / 2), 4) - 0.03;		//tuoqiuti
	//return std::pow(4 * x * x + 4 * y * y + 4 * z * z, 2) - 3 * (4 * x * x - 4 * y * y) - 0.8;		//xingqiuti1
	//return std::pow(4 * x * x + 4 * y * y + 4 * z * z, 2) - 3 * (4 * y * y - 4 * z * z) - 0.8;		//xingqiuti2
	//return (15 * x * x + 15 * y * y + 15 * z * z + 1) * (15 * x * x + 15 * y * y + 15 * z * z + 1) - 16 * (15 * x * x + 15 * y * y) - 1;		//jiaoChaHuanMian1
	//return std::pow((15 * x * x + 15 * y * y + 15 * z * z + 1), 2) - 16 * (15 * y * y + 15 * z * z) - 1;		//jiaoChaHuanMian2
	//return std::pow(x * x + y * y + z * z - 1, 2) - z * z - 2 * x - 1.5; // penut

	return (4 * x * x + 4 * y * y + 4 * z * z + 4 * y - 1) * (std::pow((4 * x * x + 4 * y * y + 4 * z * z - 4 * y - 1), 2) - 32 * z * z) +
		64 * x * z * (4 * x * x + 4 * y * y + 4 * z * z - 4 * y - 1);  // klein bottle - 适用于hole/bettle/carp_muscle model

	//return std::pow((0.5 * x * x + 0.5 * z * z), 3) - 0.25 * x * x * z * z; // QUATREFOIL CURVE - 适用于hole model
	//return std::pow((4 * x * x + 4 * y * y), 3) - 16 * x * x * y * y; // QUATREFOIL CURVE - 适用于bettle model
	//return std::pow((x * x + y * y), 3) - x * x * y * y; // QUATREFOIL CURVE - 适用于carp_muscle model

	//return -2 * y - 4 * x * z + 8 * x * x * y - 16 * x * x * z + 8 * y * y * y - 16 * y * y * z + 16 * y * z * z; // mobius - 适用于hole/bettle model
	//return -y - 2 * x * z + x * x * y - 2 * x * x * z + y * y * y - 2 * y * y * z + y * z * z; // mobius - 适用于carp_muscle model
}
static Eigen::Vector3d grad(const Eigen::Vector3d& p)
{
	//return  2 * Eigen::Vector3d(p.x(), p.y(), 0);
	//return Eigen::Vector3d(3 * p.x() * p.x() - 2 * p.x(), 2 * p.y(), 0);
	//return Eigen::Vector3d(1, -1, 0);
	//return 2 * Eigen::Vector3d((p.x() + 1), p.y(), p.z());
	//return 2 * Eigen::Vector3d(p.x(), p.y(), p.z());
	/*double x = p.x() - 1;
	return (4.0 / 9) * Eigen::Vector3d(std::pow(std::fabs(x), -1.0 / 3) * (x < 0 ? -1 : 1), std::pow(std::fabs(p.y()), -1.0 / 3) * (p.y() < 0 ? -1 : 1), 0);*/

	double x = p.x();
	double y = p.y();
	double z = p.z();
	//return Eigen::Vector3d(x / sqrt(x * x + z * z), -y / sqrt(y * y), z / sqrt(x * x + z * z));	//yuanzhuixian1
	//return Eigen::Vector3d(x / sqrt(x * x + y * y), -y / sqrt(x * x + y * y), z / sqrt(z * z));	//yuanzhuixian2

	//return Eigen::Vector3d(16 * x * (1 - 1 / std::pow((2 * x * x + 2 * z * z), 0.5)), 16 * y, 16 * z * (1 - 1 / std::pow((2 * x * x + 2 * z * z), 0.5)));		//torus1
	//return Eigen::Vector3d(8 * x * (2 - sqrt(2) / sqrt(x * x + y * y)), 8 * y * (2 - sqrt(2) / sqrt(x * x + y * y)), 16 * z);		//torus2
	//return Eigen::Vector3d(x * x * x / 4, y * y * y / 4, z * z * z / 4);															//tuoqiuti
	//return 8 * Eigen::Vector3d(x * (-3 + 8 * x * x + 8 * y * y + 8 * z * z), y * (3 + 8 * x * x + 8 * y * y + 8 * z * z), 8 * z * (x * x + y * y + z * z));	//xingqiuti1
	//return 8 * Eigen::Vector3d(8 * x * (x * x + y * y + z * z), y * (-3 + 8 * x * x + 8 * y * y + 8 * z * z), z * (3 + 8 * x * x + 8 * y * y + 8 * z * z));	//xingqiuti2
	//double ele = 15 * x * x + 15 * y * y + 15 * z * z;
	//return 60 * Eigen::Vector3d(x * (-7 + ele), y * (-7 + ele), z * (1 + ele));			//jiaoChaHuanMian1
	//return 60 * Eigen::Vector3d(x * (1 + ele), y * (-7 + ele), z * (-7 + ele));		//jiaoChaHuanMian2
	//double t = -1 + x * x + y * y + z * z;
	//return Eigen::Vector3d(-2 + 4 * x * t, 4 * y * t, -2 * z + 4 * z * t); // peanut


	return Eigen::Vector3d(8 * (48 * std::pow(x, 5) + 96 * x * x * z + 8 * z * (-1 - 4 * y + 4 * y * y + 4 * z * z) +
		8 * std::pow(x, 3) * (-3 - 4 * y + 12 * y * y + 12 * z * z) + x * (3 - 32 * std::pow(y, 3) + 48 * std::pow(y, 4) -
			56 * z * z + 48 * std::pow(z, 4) + y * (8 - 32 * z * z) + 8 * y * y * (-5 + 12 * z * z))),

		64 * x * (-4 + 8 * y) * z + 8 * (-1 + 2 * y) * (-1 + 4 * x * x - 4 * y + 4 * y * y + 4 * z * z) * (-1 + 4 * x * x + 4 * y + 4 * y * y + 4 * z * z) +
		(4 + 8 * y) * (-32 * z * z + std::pow(-1 + 4 * x * x - 4 * y + 4 * y * y + 4 * z * z, 2)),

		8 * (64 * x * z * z - 32 * std::pow(z, 3) + 8 * x * (-1 + 4 * x * x - 4 * y + 4 * y * y + 4 * z * z) +
			z * std::pow((-1 + 4 * x * x - 4 * y + 4 * y * y + 4 * z * z), 2) +
			2 * z * (-5 + 4 * x * x - 4 * y + 4 * y * y + 4 * z * z) * (-1 + 4 * x * x + 4 * y + 4 * y * y + 4 * z * z))); // klein bottle - 适用于hole/bettle/carp_muscle model

	//return Eigen::Vector3d(-0.5 * x * z * z + 0.75 * x * std::pow(x * x + z * z, 2), 0, -0.5 * x * x * z + 0.75 * z * std::pow(x * x + z * z, 2)); // QUATREFOIL CURVE - 适用于hole model
	//return Eigen::Vector3d(8 * x * (-4 * y * y + 48 * std::pow(x * x + y * y, 2), 8 * y * (-4 * x * x + 48 * std::pow(x * x + y * y, 2), 0))); // QUATREFOIL CURVE - 适用于bettle model
	//return Eigen::Vector3d(-2 * x * y * y + 6 * x * std::pow(x * x + y * y, 2), -2 * x * x * y + 6 * y * std::pow(x * x + y * y, 2), 0);  // QUATREFOIL CURVE - 适用于carp_muscle model

	//return Eigen::Vector3d(16 * x * (y - 2 * z) - 4 * z, -2 + 8 * x * x + 24 * y * y - 32 * y * z + 16 * z * z, -4 * (x + 4 * x * x + 4 * y * (y - 2 * z))); // mobius - 适用于hole/bettle model
	//return Eigen::Vector3d(2 * x * (y - 2 * z) - 2 * z, -1 + x * x + 3 * y * y - 4 * y * z + z * z, -2 * (x + x * x + y * (y - z))); // mobius - 适用于carp_muscle model
}

int main(int argc, char** argv)
{
	// parse input arguments
	CLI::App app("Mesh Scalar Cutting");
	//app.add_flag("-h,--help", "Print configuration and exit.");

	std::string modelArg = R"(test2d2_star.obj)";
	app.add_option("-f,--file", modelArg,
		"Input model's name with extension.")/*->required()*/;

	int numSamplePoints = 4;
	app.add_option("-n", numSamplePoints,
		"Specify the number of sample points at each edge of input model.");

	int iter = 20;
	app.add_option("-k,--iter", iter);

	try {
		argv = app.ensure_utf8(argv);
		app.parse(argc, argv);
	}
	catch (const CLI::ParseError& e) {
		// 输出帮助信息
		return app.exit(e);
	}

	const std::string mesh_file(MODEL_DIR + R"(\)" + modelArg);
	const std::string modelName = str_util::getFileName(mesh_file);

	ScalarFunc scalarFunc = { val, grad };

	// 传入不可导点
	std::vector<Vector3> singulars;
	singulars.emplace_back(Vector3(0, 0, 0));

	std::string meshNormVisDir;

#define MESH_NORM 0
#if MESH_NORM
	// norm版构造函数
	meshNormVisDir = str_util::concatFilePath(VIS_DIR, modelName + "_norm");
	const std::string norm_mesh_file(str_util::concatFilePath(meshNormVisDir, modelName + "_norm.obj"));
	double scaleFactor = 1.0; // 用于norm后的缩放，若不norm则无需传递
	core::MSCuttingModel mscModel(mesh_file, numSamplePoints, scalarFunc, singulars, norm_mesh_file, scaleFactor);
#else
	// 非norm版构造函数
	meshNormVisDir = str_util::concatFilePath(VIS_DIR, modelName);
	core::MSCuttingModel mscModel(mesh_file, numSamplePoints, scalarFunc, singulars);
#endif

	const std::string mae_vis_file = str_util::concatFilePath(
		meshNormVisDir, std::to_string(numSamplePoints), funcName, "circle", "mae.txt"
	);
	str_util::checkDir(mae_vis_file);
	std::ofstream mae_out(mae_vis_file);
	mae_out << funcName << "\n";

	for (int i = 20; i <= iter; ++i)
	{
		const std::string pd_vis_file = str_util::concatFilePath(
			meshNormVisDir, std::to_string(numSamplePoints), funcName, "circle", "isoline_" + std::to_string(i) + ".obj"
		);

		double mae = mscModel.launch(i, mae_out, pd_vis_file);
	}
	mae_out.close();

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