#pragma once
#include "Config.hpp"
#include "core/MSCuttingModel.hpp"

NAMESPACE_BEGIN(mscut)
namespace unit_test
{
	using namespace detail;
	using namespace core;
	using ScalarFunc = core::MSCuttingModel::ScalarFunc;

	inline bool testSamplingPoint(const std::string& mesh_file, const ScalarFunc& func, int m_numSamples, int numSamples)
	{
		if (numSamples < 0)
		{
			LOG::qpTest("The number of sample points must be larger than 0!");
			return false;
		}
		std::vector<Vector3> singulars;
		MSCuttingModel mscModel(mesh_file, func, singulars, m_numSamples, numSamples);
		mscModel.meshNorm();

		const std::string vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"sample.obj");
		str_util::checkDir(vis_file);
		std::ofstream vis_out(vis_file);
		if (!vis_out) { LOG::qpError("I/O: File ", vis_file.c_str(), " could not be opened!"); return false; }

		const std::string valid_vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"valid_sample.obj");
		str_util::checkDir(valid_vis_file);
		std::ofstream valid_vis_out(valid_vis_file);
		if (!vis_out) { LOG::qpError("I/O: File ", valid_vis_file.c_str(), " could not be opened!"); return false; }

		LOG::qpTest("Output sample points to ", std::quoted(vis_file), " ...");
		LOG::qpTest("Output valid sample points to ", std::quoted(valid_vis_file), " ...");

		bool res = mscModel.testSamplingPoints(vis_out, valid_vis_out);

		vis_out.close();
		valid_vis_out.close();
		return res;
	}

} // namespace unit_test
NAMESPACE_END(mscut)