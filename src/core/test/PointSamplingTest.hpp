#pragma once
#include "Config.hpp"
#include "core/MSCuttingModel.hpp"

NAMESPACE_BEGIN(mscut)
namespace unit_test
{
	using namespace detail;
	using namespace core;

	inline bool testSamplingPoint(const std::string& mesh_file, int numSamples)
	{
		if (numSamples <= 0)
		{
			LOG::qpTest("The number of sample points must larger than 0!");
			return false;
		}
		MSCuttingModel mscModel(mesh_file, numSamples);

		const std::string vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"sample.obj");
		str_util::checkDir(vis_file);
		std::ofstream vis_out(vis_file);
		if (!vis_out) { LOG::qpError("I/O: File ", vis_file.c_str(), " could not be opened!"); return false; }

		LOG::qpTest("Output sample points to ", std::quoted(vis_file), " ...");

		bool res = mscModel.testSamplingPoints(vis_out);

		vis_out.close();
		return res;
	}

} // namespace unit_test
NAMESPACE_END(mscut)