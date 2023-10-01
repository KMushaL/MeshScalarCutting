#pragma once
#include "Config.hpp"
#include "core/MSCuttingModel.hpp"

NAMESPACE_BEGIN(mscut)
namespace unit_test
{
	using namespace detail;
	using namespace core;

	inline bool testComputeAD(const std::string& mesh_file, int numSamples, int facetIdx = -1)
	{
		if (numSamples <= 0)
		{
			LOG::qpTest("The number of sample points must larger than 0!");
			return false;
		}
		MSCuttingModel mscModel(mesh_file, numSamples);

		const std::string vis_file = str_util::concatFilePath(VIS_DIR, mscModel.modelName, (std::string)"apollonius_diagram.obj");
		str_util::checkDir(vis_file);
		std::ofstream vis_out(vis_file);
		if (!vis_out) { LOG::qpError("I/O: File ", vis_file.c_str(), " could not be opened!"); return false; }

		LOG::qpTest("Output Apollonius Diagram to ", std::quoted(vis_file), " ...");

		bool res = false;
		if (facetIdx <= -1)
		{
			int numFaces = mscModel.numFaces();
			/*for (int i = 0; i < numFaces; ++i)
			{
				mscModel.testComputeADForFacet(i, vis_out);
			}*/
			res = mscModel.testComputeADForFacet(0, vis_out);
		}
		else
			res = mscModel.testComputeADForFacet(facetIdx, vis_out);

		vis_out.close();
		return res;
	}

} // namespace unit_test
NAMESPACE_END(mscut)