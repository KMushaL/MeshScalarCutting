#pragma once
#include "Config.hpp"
#include "core/MSCuttingModel.hpp"

NAMESPACE_BEGIN(mscut)
namespace unit_test
{
	using namespace detail;
	using namespace core;

	inline bool testLocalGlobalTransoform(const std::string& mesh_file, int numSamples, int facetIdx = -1)
	{
		if (numSamples <= 0)
		{
			LOG::qpTest("The number of sample points must larger than 0!");
			return false;
		}
		MSCuttingModel mscModel(mesh_file, numSamples);

		if (facetIdx <= -1)
		{
			int numFaces = mscModel.numFaces();
			/*for (int i = 0; i < numFaces; ++i)
			{
				mscModel.testLocalGlobalTransform(i);
			}*/
			return (mscModel.testLocalGlobalTransform(0));
		}
		else
			return (mscModel.testLocalGlobalTransform(facetIdx));
	}

} // namespace unit_test
NAMESPACE_END(mscut)