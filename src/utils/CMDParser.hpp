#pragma once
#include <cstring>
#include <string>

namespace cmd_parser
{
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
	// 执行不区分大小写的字符串比较
	int strcasecmp(const char* c1, const char* c2) { return _stricmp(c1, c2); }
#endif // WIN

#ifndef EMPTY_ARG
#  define EMPTY_ARG   -1
#endif 
#ifndef INVALID_ARG
#  define INVALID_ARG -2
#endif

	template <typename T>
	inline void cmdLineCleanUp(T* t);

	template <typename T>
	inline T cmdLineInitialize();

	template <typename T>
	inline T cmdLineCopy(T t);

	template <typename T>
	inline T cmdLineStringToT(const char* str);

	class ParserInterface
	{
	public:
		bool set;
		char* name;

		ParserInterface(const char* name);

		virtual ~ParserInterface();

		virtual int read(int argc, char** argv);
	};

	template <typename T>
	class CmdLineParameter : public ParserInterface
	{
	public:
		T value;

		CmdLineParameter(const char* name);

		CmdLineParameter(const char* name, T v);

		~CmdLineParameter() override;

		inline int read(int argc, char** argv) override;
	};

	inline void cmdLineParse(int argc, char** argv, ParserInterface** params);

}// namespace cmd_parser

#include "CMDParser.inl"