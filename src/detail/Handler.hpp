#pragma once
#include "Config.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>

NAMESPACE_BEGIN(mscut)
namespace detail
{
#define OFFSET_FMT_CSTR(description, ...) description
#define OFFSET_FMT_STR(description, ...) std::string(description)
#define OFFSET_FMT_PRINT(description, ...) /*std::printf("%s\n", description)*/
#define OFFSET_FMT_ARG(arg)

#if defined(OFFSET_DISABLE_ENSURES)

#define OFFSET_ENSURE(expr, ...) ((void)0)

#elif defined(OFFSET_ENABLE_ENSURE_HANDLER)

	void ensureFailed(char const* function, char const* file, int line,
		char const* description);

#define OFFSET_ENSURE(expr, description, ...)                                  \
  ((expr)                                                                      \
       ? ((void)0)                                                             \
       : ::OFFSET::ensureFailed(OFFSET_FUNCTION, __FILE__, __LINE__,           \
                                OFFSET_FMT_CSTR(description, ##__VA_ARGS__)))

#else

#define OFFSET_DEDAULT_ENSURE_FAILURE_IMPL(function, file, line, description, \
                                           ...)                                \
  do {                                                                         \
    std::printf("Assert failed in function '%s', "                             \
                "file '%s', line %d.\n",                                       \
                function, file, line);                                         \
    OFFSET_FMT_PRINT(description, ##__VA_ARGS__);                              \
    std::abort();                                                              \
  } while (0)

#ifdef __CUDACC__
#define OFFSET_ENSURE(expr, description, ...)                                  \
  do {                                                                         \
    if (!(expr)) {                                                             \
      std::printf("Assert failed in function '%s', file '%s', line %d.\n",     \
                  OFFSET_FUNCTION, __FILE__, __LINE__);                        \
      std::printf("%s", description);                                          \
      /* there is no std::abort in cuda kernels, hence we just print the error \
       * message here*/                                                        \
    }                                                                          \
  } while (0)
#else
#define OFFSET_ENSURE(expr, ...)                                               \
  do {                                                                         \
    if (!(expr)) {                                                             \
      OFFSET_DEDAULT_ENSURE_FAILURE_IMPL(OFFSET_FUNCTION, __FILE__, __LINE__,  \
                                         ##__VA_ARGS__);                       \
    }                                                                          \
  } while (0)
#endif

#endif

} // namespace detail
NAMESPACE_END(mscut)