#pragma once

#if defined(__clang__) || defined(__GNUC__)
#  define FORCE_INLINE __attribute__((always_inline)) inline
#elif defined(_MSC_VER)
#  define FORCE_INLINE __forceinline
#endif

/* namespace macro */
#if !defined(NAMESPACE_BEGIN)
#  define NAMESPACE_BEGIN(name) namespace name {
#endif

#if !defined(NAMESPACE_END)
#  define NAMESPACE_END(name) }
#endif

#ifdef __GNUC__
#define OFFSET_FUNCTION __PRETTY_FUNCTION__
#elif defined(__clang__) || (_MSC_VER >= 1310)
#define OFFSET_FUNCTION __FUNCTION__
#else
#define OFFSET_FUNCTION "unknown"
#endif

/* CUDA macro */
#ifdef __CUDACC__
#  define _CUDA_HOST_CALL_ __host__
#else
#  define _CUDA_HOST_CALL_
#endif

#ifdef __CUDACC__
#  define _CUDA_DEVICE_CALL_ __device__
#else
#  define _CUDA_DEVICE_CALL_
#endif

#ifdef __CUDACC__
#  define _CUDA_GENERAL_CALL_ __host__ __device__
#else
#  define _CUDA_GENERAL_CALL_
#endif

#ifdef __CUDACC__
#  define _CUDA_KERNEL_CALL_ __global__
#else
#  define _CUDA_KERNEL_CALL_
#endif

#ifndef MODEL_DIR
#  define MODEL_DIR (std::string)".\\model"
#endif
#ifndef SHARED_PATH
#  define SHARED_PATH (std::string)".\\data"
#endif
#ifndef VIS_DIR
#  define VIS_DIR (std::string)".\\vis"
#endif
#ifndef OUT_DIR
#  define OUT_DIR (std::string)".\\output"
#endif