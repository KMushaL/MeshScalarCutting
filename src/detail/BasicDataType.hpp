#pragma once

#include "Config.hpp"
#include <numeric>
#include <array>
#include <Eigen/Dense>
#include <unordered_set>

namespace std
{
	template<typename T, size_t N>
	struct hash<array<T, N> >
	{
		typedef array<T, N> argument_type;
		typedef size_t result_type;

		result_type operator()(const argument_type& a) const
		{
			hash<T> hasher;
			result_type h = 0;
			for (result_type i = 0; i < N; ++i)
			{
				h = h * 31 + hasher(a[i]);
			}
			return h;
		}
	};
}

NAMESPACE_BEGIN(mscut)
namespace detail
{
	typedef unsigned int uint;

#ifdef USING_FLOAT
	typedef float Scalar;
	constexpr Scalar EQUAL_EPSILON = 1e-6;
	constexpr Scalar ZERO_EPSILON = 1e-6;
	const Scalar DINF = std::numeric_limits<float>::max();
#else
	typedef double Scalar;
	constexpr Scalar EQUAL_EPSILON = 1e-8;
	constexpr Scalar ZERO_EPSILON = 1e-8;
	const Scalar DINF = std::numeric_limits<double>::max();
#endif

	template<typename _Scalar, int Dim>
	using Vector = typename Eigen::Matrix<_Scalar, Dim, 1>;

	using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

	using VectorXi = Eigen::VectorXi;

	using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

	using Vector3i = Eigen::Vector3i;

	using Array3 = Eigen::Array<Scalar, 3, 1>;

	template<typename _Scalar, int Rows, int Cols>
	using Matrix = typename Eigen::Matrix<_Scalar, Rows, Cols>;

	using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

	using MatrixXi = Eigen::MatrixXi;

	template<typename _Scalar, int Num>
	using ArrayX = typename Eigen::Array<_Scalar, Num, 1>;

	template<typename _Scalar, int Rows, int Cols>
	using ArrayXX = typename Eigen::Array<_Scalar, Rows, Cols>;
}
NAMESPACE_END(mscut)