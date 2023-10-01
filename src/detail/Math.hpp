#pragma once
#include "Config.hpp"
#include "Handler.hpp"
#include "BasicDataType.hpp"
#include <type_traits>

NAMESPACE_BEGIN(mscut)
namespace detail 
{
	template<class _Scalar>
	class MaxMetric {
	public:
		static _Scalar impl(_Scalar s0, _Scalar s1) {
			using std::abs;
			return abs(s0 - s1);
		}
	};

	template<class _Scalar, int M, int N>
	class MaxMetric<Matrix<_Scalar, M, N>> {
	public:
		static _Scalar impl(Matrix<_Scalar, M, N> const& p0,
			Matrix<_Scalar, M, N> const& p1) {
			return (p0 - p1).template lpNorm<Eigen::Infinity>();
		}
	};

	template<class _Scalar>
	class SetToZero {
	public:
		static void impl(_Scalar& s) { s = _Scalar(0); }
	};

	template<class _Scalar, int M, int N>
	class SetToZero<Matrix<_Scalar, M, N>> {
	public:
		static void impl(Matrix<_Scalar, M, N>& v) { v.setZero(); }
	};

	template<class T1, class _Scalar>
	class SetElementAt;

	template<class _Scalar>
	class SetElementAt<_Scalar, _Scalar> {
	public:
		static void impl(_Scalar& s, _Scalar value, int at) {
			OFFSET_ENSURE(at == 0, "is {}", OFFSET_FMT_ARG(at));
			s = value;
		}
	};

	template<class _Scalar, int N>
	class SetElementAt<Vector<_Scalar, N>, _Scalar> {
	public:
		static void impl(Vector<_Scalar, N>& v, _Scalar value, int at) {
			//OFFSET_ENSURE(at >= 0 && at < N, "is {}", OFFSET_FMT_ARG(at));
			v[at] = value;
		}
	};

	template<class _Scalar>
	class SquaredNorm {
	public:
		static _Scalar impl(_Scalar const& s) { return s * s; }
	};

	template<class _Scalar, int N>
	class SquaredNorm<Matrix<_Scalar, N, 1>> {
	public:
		static _Scalar impl(Matrix<_Scalar, N, 1> const& s) { return s.squaredNorm(); }
	};

	template<class _Scalar>
	class Transpose {
	public:
		static _Scalar impl(_Scalar const& s) { return s; }
	};

	template<class _Scalar, int M, int N>
	class Transpose<Matrix<_Scalar, M, N>> {
	public:
		static Matrix<_Scalar, M, N> impl(Matrix<_Scalar, M, N> const& s) {
			return s.transpose();
		}
	};

	/// Returns maximum metric between two points ``p0`` and ``p1``, with ``p0, p1``
	/// being matrices or a _Scalars.
	///
	template<class T>
	auto maxMetric(T const& p0, T const& p1)
		-> decltype(detail::MaxMetric<T>::impl(p0, p1)) {
		return detail::MaxMetric<T>::impl(p0, p1);
	}

	/// Sets point ``p`` to zero, with ``p`` being a matrix or a _Scalar.
	///
	template<class T>
	void setToZero(T& p) {
		return detail::SetToZero<T>::impl(p);
	}

	/// Sets ``i``th component of ``p`` to ``value``, with ``p`` being a
	/// matrix or a _Scalar. If ``p`` is a _Scalar, ``i`` must be ``0``.
	///
	template<class T, class _Scalar>
	void setElementAt(T& p, _Scalar value, int i) {
		return detail::SetElementAt<T, _Scalar>::impl(p, value, i);
	}

	/// Returns the squared 2-norm of ``p``, with ``p`` being a vector or a _Scalar.
	///
	template<class T>
	auto squaredNorm(T const& p) -> decltype(detail::SquaredNorm<T>::impl(p)) {
		return detail::SquaredNorm<T>::impl(p);
	}

	/// Returns ``p.transpose()`` if ``p`` is a matrix, and simply ``p`` if m is a
	/// _Scalar.
	///
	template<class T>
	auto transpose(T const& p) -> decltype(detail::Transpose<T>::impl(T())) {
		return detail::Transpose<T>::impl(p);
	}

	template<class _Scalar>
	struct IsFloatingPoint {
		static bool const value = std::is_floating_point<_Scalar>::value;
	};

	template<class _Scalar, int M, int N>
	struct IsFloatingPoint<Matrix<_Scalar, M, N>> {
		static bool const value = std::is_floating_point<_Scalar>::value;
	};

	template<class _Scalar_>
	struct Get_Scalar {
		using _Scalar = _Scalar_;
	};

	template<class _Scalar_, int M, int N>
	struct Get_Scalar<Matrix<_Scalar_, M, N>> {
		using _Scalar = _Scalar_;
	};

	/// If the Vector type is of fixed size, then IsFixedSizeVector::value will be
	/// true.
	template<typename Vector, int NumDimensions,
		typename = typename std::enable_if<
		Vector::RowsAtCompileTime == NumDimensions &&
		Vector::ColsAtCompileTime == 1>::type>
	struct IsFixedSizeVector : std::true_type {
	};

	namespace {
		static constexpr Scalar GEO_EPSILON = 1E-6;
	}

	template<typename T>
	int dcmp(const T& x, const T& epsilon) {
		if (std::abs(x) < epsilon) return 0;
		else return x < 0 ? -1 : 1;
	}

} // namespace detail
NAMESPACE_END(mscut)