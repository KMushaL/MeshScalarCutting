#pragma once
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

namespace donut {

	template <class T, class F, T... inds>
	constexpr inline void loop(std::integer_sequence<T, inds...>, F&& f) {
		(f(std::integral_constant<T, inds>{}), ...);
	}
	template <class T, T count, class F>
	constexpr inline void Loop(F&& f) {
		loop(std::make_integer_sequence<T, count>{}, std::forward<F>(f));
	}

	template<typename T>
	inline T max(const T& value)
	{
		return value;
	}
	template<typename T, typename...Args>
	inline T max(const T& value, Args&&...args)
	{
		return std::max(value, max(std::forward<Args>(args)...));
	}

	template<typename T>
	inline T min(const T& value)
	{
		return value;
	}
	template<typename T, typename...Args>
	inline T min(const T& value, Args&&...args)
	{
		return std::min(value, min(std::forward<Args>(args)...));
	}

	template<typename T, typename...Args>
	inline bool isSandwitchVal(const T& value, Args&&...args)
	{
		T minArg = min(std::forward<Args>(args)...);
		T maxArg = max(std::forward<Args>(args)...);
		//return (minArg < value && value < maxArg);
		return (minArg <= value && value <= maxArg);
	}

	template<typename T>
	inline void removeDupicates(std::vector<T>& vec) {
		auto endIter{ vec.end() };
		for (auto iter = vec.begin(); iter != endIter; ++iter)
			endIter = std::remove(iter + 1, endIter, *iter);
		vec.erase(endIter, vec.end());
	}

} // namespace donut