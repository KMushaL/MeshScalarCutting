#pragma once
#include "Config.hpp"
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>

NAMESPACE_BEGIN(mscut)

namespace donut {
	template <class T, class F, T... inds>
	constexpr _FORCE_INLINE_ inline void loop(std::integer_sequence<T, inds...>, F&& f) {
		(f(std::integral_constant<T, inds>{}), ...);
	}
	template <class T, T count, class F>
	constexpr _FORCE_INLINE_ inline void Loop(F&& f) {
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

namespace LOG {

	// word type
	namespace wType {
		const std::string NONE("0"), BOLD("1"), DIM("2"), UNDERLINE("4"), BLINK("5"),
			INVERSE("7"), HIDDEN("8");
	}
	// word color
	namespace wColor {
		const std::string BLACK("30"), RED("31"), GREEN("32"), YELLOW("33"), BLUE("34"),
			MAGENTA("35"), CYAN("36");
	}

	template <typename T, typename... Ts> inline void qpCtrl(T v, Ts... vl) {
		std::string ctrl("\033[");
		ctrl += std::string(v);
		if constexpr (sizeof...(vl) > 0) {
			std::array cl = { std::string(vl)... };
			for (auto& c : cl)
				ctrl += ";" + c;
		}
		ctrl += "m";
		std::cout << ctrl;
	}
	inline void qpCtrl() { std::cout << "\033[0m"; }

	// Print the values without line break.
	template <typename T, typename... Ts> inline void qprint_nlb(T v, Ts... vl) {
		std::cout << v;
		if constexpr (sizeof...(vl) > 0) {
			qprint_nlb(vl...);
			return;
		}
		// std::cout << std::endl;
	}

	// Print the values with line break.
	template <typename T, typename... Ts> inline void qprint(T v, Ts... vl) {
		std::cout << v;
		if constexpr (sizeof...(vl) > 0) {
			qprint(vl...);
			return;
		}
		std::cout << std::endl;
	}

	inline void qprint() { printf("\n"); }

	template <typename T, typename... Ts>
	inline void	qpDebug(T v, Ts... vl) {
		LOG::qpCtrl(LOG::wColor::BLUE);
		LOG::qpCtrl(LOG::wType::BOLD);
		std::cout << "[DEBUG] ";
		LOG::qpCtrl();
		LOG::qprint(v, vl...);
	}

	template <typename T, typename... Ts>
	inline void	qpInfo(T v, Ts... vl) {
		LOG::qpCtrl(LOG::wColor::CYAN);
		LOG::qpCtrl(LOG::wType::BOLD);
		std::cout << "[INFO] ";
		LOG::qpCtrl();
		LOG::qprint(v, vl...);
	}

	template <typename T, typename... Ts>
	inline void	qpWarn(T v, Ts... vl) {
		LOG::qpCtrl(LOG::wColor::YELLOW);
		LOG::qpCtrl(LOG::wType::BOLD);
		std::cout << "[WARN] ";
		LOG::qpCtrl();
		LOG::qprint(v, vl...);
	}

	template <typename T, typename... Ts>
	inline void	qpError(T v, Ts... vl) {
		LOG::qpCtrl(LOG::wColor::RED);
		LOG::qpCtrl(LOG::wType::BOLD);
		std::cout << "[ERROR] ";
		LOG::qpCtrl();
		LOG::qprint(v, vl...);
	}

	template <typename T, typename... Ts>
	inline void	qpTest(T v, Ts... vl) {
		LOG::qpCtrl(LOG::wColor::GREEN);
		LOG::qpCtrl(LOG::wType::BOLD);
		std::cout << "[TEST] ";
		LOG::qpCtrl();
		LOG::qprint(v, vl...);
	}

	template <typename T, typename... Ts>
	inline void	qpNormal(T v, Ts... vl) {
		std::cout << "-- ";
		LOG::qpCtrl();
		LOG::qprint(v, vl...);
	}

	inline void qpSplit()
	{
		LOG::qprint("==========");
	}

} // namespace LOG
NAMESPACE_END(mscut)