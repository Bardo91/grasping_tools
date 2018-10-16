/*
 * utilsMisc.h
 *
 *  Created on: 10 January 2016
 *      Author: rfitch
 *
 *      static utilities for miscellaneous things like throwing errors and debugging
 */

#ifndef UTILSMISC_H
#define UTILSMISC_H

#include <iostream>
#include <chrono>

namespace gpis {
	template<class T> inline void debugPrint(const T& data) {
		std::cout << data << std::endl;
	}

	template<class T> inline void error(const T& data) {
		std::cerr << data << std::endl;
	}

	inline void error(const std::ostringstream &errorString) {
		std::cerr << errorString.str() << std::endl;
	}

	inline double timeDifferenceMillis(const clock_t time0, const clock_t time1) {
		return ((time1 - time0) / (double)CLOCKS_PER_SEC) * 1000.0;
	}

	inline double getUnixTimeMillis() {
		using namespace std::chrono;

		return time_point_cast<microseconds>(steady_clock::now()).time_since_epoch().count()/1000;

		//    // On the Mac
		//    clock_t currentTime = clock();
		//    return double(currentTime) / CLOCKS_PER_SEC * 1000.0;
	}
}	// namespace gpis
	
#endif // UTILSMISC_H
