#ifndef UTILITY_H
#define UTILITY_H
#include <sstream>
#include <vector>
#include <string>

std::vector<std::string> split(const std::string& str, char sep = ' ') {
	std::vector<std::string> res;
	std::stringstream iss(str);
	std::string word;
	while (std::getline(iss, word, sep)) {
		if (word.empty())
			continue;
		res.push_back(word);
	}
	return res;
}

inline double clamp(double val, double low, double high) {
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

#endif
