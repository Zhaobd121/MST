#pragma once
#include <cstdio>
#include <cstdlib>

class Utility {
public:
	static FILE *open_file(const char *file_name, const char *mode) {
		FILE *f = fopen(file_name, mode);
		if(f == nullptr) {
			printf("Can not open file: %s\n", file_name);
			exit(1);
		}

		return f;
	}

	static std::string integer_to_string(long long number) {
		std::vector<uint32_t> sequence;
		if(number == 0) sequence.push_back(0);
		while(number > 0) {
			sequence.push_back(number%1000);
			number /= 1000;
		}

		char buf[5];
		std::string res;
		for(uint32_t i = sequence.size();i > 0;i --) {
			if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
			else sprintf(buf, ",%03u", sequence[i-1]);
			res += std::string(buf);
		}
		return res;
	}
};