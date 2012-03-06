
#ifndef __STRINGUTILS_H_
#define __STRINGUTILS_H_

#include <string>
#include <sstream>
#include <vector>

//copied from Paul J. Weiss and modified code!
int SplitString(const std::string& input, 
		const std::string& delimiter,
		std::vector<std::string>& results, 
		bool includeEmpties = true);

int str2int(std::string &s);
double str2double(std::string &s);
long str2long(std::string &s);
void chomp(std::string& s);
void strim(std::string& s);

#endif
