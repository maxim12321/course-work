#pragma once

#include <vector>
#include <string>
#include <cassert>
#include <fstream>

std::vector<std::string> SplitString(std::string str, std::string delim);

// Read line from stream and check if stream is not eof
std::string ReadLine(std::istream& in);
