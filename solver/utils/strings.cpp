#include "strings.h"

std::vector<std::string> SplitString(std::string str, std::string delim) {
    std::vector<std::string> split;
    size_t pos = 0;
    while ((pos = str.find(delim)) != std::string::npos) {
        split.push_back(str.substr(0, pos));
        str.erase(0, pos + delim.length());
    }
    split.push_back(str);

    return std::move(split);
}

// Read line from stream and check if stream is not eof
std::string ReadLine(std::istream& in) {
    assert(in.peek() != EOF);
    std::string line;
    std::getline(in, line);
    return std::move(line);
}