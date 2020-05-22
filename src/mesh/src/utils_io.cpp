// utils_io.cpp

#include "utils_io.hpp"

namespace voromesh
{
namespace io
{

void error(const std::string & which_function, const std::string & description)
{
    std::cout << "ERROR in " << which_function << std::endl;
    std::cout << "| " << description << std::endl;
    std::cout << std::endl;
    exit(-1);
}

void warning(const std::string & which_function, const std::string & description)
{
    std::cout << "WARNING in " << which_function << std::endl;
    std::cout << "| " << description << std::endl;
    std::cout << std::endl;
}

}
}