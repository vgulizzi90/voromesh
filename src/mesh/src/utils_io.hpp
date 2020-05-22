// utils_io.hpp

#ifndef VOROMESH_IO_HPP_
#define VOROMESH_IO_HPP_

#include <iostream>
#include <string>

namespace voromesh
{
namespace io
{

void error(const std::string &, const std::string &);
void warning(const std::string &, const std::string &);

}
}

#endif