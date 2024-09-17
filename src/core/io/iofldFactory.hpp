#ifndef IOFLDFACTORY_HPP
#define IOFLDFACTORY_HPP

#include "iofld.hpp"
#include <memory>
#include <string>

class iofldFactory {
public:
    static std::unique_ptr<iofld> create(const std::string& engineType = "");
};

#endif // IOFLDFACTORY_HPP
