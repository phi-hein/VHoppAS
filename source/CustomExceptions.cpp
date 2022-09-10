#include "CustomExceptions.hpp"

#include <string>

MC::EX::TInvalidInput::TInvalidInput(const std::string& msg)
    : std::invalid_argument(msg)
{

}

MC::EX::TFileAccess::TFileAccess(const std::string& msg)
    : std::runtime_error(msg)
{

}

MC::EX::TInvalidStatus::TInvalidStatus(const std::string& msg, const std::string& func_name)
    : std::logic_error(msg + "\n[in function: " + func_name + "]")
{

}

MC::EX::TOutOfMemory::TOutOfMemory(const std::string& msg, const std::string& func_name, const std::string& ex_msg)
    : std::runtime_error(msg + "\n[in function: " + func_name + "]\n" + ex_msg)
{

}