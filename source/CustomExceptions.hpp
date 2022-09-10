#ifndef MC_CustomExceptions_H_
#define MC_CustomExceptions_H_

#include <stdexcept>

namespace MC::EX
{
    // Exception for reporting invalid user input
    class TInvalidInput : public std::invalid_argument
    {
        public:
            TInvalidInput(const std::string& msg);
    };

    // Exception for reporting file access errors
    class TFileAccess : public std::runtime_error
    {
        public:
            TFileAccess(const std::string& msg);
    };

    // Exception for reporting inconsistencies in program flow
    class TInvalidStatus : public std::logic_error
    {
        public:
            TInvalidStatus(const std::string& msg, const std::string& func_name);
    };

    // Exception for reporting memory allocation errors
    class TOutOfMemory : public std::runtime_error
    {
        public:
            TOutOfMemory(const std::string& msg, const std::string& func_name, const std::string& ex_msg);
    };
}

#endif