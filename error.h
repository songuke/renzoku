#ifndef _ERROR_H_
#define _ERROR_H_

#include <exception>
#include <string>

using namespace std;

namespace Renzoku {

void throw_error(const char *error, const char *file, int line);
void throw_warning(const char *warn, const char *file, int line);

void error_alloc(const char *file, int line);
void error_file(const char *file, int line, const char *filename);

class BaseException : public std::exception {
public:
    BaseException() : s("") {}
    BaseException(string s) : s(s) {}
    virtual ~BaseException() throw() {}

    virtual const char *what() const throw() {
        return s.c_str();
    }

protected:
    string s;
};

/**
 * TypeException is thrown when a type is passed into a template class 
 * but it is not fully supported. 
 * 
 * For example, when an arbitrary object 
 * is passed into Matrix<T> class, the load/save function does not work
 * as it does not know how to load/save the type T. 
 * Therefore, this exception must be thrown when the load/save is called 
 * on the unknown type. 
 */ 
class TypeException : public BaseException {
public:
    TypeException(string s) : BaseException(s) {}
};

class FileException : public BaseException {
public:
    FileException(string s) : BaseException(s) {}
};

/**
 * NotEnoughCapacityException is thrown when a matrix fails to expand its rows 
 * due to the maximum capacity of the matrix is reached. 
 * 
 * New allocation of matrix internal data is required.
 */
class NotEnoughCapacityException : public BaseException {
public:
    NotEnoughCapacityException(string s) : BaseException(s) {}
};

/**
 * IncompatibleSizeException is thrown when an operation that requires two vectors
 * with the same size is called but the input is not.
 */
class IncompatibleSizeException : public BaseException {
public:
    IncompatibleSizeException() {}
};

} // end namespace 

#endif
