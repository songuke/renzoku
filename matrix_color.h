#ifndef _MATRIX_COLOR_H_
#define _MATRIX_COLOR_H_

#include <cstdio>
#include <memory>
#include <typeinfo>
#include <exception>

#include "error.h"

using namespace std;

namespace Renzoku {

/** 
 * Matrix with column-major internal storage. The third dimension stores color.
 */
template <class T>
class MatrixColor {
public:
	MatrixColor();
	MatrixColor(const MatrixColor<T> &m);

    /**
     * Defautl constructor of type T is called.
     * In case of int, float, double: what is the default value of int()?
     */
	MatrixColor(int rows, int cols);	
    MatrixColor(int rows, int cols, T val);
	~MatrixColor();

public:
    int get_rows() const;
    int get_cols() const;
	int get_colors() const;
	
public:
    // The array operator [] only accepts a single argument so it cannot be used to accept row and column index.
    // To be consistent, we use operator() to access elements of the matrix.
	T operator()(int i) const;
	T& operator()(int i);
	T operator()(int i, int j, int k) const;
	T& operator()(int i, int j, int k);

public:
	void load(const char *file);

    /**     
     * RGB matrix format:
     * 
     * <width>(space)<height>\n
     * binary data with RGB of pixel 1, RGB of pixel 2, in column major order.
     */
	void save(const char *file);
	
public:
    /**
     * Modify the storage of a matrix by setting a new number of rows.
     */
    void resize(int new_rows);

    /**
     * Resize the matrix to new rows/columns. 
     * 
     * The internal memory is reallocated if necessary. Old matrix data becomes invalid.
     */
    void resize(int new_rows, int new_cols);

protected:
    /**
     * Ensure only compatible types are used with the matrix.
     */
    static bool check_typeid() throw();

protected:
	T *data;
    int max_rows, max_cols;
	int rows, cols, colors;
};

template <class T>
bool MatrixColor<T>::check_typeid() throw() {
    if (typeid(T) == typeid(char)  || typeid(T) == typeid(int)    || typeid(T) == typeid(long) || typeid(T) == typeid(long long) ||
        typeid(T) == typeid(float) || typeid(T) == typeid(double) ||         
        typeid(T) == typeid(unsigned int)  || typeid(T) == typeid(unsigned char) || 
        typeid(T) == typeid(unsigned long) || typeid(T) == typeid(unsigned long long)
        ) return true;

    throw TypeException("Matrix type not supported.");
    
    return false;
}

template <class T>
MatrixColor<T>::MatrixColor() : data(NULL), rows(0), cols(0), colors(3), max_rows(0), max_cols(0) {	

}

template <class T>
MatrixColor<T>::MatrixColor(int rows, int cols) : rows(rows), cols(cols), colors(3), max_rows(rows), max_cols(cols) {
	data = new T[max_rows * max_cols * colors];
	if (!data) error_alloc(__FILE__, __LINE__);
	for (int i = 0; i < rows * cols * colors; ++i)
		data[i] = T();

}

template <class T>
MatrixColor<T>::MatrixColor(int rows, int cols, T val) : rows(rows), cols(cols), colors(3), max_rows(rows), max_cols(cols) {
	data = new T[max_rows * max_cols * colors];
	if (!data) error_alloc(__FILE__, __LINE__);
	for (int i = 0; i < rows * cols * colors; ++i) 
		data[i] = val;
}

template <class T>
MatrixColor<T>::MatrixColor(const MatrixColor<T> &m) : rows(m.rows), cols(m.cols), colors(m.colors), max_rows(m.max_rows), max_cols(m.max_cols) {
	data = new T[max_rows * max_cols * colors];
	if (!data) error_alloc(__FILE__, __LINE__);
	memcpy(data, m.data, sizeof(T) * rows * cols * colors);
}

template <class T>
MatrixColor<T>::~MatrixColor() {
	if (data)
		delete [] data;	
}

template <class T>
int MatrixColor<T>::get_rows() const { return rows; }

template <class T>
int MatrixColor<T>::get_cols() const { return cols; }

template <class T>
int MatrixColor<T>::get_colors() const { return colors; }

template <class T>
inline T MatrixColor<T>::operator()(int i) const {
	return data[i];
}

template <class T>
inline T& MatrixColor<T>::operator()(int i) {
	return data[i];
}

template <class T>
inline T MatrixColor<T>::operator()(int i, int j, int k) const {
	return data[k * cols * rows + j * rows + i];
}

template <class T>
inline T& MatrixColor<T>::operator()(int i, int j, int k) {
	return data[k * cols * rows + j * rows + i];
}

template <class T>
void MatrixColor<T>::save(const char *file) {
    check_typeid();

	FILE *f = fopen(file, "wb");
    if (f == NULL) 
        throw FileException("Cannot open file.");

	fprintf(f, "%d %d\n", rows, cols);    
	fwrite(data, sizeof(T), rows * cols * colors, f);
    fclose(f);
}

template <class T>
void MatrixColor<T>::load(const char *file) {
    check_typeid();

	FILE *f = fopen(file, "rb");
    if (f == NULL)
        throw FileException("Cannot open file.");

	fscanf(f, "%d %d\n", &rows, &cols);
	max_rows = rows;
	max_cols = cols;
	colors = 3;
	data = new T[rows * cols * colors];
	fread(data, sizeof(T), rows * cols, f);
	fclose(f);
}

template <class T>
void MatrixColor<T>::resize(int new_rows) {
    resize(new_rows, cols);
}

template <class T>
void MatrixColor<T>::resize(int new_rows, int new_cols) {
    if (new_rows * new_cols <= max_rows * max_cols) {
        rows = new_rows;
        cols = new_cols;
    } else {
        // reallocate
        rows = new_rows;
        cols = new_cols;
        max_rows = new_rows;
        max_cols = new_cols;
        if (data) delete [] data;
        data = new T[max_rows * max_cols * colors];
    }
}

} // end namespace

#endif
