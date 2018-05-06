#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <cstdio>
#include <memory>
#include <typeinfo>
#include <exception>

#include "rgb.h"
#include "error.h"

using namespace std;

namespace Renzoku {

/**
 * Matrix with row-major internal storage.
 */
template <class T>
class Matrix {
public:
	Matrix();
	Matrix(const Matrix<T> &m);

    /**
     * Defautl constructor of type T is called.
     * In case of int, float, double: what is the default value of int()?
     */
	Matrix(int rows, int cols);
	
    Matrix(int rows, int cols, T val);
	~Matrix();

public:
    int get_rows() const;
    int get_cols() const;

public:
    // The array operator [] only accepts a single argument so it cannot be used to accept row and column index.
    // To be consistent, we use operator() to access elements of the matrix.
	T operator()(int i) const;
	T& operator()(int i);
	T operator()(int i, int j) const;
	T& operator()(int i, int j);

public:
	void load(const char *file);

    /**     
     * RGB matrix format:
     * 
     * <width>(space)<height>\n
     * binary data with RGB of pixel 1, RGB of pixel 2, in row major order.
     */
	void save(const char *file) const;
	
public:
	void mul(T *in, T *out);

    /**
     * Convert 3-channel matrix into a flat matrix.
     */
    static void convert(const Matrix<Rgb> &m, Matrix<Float> &d); 
    
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
	int rows, cols;
};

// template specialization for Rgb data type. These functions are stored in matrix_c.cc to avoid duplicate definition.
template <>
void Matrix<Rgb>::load(const char *file);

template <>
void Matrix<Rgb>::save(const char *file) const;

template <class T>
Matrix<T>::Matrix() : data(NULL), rows(0), cols(0), max_rows(0), max_cols(0) {	

}

template <class T>
Matrix<T>::Matrix(int rows, int cols) : rows(rows), cols(cols), max_rows(rows), max_cols(cols) {
	data = new T[max_rows * max_cols];
	if (!data) error_alloc(__FILE__, __LINE__);
	for (int i = 0; i < rows * cols; ++i)
		data[i] = T();

}

template <class T>
Matrix<T>::Matrix(int rows, int cols, T val) : rows(rows), cols(cols), max_rows(rows), max_cols(cols) {
	data = new T[max_rows * max_cols];
	if (!data) error_alloc(__FILE__, __LINE__);
	for (int i = 0; i < rows * cols; ++i) 
		data[i] = val;
}

template <class T>
Matrix<T>::Matrix(const Matrix<T> &m) : rows(m.rows), cols(m.cols), max_rows(m.max_rows), max_cols(m.max_cols) {
	data = new T[max_rows * max_cols];
	if (!data) error_alloc(__FILE__, __LINE__);
	memcpy(data, m.data, sizeof(T) * rows * cols);
}

template <class T>
Matrix<T>::~Matrix() {
	if (data)
		delete [] data;	
}

template <class T>
int Matrix<T>::get_rows() const { return rows; }

template <class T>
int Matrix<T>::get_cols() const { return cols; }

template <class T>
T Matrix<T>::operator()(int i) const {
	return data[i];
}

template <class T>
T& Matrix<T>::operator()(int i) {
	return data[i];
}

template <class T>
T Matrix<T>::operator()(int i, int j) const {
	return data[i * cols + j];
}

template <class T>
T& Matrix<T>::operator()(int i, int j) {
	return data[i * cols + j];
}

template <class T>
void Matrix<T>::save(const char *file) const {
	FILE *f = fopen(file, "wb");
    if (f == NULL) 
        throw FileException("Cannot open file.");

	fprintf(f, "%d %d\n", rows, cols);
    fwrite(data, sizeof(T), rows * cols, f);
    fclose(f);
}

template <class T>
void Matrix<T>::load(const char *file) {
    FILE *f = fopen(file, "rb");
    if (f == NULL)
        throw FileException("Cannot open file.");

	fscanf(f, "%d %d\n", &rows, &cols);
    max_rows = rows;
    max_cols = cols;
	data = new T[rows * cols];
	if (!data) error_alloc(__FILE__, __LINE__);
    fread(data, sizeof(T), rows * cols, f);
	fclose(f);
}

template <class T>
void Matrix<T>::mul(T *in, T *out) {
	for (int i = 0; i < rows; ++i) {
		T sum = (T)0.;
		for (int j = 0; j < cols; ++j) {
			sum += data[i * cols + j] * in[j];
		}
		out[i] = sum;
	}
};

template <class T>
void Matrix<T>::convert(const Matrix<Rgb> &m, Matrix<Float> &d) {
    for (int i = 0; i < m.get_rows(); ++i) {
        for (int j = 0; j < m.get_cols(); ++j) {
            d(i, j) = m(i, j).luminance();
        }
    }
}

template <class T>
void Matrix<T>::resize(int new_rows) {
    resize(new_rows, cols);
}

template <class T>
void Matrix<T>::resize(int new_rows, int new_cols) {
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
        data = new T[max_rows * max_cols];
    }
}

} // end namespace

#endif
