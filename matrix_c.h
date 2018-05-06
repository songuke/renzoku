#ifndef _MATRIX_C_H_
#define _MATRIX_C_H_

#include <cstdio>
#include <memory>
#include <typeinfo>
#include <exception>

#include "rgb.h"
#include "error.h"

using namespace std;

namespace Renzoku {

/** 
 * Matrix with column-major internal storage.
 */
template <class T>
class MatrixC {
public:
	MatrixC();
	MatrixC(const MatrixC<T> &m);

    /**
     * Defautl constructor of type T is called.
     * In case of int, float, double: what is the default value of int()?
     */
	MatrixC(int rows, int cols);
	
    MatrixC(int rows, int cols, T val);
	~MatrixC();

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
    void transpose(MatrixC<T> &t);

    /**
     * Matrix-vector multiplication.
     */
	void mul(const vector<T> &in, vector<T> &out);

    /**
     * Sum all values of a row.
     */
    void row_sum(vector<T> &out);

    /**
     * Convert 3-channel matrix into a flat matrix.
     */
    void flatten_magnitude(MatrixC<Float> &d);

    void flatten(MatrixC<Float> &d) const;

    void get_column(int c, vector<T> &col) const;
    void get_row(int r, vector<T> &row) const;

    void get_column_norm(vector<Float> &norm) const;

    void set(const T &val);

    /**
     * Setting a new number of rows.
     */
    void resize(int new_rows);

    /**
     * Resize the matrix to new rows/columns. 
     * 
     * The internal memory is reallocated if necessary. Old matrix data becomes invalid.
     */
    void resize(int new_rows, int new_cols);

    /**
     * Preallocate some memory but keep existing matrix size. 
     * Existing data will be invalidated.
     */
    void reserve(int max_rows);
    void reserve(int max_rows, int max_cols);

    /**
     * Extract sub-matrix using row or column indices.
     */
    void subrow(vector<int> row_indices, MatrixC<T> &mat) const;
    void subcol(vector<int> col_indices, MatrixC<T> &mat) const;

protected:
	T *data;
    int max_rows, max_cols;
	int rows, cols;
};

// template specialization for Rgb data type. These functions are stored in matrix_c.cc to avoid duplicate definition.
template <>
void MatrixC<Rgb>::load(const char *file);

template <>
void MatrixC<Rgb>::save(const char *file) const;

template <>
void MatrixC<Rgb>::flatten(MatrixC<Float> &mat) const;

template <class T>
MatrixC<T>::MatrixC() : data(NULL), rows(0), cols(0), max_rows(0), max_cols(0) {	

}

template <class T>
MatrixC<T>::MatrixC(int rows, int cols) : rows(rows), cols(cols), max_rows(rows), max_cols(cols) {
	data = new T[max_rows * max_cols];
	if (!data) error_alloc(__FILE__, __LINE__);
	for (int i = 0; i < rows * cols; ++i)
		data[i] = T();

}

template <class T>
MatrixC<T>::MatrixC(int rows, int cols, T val) : rows(rows), cols(cols), max_rows(rows), max_cols(cols) {
	data = new T[max_rows * max_cols];
	if (!data) error_alloc(__FILE__, __LINE__);
	for (int i = 0; i < rows * cols; ++i) 
		data[i] = val;
}

template <class T>
MatrixC<T>::MatrixC(const MatrixC<T> &m) : rows(m.rows), cols(m.cols), max_rows(m.max_rows), max_cols(m.max_cols) {
	data = new T[max_rows * max_cols];
	if (!data) error_alloc(__FILE__, __LINE__);
	memcpy(data, m.data, sizeof(T) * rows * cols);
}

template <class T>
MatrixC<T>::~MatrixC() {
	if (data)
		delete [] data;	
}

template <class T>
int MatrixC<T>::get_rows() const { return rows; }

template <class T>
int MatrixC<T>::get_cols() const { return cols; }

template <class T>
inline T MatrixC<T>::operator()(int i) const {
	return data[i];
}

template <class T>
inline T& MatrixC<T>::operator()(int i) {
	return data[i];
}

template <class T>
inline T MatrixC<T>::operator()(int i, int j) const {
	return data[j * rows + i];
}

template <class T>
inline T& MatrixC<T>::operator()(int i, int j) {
	return data[j * rows + i];
}

template <class T>
void MatrixC<T>::save(const char *file) const {
	FILE *f = fopen(file, "wb");
    if (f == NULL) 
        throw FileException("Cannot open file.");

	fprintf(f, "%d %d\n", rows, cols);    
    fwrite(data, sizeof(T), rows * cols, f);    
    fclose(f);
}

template <class T>
void MatrixC<T>::load(const char *file) {
	FILE *f = fopen(file, "rb");
    if (f == NULL)
        throw FileException("Cannot open file.");

	fscanf(f, "%d %d\n", &rows, &cols);
	max_rows = rows;
    max_cols = cols;
    data = new T[rows * cols];    
	if (!data) error_alloc(__FILE__, __LINE__);
    if (typeid(T) == typeid(Rgb)) {
		for (int j = 0; j < cols; ++j) {
			for (int i = 0; i < rows; ++i) {            
                Float v[3];
                fread(data, sizeof(Float), 3, f);

                data[j * rows + i] = Rgb(v[0], v[1], v[2]);
            }
        }    
    } else {
	    fread(data, sizeof(T), rows * cols, f);
    }
	fclose(f);
}

template <class T>
void MatrixC<T>::mul(const vector<T> &in, vector<T> &out) {
    out.resize(rows);
	for (int i = 0; i < rows; ++i) out[i] = 0.f;
	
	for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {
			out[i] += data[j * rows + i] * in[j];
		}		
	}
}

template <class T>
void MatrixC<T>::row_sum(vector<T> &out) {
    out.resize(rows);
    for (int i = 0; i < rows; ++i) out[i] = 0.f;

    for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {
			out[i] += data[j * rows + i];
		}		
	}
}

template <class T>
void MatrixC<T>::transpose(MatrixC<T> &out) {
    out.resize(cols, rows);

	for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {
			out(j, i) = (*this)(i, j);
		}		
	}
}

template <>
void MatrixC<Rgb>::flatten_magnitude(MatrixC<Float> &d);

template <class T>
void MatrixC<T>::resize(int new_rows) {
    resize(new_rows, cols);
}

template <class T>
void MatrixC<T>::resize(int new_rows, int new_cols) {
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

template <class T>
void MatrixC<T>::reserve(int max_rows) {    
    reserve(max_rows, cols);
}

template <class T>
void MatrixC<T>::reserve(int max_rows, int max_cols) {
    if (this->max_rows * this->max_cols <= max_rows * max_cols) {
        this->max_rows = max_rows;
        this->max_cols = max_cols;
    
        // reallocate
        if (data) delete [] data;
        data = new T[max_rows * max_cols];
    }
}

template <class T>
void MatrixC<T>::subrow(vector<int> row_indices, MatrixC<T> &mat) const {
    mat.resize(row_indices.size(), cols);
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < row_indices.size(); ++i) {
            mat(i, j) = (*this)(row_indices[i], j);
        }
    }
}

template <class T>
void MatrixC<T>::subcol(vector<int> col_indices, MatrixC<T> &mat) const {
    mat.resize(rows, col_indices.size());
    for (int j = 0; j < col_indices.size(); ++j) {
        for (int i = 0; i < rows; ++i) {
            mat(i, j) = (*this)(i, col_indices[j]);
        }
    }
}

template <class T>
void MatrixC<T>::set(const T& val) {
    for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {
			(*this)(i, j) = val;
		}		
	}
}

template <class T>
void MatrixC<T>::get_column(int c, vector<T> &col) const {
    col.resize(rows);
    for (int i = 0; i < rows; ++i)
        col[i] = (*this)(i, c);
}

template <class T>
void MatrixC<T>::get_row(int r, vector<T> &row) const {
    row.resize(cols);
    for (int i = 0; i < cols; ++i)
        row[i] = (*this)(r, i);
}

template <>
void MatrixC<Rgb>::get_column_norm(vector<Float> &norm) const; 

template <>
void MatrixC<Float>::get_column_norm(vector<Float> &norm) const; 

} // end namespace

#endif
