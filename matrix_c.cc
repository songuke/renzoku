#include "matrix_c.h"

namespace Renzoku {

template <>
void MatrixC<Rgb>::save(const char *file) const {
	FILE *f = fopen(file, "wb");
    if (f == NULL) 
        throw FileException("Cannot open file.");

	fprintf(f, "%d %d\n", rows, cols);
    
	for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {            
			Float v[3] = { data[j * rows + i].red(), 
						   data[j * rows + i].green(), 
						   data[j * rows + i].blue() };
			fwrite(v, sizeof(Float), 3, f);
		}
	}
    
	fclose(f);
}

template <>
void MatrixC<Rgb>::load(const char *file) {
    FILE *f = fopen(file, "rb");
    if (f == NULL)
        throw FileException("Cannot open file.");

	fscanf(f, "%d %d\n", &rows, &cols);
	max_rows = rows;
    max_cols = cols;
    data = new Rgb[rows * cols];    
	if (!data) error_alloc(__FILE__, __LINE__);
    
	for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {            
            Float v[3];
            fread(data, sizeof(Float), 3, f);

            data[j * rows + i] = Rgb(v[0], v[1], v[2]);
        }
    }        
	fclose(f);
}

template <>
void MatrixC<Rgb>::flatten_magnitude(MatrixC<Float> &d) {
    d.resize(rows, cols);
    for (int j = 0; j < cols; ++j) {
		for (int i = 0; i < rows; ++i) {        
            d(i, j) = sqrt((*this)(i, j).r * (*this)(i, j).r + 
                           (*this)(i, j).g * (*this)(i, j).g +
                           (*this)(i, j).b * (*this)(i, j).b);
        }
    }
}


template <>
void MatrixC<Rgb>::flatten(MatrixC<Float> &d) const {
    d.resize(3 * rows, cols);
    for (int j = 0; j < cols; ++j) {
        for (int i = 0; i < rows; ++i) {
            d(3 * i,     j) = (*this)(i, j).r;
            d(3 * i + 1, j) = (*this)(i, j).g;
            d(3 * i + 2, j) = (*this)(i, j).b;
        }
    }
}

template <>
void MatrixC<Rgb>::get_column_norm(vector<Float> &norms) const {
    norms.resize(cols);
    
    for (int c = 0; c < cols; ++c) {
        Float magnitude = 0.f;

        for (int r = 0; r < rows; ++r) {
            Rgb v = (*this)(r, c);
            magnitude += v.red() * v.red() + v.green() * v.green() + v.blue() * v.blue();                        
        }    

        norms[c] = sqrt(magnitude);
    }
}

template <>
void MatrixC<Float>::get_column_norm(vector<Float> &norms) const {
    norms.resize(cols);
    
    for (int c = 0; c < cols; ++c) {
        Float magnitude = 0.f;

        for (int r = 0; r < rows; ++r) {            
            magnitude += (*this)(r, c) * (*this)(r, c);                        
        }

        norms[c] = sqrt(magnitude);
    }
}

} // end namespace
