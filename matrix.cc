#include "matrix.h"

namespace Renzoku {

template <>
void Matrix<Rgb>::save(const char *file) const {
	FILE *f = fopen(file, "wb");
    if (f == NULL) 
        throw FileException("Cannot open file.");

	fprintf(f, "%d %d\n", rows, cols);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Float v[3] = { data[i * cols + j].red(), 
                            data[i * cols + j].green(), 
                            data[i * cols + j].blue() };
            fwrite(v, sizeof(Float), 3, f);
        }
    }    
    fclose(f);
}

template <>
void Matrix<Rgb>::load(const char *file) {
    FILE *f = fopen(file, "rb");
    if (f == NULL)
        throw FileException("Cannot open file.");

	fscanf(f, "%d %d\n", &rows, &cols);
    max_rows = rows;
    max_cols = cols;
	data = new Rgb[rows * cols];
	if (!data) error_alloc(__FILE__, __LINE__);
    
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Float v[3];
            fread(data, sizeof(Float), 3, f);

            data[i * cols + j] = Rgb(v[0], v[1], v[2]);
        }
    }    
	fclose(f);
}

} // end namespace
