#include "pdf.h"
#include "random.h"
#include "image.h"

namespace Renzoku {

// ----------------------------------------------------------------------------
// 1D distribution
// ----------------------------------------------------------------------------

DiscretePdf::DiscretePdf() : is_zero(true) {
}

DiscretePdf::DiscretePdf(const vector<Float> &arr) {    
    set_distribution(arr);
}

void DiscretePdf::set_distribution(const vector<Float> &arr) {
    is_zero = false;

    // vector<Float>::size_type to fix signed/unsigned warning during compilation
    // http://stackoverflow.com/questions/275853/acceptable-fix-for-majority-of-signed-unsigned-warnings
    if (arr.size() == (vector<Float>::size_type)0) {
        is_zero = true;
        return;
    }

    // normalize to create a probability distribution function
	Float sum = 0.f;
	for (vector<Float>::size_type i = 0; i < arr.size(); ++i)
		sum += arr[i];
	
    pdf.resize(arr.size());
    cdf.resize(arr.size());
    
    if (sum == 0.0f) {
        // not sampleable
        for (vector<Float>::size_type i = 0; i < arr.size(); ++i) {
		    pdf[i] = 0.0f;		
            cdf[i] = 0.0f;
        }
        is_zero = true;

    } else {
	    for (vector<Float>::size_type i = 0; i < arr.size(); ++i)
		    pdf[i] = arr[i] / sum;		

        cdf[0] = pdf[0];
	    for (vector<Float>::size_type i = 1; i < pdf.size(); ++i)
		    cdf[i] = cdf[i - 1] + pdf[i];
        cdf[pdf.size() - 1] = 1.0f;
    }
}

bool DiscretePdf::is_valid() const {
    return ! is_zero;
}

void DiscretePdf::sample(Random &rd, int &k, Float &pdf_k) const {
	Float v = rd();    // FIXME: only generate random sample when the distribution is valid.
	
    // TODO: use binary search here

    vector<Float>::size_type i;
	for (i = 0; i < cdf.size(); ++i) {
		if (v < cdf[i]) {
			k = i;
			pdf_k = pdf[i];
			return;
		}
	}
	
    // no sampleable cases
	k = -1;
	pdf_k = 0.0f;
	return;
}

void DiscretePdf::sample(Float rand, int &k, Float &pdf_k) const {
    vector<Float>::size_type i;
	for (i = 0; i < cdf.size(); ++i) {
		if (rand < cdf[i]) {    // no = here as cdf[i] = 0 should indicate not sampleable (rand can have value 0).
			k = i;
			pdf_k = pdf[i];
			return;
		}
	}
	
	k = -1;
	pdf_k = 0.0f;
	return;
}

Float DiscretePdf::probability(int k) const {
	if (k >= 0 && k < pdf.size())
    	return pdf[k];
    return 0.0f;
}

void DiscretePdf::sample(Float rand, vector<Float> pdf, vector<Float> cdf, int &k, Float &pdf_k) {
	vector<Float>::size_type i;
	for (i = 0; i < cdf.size(); ++i) {
		if (rand < cdf[i]) {
			k = i;
			pdf_k = pdf[i];
			return;
		}
	}
	
	k = -1;
	pdf_k = 0.0f;
	return;
}

void DiscretePdf::sample(Float rand, vector<Float> cdf, int &k) {
    vector<Float>::size_type i;
	for (i = 0; i < cdf.size(); ++i) {
		if (rand < cdf[i]) {
			k = i;			
			return;
		}
	}
	
	k = -1;
	return;
}

void DiscretePdf::sample(Float rand, Float *pdf, Float *cdf, int n, int &k, Float &pdf_k) {
	int i = 0;
    for (i = 0; i < n; ++i) {
		if (rand < cdf[i]) {
			k = i;
			pdf_k = pdf[i];
			return;
		}
	}
	
	k = -1;
	pdf_k = 0.0f;
	return;
}

void DiscretePdf::sample(Float rand, Float *cdf, int n, int &k) {
    int i = 0;
    for (i = 0; i < n; ++i) {
		if (rand < cdf[i]) {
			k = i;			
			return;
		}
	}
	
	k = -1;
	return;
}
    

// ----------------------------------------------------------------------------
// 2D distribution
// ----------------------------------------------------------------------------

DiscretePdf2D::DiscretePdf2D() : dist(NULL), dist_total(0.f), cdf_col(NULL) {

}
	
DiscretePdf2D::DiscretePdf2D(const ImageFloat &arr) 
    : dist(NULL), dist_total(0.f), cdf_col(NULL) {

    set_distribution(arr);
}
	    
void DiscretePdf2D::sample(Random &rd, int &row, int &col, Float &pdf) const {
    Vec2 uv;
    uv.random(rd);

    int r, c;
    Float pdf_r;
    dist_row.sample(uv[0], r, pdf_r);

    int n = dist->get_width();    // should be pitch if aligned
    Float *local_cdf = cdf_col->get_data() + r * n;        
    DiscretePdf::sample(uv[1], local_cdf, n, c);
    
    row = r;
    col = c;
    pdf = (*dist)(r, c, 0) / dist_total;
}

Float DiscretePdf2D::probability(int row, int col) const {
    if (row < 0 || row >= dist->get_height()) return 0.f;
    if (col < 0 || col >= dist->get_width()) return 0.f;

    return (*dist)(row, col, 0) / dist_total;
}

void DiscretePdf2D::set_distribution(const ImageFloat &arr) {
    dist = &arr;
    dist_total = 0.f;

    // sum over row entries to build distribution of rows
    Float *data = dist->get_data();
    int height = dist->get_height();
    int width = dist->get_width();

    // build cdf
    vector<Float> row_sum(height);

    if (cdf_col != NULL) delete cdf_col;
    cdf_col = new ImageFloat(0.f, height, width);

    for (int i = 0; i < height; ++i) {

        // for row sampling
        row_sum[i] = 0.f;
        for (int j = 0; j < width; ++j) {
            row_sum[i] += (*dist)(i, j, 0);
        }

        // for normalization
        dist_total += row_sum[i];

        // for column sampling given a row
        Float sum = 0.f;
        for (int j = 0; j < width; ++j) {
            sum += (*dist)(i, j, 0);
            
            (*cdf_col)(i, j) = sum / row_sum[i];            
        }
    }
    
    dist_row.set_distribution(row_sum);

    // TODO: consider the pixel size
}

bool DiscretePdf2D::is_valid() const {
    return (dist_total > 0.0f);
}

} // end namespace