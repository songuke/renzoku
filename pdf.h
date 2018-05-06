#ifndef _PROBABILITY_DISTRIBUTION_FUNCTION_H_
#define _PROBABILITY_DISTRIBUTION_FUNCTION_H_

#include "common.h"

#include <vector>
using namespace std;

namespace Renzoku {

/**
 * Probability distribution function 
 */
class Pdf {
};

/**
 * Discrete probability distribution function given by an array
 */
class DiscretePdf : public Pdf {
public:
    DiscretePdf();
	DiscretePdf(const vector<Float> &arr);
	
    void sample(Float rand, int &k, Float &pdf_k) const;
    void sample(Random &rd, int &k, Float &pdf_k) const;
    
    static void sample(Float rand, Float *pdf, Float *cdf, int n, int &k, Float &pdf_k);
    static void sample(Float rand, vector<Float> pdf, vector<Float> cdf, int &k, Float &pdf_k);
    static void sample(Float rand, Float *cdf, int n, int &k);
    static void sample(Float rand, vector<Float> cdf, int &k);

    /**
     * Return the probability of element k.
     */
    Float probability(int k) const;

    void set_distribution(const vector<Float> &arr);
    
    /**
     * Return true when all probabilities are zero.
     */
    bool is_valid() const;

protected:
	vector<Float> pdf;
    vector<Float> cdf;

    bool is_zero;
};

/**
 * Discrete probability distribution function given by 2D array.
 */
class DiscretePdf2D : public Pdf {
public:
    DiscretePdf2D();
	DiscretePdf2D(const ImageFloat &arr);
	    
    /**
     * Return row and column.
     */
    void sample(Random &rd, int &row, int &col, Float &pdf) const;    
    Float probability(int row, int col) const;

    void set_distribution(const ImageFloat &arr);

    bool is_valid() const;

protected:
	const ImageFloat *dist;     /// the raw distribution, not normalized.
    Float dist_total;
    
    DiscretePdf dist_row;       /// pdf for sampling a row.

    ImageFloat *cdf_col;        /// cdf for sampling a column given a row.    
};

} // end namespace

#endif