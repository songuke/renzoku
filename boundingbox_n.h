#ifndef _BOUNDINGBOX_N_H_
#define _BOUNDINGBOX_N_H_

namespace Renzoku {

template <int dim>
struct VecN {
    Float p[dim];

    VecN() {
    }

    VecN(const VecN &v) {
        memcpy(p, v.p, sizeof(Float) * dim);
    }

    VecN &operator=(const VecN &v) {
        memcpy(p, v.p, sizeof(Float) * dim);
        return *this;
    }

    Float &operator[](int index) {
        return p[index];
    }

    Float operator[](int index) const {
        return p[index];
    }

    VecN operator+(const VecN &v) const {
        VecN u;
        for (int i = 0; i < dim; ++i)
            u.p[i] = p[i] + v.p[i];
        return u;
    }

    VecN operator-(const VecN &v) const {
        VecN u;
        for (int i = 0; i < dim; ++i)
            u.p[i] = p[i] - v.p[i];
        return u;
    }

    VecN operator*(const Float &v) const {
        VecN u;
        for (int i = 0; i < dim; ++i)
            u.p[i] = p[i] * v;
        return u;
    }

    Float squared_length() const {
        Float sum = 0.0f;
        for (int i = 0; i < dim; ++i)
            sum += p[i] * p[i];
        return sum;
    }
};


/**
* N-dimensional bounding box
*/
template <int dim>
struct BoundingBoxN {
    VecN<dim> vmin, vmax;

    BoundingBoxN() {
        for (int i = 0; i < dim; ++i) {
            vmin.p[i] = FLT_MAX;
            vmax.p[i] = -FLT_MAX;
        }
    }

    int get_longest_axis() const {
        VecN<dim> axes = vmax - vmin;

        int index = 0;
        Float longest = 0.0f;
        for (int i = 0; i < dim; ++i) {
            if (fabs(axes.p[i]) > longest) {
                index = i;
                longest = fabs(axes.p[i]);
            }
        }
        return index;
    }

    void merge(const VecN<dim> &p) {
        for (int i = 0; i < dim; ++i) {
            vmin.p[i] = std::min(vmin.p[i], p.p[i]);
            vmax.p[i] = std::max(vmax.p[i], p.p[i]);
        }
    }

    VecN<dim> clamp(const VecN<dim> &p, const VecN<dim> &a, const VecN<dim> &b) const {
        VecN<dim> u;
        for (int i = 0; i < dim; ++i) {
            u.p[i] = std::min(std::max(p[i], a[i]), b[i]);
        }
        return u;
    }

    Float min_squared_distance(const VecN<dim> &p) const {
        VecN<dim> clamped_p = clamp(p, vmin, vmax);
        return (clamped_p - p).squared_length();
    }
};

}  // end namespace

#endif
