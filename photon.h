#ifndef _PHOTON_H_
#define _PHOTON_H_

#include "common.h"
#include "vec3.h"
#include "rgb.h"
#include "onb.h"

#include <queue>
using namespace std;

namespace Renzoku {
struct Photon {
    Rgb flux;       // power
    Vec3 p;         // location
    Vec3 wi;        // incoming direction    
    int id;

    Photon(Rgb flux, const Vec3 &p, const Vec3 &wi) {
        this->flux = flux;
        this->p = p;
        this->wi = wi;
    }
};

typedef vector<Photon*> Photons;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~add a new struct for ppm~~~~~~~~~~~~~~~~~~~~~
struct Hitpoint
{
	Vec3 p;			//hit position
	Vec3 pn;		//hit point normal
	Vec3 dir;		//ray direction
	Material* m;	//hit point material
	Onb uvn;
	int pixel_x;
	int pixel_y;
	Rgb pixel_weight;
	Float radius;		//current photon radius
	int accum_photons;	//accumulated photon count
	Rgb accum_flux;		//accumulated reflected flux
	//Float distance;		//distance between the hitpoint and a photon

	Hitpoint(const Vec3 &p, const Vec3 &pn, const Vec3 &dir, Material *m, int p_x, int p_y, Onb &uvn) {
		this->p = p;
		this->pn = pn;
		this->dir = dir;
		this->m = m;
		pixel_x = p_x;
		pixel_y = p_y;
		this->uvn = uvn;
		//this->radius = radius;
		//this->distance = -1;
	}

};
typedef vector<Hitpoint*> Hitpoints;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
static bool operator<(const Photon &p1, const Photon &p2) {
    // We want to put the photon with furthest distance to p to the top of the heap 
    // (so we can quickly discard it if needed).
    if (p1.dist == -1 || p2.dist == -1) {
        cout << "Error photon comparison." << endl;
    }
    return p1.dist < p2.dist;
}*/

/**
 * Functor to compare two "pointers". 
 * Priority given to "larger" pointer (the top of the heap is greatest).
 */
/*
template <class T>
struct ptr_less : public binary_function<T, T, bool> {
    bool operator()(const T &a, const T &b) {
        return (*a) < (*b); // this will call operator>() function of type T
    }
};
typedef priority_queue<Photon*, Photons, ptr_less<Photon*> > PhotonHeap;
*/

/*
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~for ppm~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static bool operator<(const Hitpoint &hp1, const Hitpoint &hp2) {
	if(hp1.distance == -1 || hp2.distance == -1) {
		cout << "Error hitpoint comparison." << endl;
	}
	return hp1.distance < hp2.distance;
}
typedef priority_queue<Hitpoint*, Hitpoints, ptr_less<Hitpoint*> > HitPointHeap;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
struct ptr_greater {
    bool operator()(const Photon *a, const Photon *b) {
        return a->dist > b->dist;
    }
};
typedef priority_queue<Photon*, Photons, ptr_greater> PhotonHeap;
*/

};

#endif 
