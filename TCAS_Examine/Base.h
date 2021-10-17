#pragma once

struct vect2_t {
	double	x;
	double	y;

	float mod;
	float modsquared;
	float bearing;

	vect2_t() {}

	vect2_t(double _x, double _y) :
		x(_x), y(_y) {
		updateangular();
	}

	void updateangular() {
		modsquared = (float(x) * float(x)) + (float(y) * float(y));
		mod = sqrt(modsquared);
		bearing = atan2(y, x);
	}
};

struct vect3_t {
	double	x;
	double	y;
	double	z;

	vect2_t xy;

	vect3_t() {}

	vect3_t(double _x, double _y, double _z) :
		xy(_x, _y), z(_z) {
	}
};
