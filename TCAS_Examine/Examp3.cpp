#include "Exam3.h"

static const int TAU = 20;
static const double DMOD = 0.3 * 1852;
static const int ZTHR = 850;
static const int HMD = 0;

bool Vertical_RA(const double& s1, const double& v1)
{
	double cpatime = -s1 / v1;
	printf("cpatime: %f\n", cpatime);

	return abs(s1) < ZTHR || (cpatime > 0 && cpatime < TAU);
}

float dot(const vect2_t& p1, const vect2_t& p2)
{
	return p1.x * p2.x + p1.y * p2.y;
}

bool Horizontal_RA(const vect2_t& s2, const vect2_t& v2)
{
	double TAUmod = 0;
	if (s2.modsquared > 10)
		TAUmod = ((DMOD * DMOD) - s2.modsquared) / dot(s2, v2);

	printf("TAUmod: %f\n", TAUmod);

	return s2.mod < DMOD || (TAUmod > 0 && TAUmod <= TAU);
}

double Delta(const vect2_t s2, const vect2_t v2, double D)
{
	vect2_t Perpendicular(v2.y, -v2.x);

	return (D * D * v2.modsquared) - pow(dot(s2, Perpendicular), 2);
}

double root(const double& a, const double& b, const double& c, int epsilon)
{
	double n = ((b * b) - (4 * a * c));
	if (a == 0 || n < 0)
		return 0;
	return (sqrt(n) * epsilon - b) / (2 * a);
}

double Phi(const vect2_t& s2, const vect2_t& v2, int D, int epsilon)
{
	return root(v2.modsquared, 2 * dot(s2, v2), s2.modsquared - (D * D), epsilon);
}

bool CD2D(const vect2_t& s2, const vect2_t& v2, int D, int B)
{
	//if (v2.iszero() && s2.mod < D)
	//	return true;
	double delta = Delta(s2, v2, D);
	double phi = Phi(s2, v2, D, 1);
	//if (!v2.iszero() && delta >= 0 && phi >= B)
	if (delta >= 0 && phi >= B)
		return true;
	return false;
}

bool TCASII_RA(vect3_t& s3, vect3_t& v3)
{
	bool v = Vertical_RA(s3.z, v3.z);
	bool h = Horizontal_RA(s3.xy, v3.xy);
	return v && h && (HMD == 0 || CD2D(s3.xy, v3.xy, HMD, 0));
}

void Example3()
{
	vect3_t s(-2274258.3676679074, -5961969.6889963392, 974.30783268821199);
	vect3_t v(-2274252.1511302376, -5961953.3923390713, 1024.3630242545648);
	bool res = TCASII_RA(s, v);
}