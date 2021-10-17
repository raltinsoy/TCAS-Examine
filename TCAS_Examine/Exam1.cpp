#include "Exam1.h"

#define m_to_nm 0.000539956803
#define nm_to_m 1852
#define PI 3.14159265359
#define m_to_ft 3.2808399
#define th_d_alt 0.5 // m

// calculates distance between 2 points
double Distance(double X0, double Y0, double Z0, double X1, double Y1, double Z1)
{
	return sqrt(pow(X1 - X0, 2) + pow(Y1 - Y0, 2) + pow(Z1 - Z0, 2));
}

// calculates range tau: time to CPA
double calculate_tau_r(double vx, double  vy, double vz, double vx_i, double vy_i, double vz_i, double distance)
{
	float closure_rate = sqrt(pow(vx - vx_i, 2) + pow(vy - vy_i, 2) + pow(vz - vz_i, 2));

	return (distance / closure_rate);
}

void wgs_to_geo(double x, double y, double z, float* lat_back, float* lon_back, float* alt_back)
{
	// WGS84 ellipsoid constants
	double a = 6378137; // radius
	double e = 8.1819190842622e-2;  // eccentricity

	double asq = pow(a, 2);
	double esq = pow(e, 2);

	double b = sqrt(asq * (1 - esq));
	double bsq = pow(b, 2);
	double ep = sqrt((asq - bsq) / bsq);
	double p = sqrt(pow(x, 2) + pow(y, 2));
	double th = atan2(a * z, b * p);

	double lon = atan2(y, x);
	double lat = atan2((z + pow(ep, 2) * b * pow(sin(th), 3)), (p - esq * a * pow(cos(th), 3)));
	double N = a / (sqrt(1 - esq * pow(sin(lat), 2)));
	double alt = p / cos(lat) - N;

	// mod lat to 0-2pi
	lon = lon - floor((lon + PI) / (2 * PI)) * (2 * PI);

	// correction for altitude near poles left out.

	*lon_back = lon;
	*lat_back = lat;
	*alt_back = alt;
}

void wgs_to_enu(double x, double y, double z, double* x_enu, double* y_enu, double* z_enu, float lat, float lon)
{
	*x_enu = x * -sin(lon) + y * cos(lon);
	*y_enu = x * -cos(lon) * sin(lat) + y * -sin(lon) * sin(lat) + z * cos(lat);
	*z_enu = x * cos(lon) * cos(lat) + y * sin(lon) * cos(lat) + z * sin(lat);
}

// returns vertical tau: time to CPA on vertical plane
double calculate_tau_vert(double vx, double  vy, double vz, double vx_i, double vy_i, double vz_i, double lat, double lon, double lat_i, double lon_i, double d_alt, double* bearing, double* v_speed)
{
	double vx_enu, vy_enu, vz_enu, vx_enu_i, vy_enu_i, vz_enu_i;
	double mod, bearing_o, bearing_i;

	wgs_to_enu(vx, vy, vz, &vx_enu, &vy_enu, &vz_enu, lat, lon);

	wgs_to_enu(vx_i, vy_i, vz_i, &vx_enu_i, &vy_enu_i, &vz_enu_i, lat_i, lon_i);

	// printf("%f,%f,%f\n",vx_enu_i,vy_enu_i,vz_enu_i);

	bearing_o = atan(vy_enu / vx_enu) * 180 / PI;
	bearing_i = atan(vy_enu_i / vx_enu_i) * 180 / PI;

	*bearing = bearing_i - bearing_o;
	*v_speed = vz_enu_i;

	if (abs(-vz_enu + vz_enu_i) < th_d_alt)
	{
		return 0.0;
	}

	return abs(d_alt / (-vz_enu + vz_enu_i));
}

#define _TORAD (PI/180)

void geo_to_wgs(float lat, float lon, float alt, double* x, double* y, double* z)
{
	// Input em graus
	lon = lon * _TORAD;
	lat = lat * _TORAD;

	// WGS84 ellipsoid constants
	double a = 6378137; // radius
	double b = 0.00669438;  // eccentricity

	double x_aux, y_aux, z_aux;
	double N;

	N = a * a / (sqrt(a * a * cos(lat) * cos(lat) + b * b * sin(lat) * sin(lat)));

	x_aux = (N + alt) * cos(lon) * cos(lat);
	y_aux = (N + alt) * sin(lon) * cos(lat);
	z_aux = ((N * b * b / (a * a)) + alt) * sin(lat);

	*x = x_aux;
	*y = y_aux;
	*z = z_aux;
}

void get_thresold(int* th_tau_TA, int* th_tau_RA, int* th_alt_TA, int* th_alt_RA, int alt)
{
	if (alt < 1000)
	{
		*th_tau_TA = 20;
		*th_tau_RA = 0; // N/A
		*th_alt_TA = 850;
		*th_alt_RA = 0; // N/A
	}

	else if (alt < 2350)
	{
		*th_tau_TA = 25;
		*th_tau_RA = 15;
		*th_alt_TA = 850;
		*th_alt_RA = 600;
	}

	else if (alt < 5000)
	{
		*th_tau_TA = 30;
		*th_tau_RA = 20;
		*th_alt_TA = 850;
		*th_alt_RA = 600;
	}

	else if (alt < 10000)
	{
		*th_tau_TA = 40;
		*th_tau_RA = 25;
		*th_alt_TA = 850;
		*th_alt_RA = 600;
	}

	else if (alt < 20000)
	{
		*th_tau_TA = 45;
		*th_tau_RA = 30;
		*th_alt_TA = 850;
		*th_alt_RA = 600;
	}

	else if (alt < 42000)
	{
		*th_tau_TA = 48;
		*th_tau_RA = 35;
		*th_alt_TA = 850;
		*th_alt_RA = 700;
	}

	else
	{
		*th_tau_TA = 48;
		*th_tau_RA = 35;
		*th_alt_TA = 1200;
		*th_alt_RA = 800;
	}
}

bool lla2xyz(double xyz[3], const double lla[3]) {

	if ((lla[0] < -90.0) || (lla[0] > +90.0) || (lla[1] < -180.0) || (lla[1] > +360.0)) {
		std::cout << "WGS lat or WGS lon out of range" << std::endl;
		return 0;
	}

	double A_EARTH = 6378137.0;
	double flattening = 1.0 / 298.257223563;
	double NAV_E2 = (2.0 - flattening) * flattening; // also e^2
	double deg2rad = PI / 180.0;

	double slat = sin(lla[0] * deg2rad);
	double clat = cos(lla[0] * deg2rad);
	double r_n = A_EARTH / sqrt(1.0 - NAV_E2 * slat * slat);
	xyz[0] = (r_n + lla[2]) * clat * cos(lla[1] * deg2rad);
	xyz[1] = (r_n + lla[2]) * clat * sin(lla[1] * deg2rad);
	xyz[2] = (r_n * (1.0 - NAV_E2) + lla[2]) * slat;

	return 1;
}

void matrixMultiply(double C[3][3], const double A[3][3], const double B[3][3]) {
	C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
	C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
	C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];
	C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
	C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
	C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];
	C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
	C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
	C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}

//------------------------------------------------------------------------------------------------
// WgsConversions::matrixMultiply [Private]  --- Multiply 3x3 matrix times a 3x1 vector c=Ab
//------------------------------------------------------------------------------------------------
void matrixMultiply(double c[3], const double A[3][3], const double b[3]) {
	c[0] = A[0][0] * b[0] + A[0][1] * b[1] + A[0][2] * b[2];
	c[1] = A[1][0] * b[0] + A[1][1] * b[1] + A[1][2] * b[2];
	c[2] = A[2][0] * b[0] + A[2][1] * b[1] + A[2][2] * b[2];
}

void rot(double R[3][3], const double angle, const int axis) {
	double cang = cos(angle * PI / 180);
	double sang = sin(angle * PI / 180);

	if (axis == 1) {
		R[0][0] = 1;
		R[0][1] = 0;
		R[0][2] = 0;
		R[1][0] = 0;
		R[2][0] = 0;
		R[1][1] = cang;
		R[2][2] = cang;
		R[1][2] = sang;
		R[2][1] = -sang;
	}
	else if (axis == 2) {
		R[0][1] = 0;
		R[1][0] = 0;
		R[1][1] = 1;
		R[1][2] = 0;
		R[2][1] = 0;
		R[0][0] = cang;
		R[2][2] = cang;
		R[0][2] = -sang;
		R[2][0] = sang;
	}
	else if (axis == 3) {
		R[2][0] = 0;
		R[2][1] = 0;
		R[2][2] = 1;
		R[0][2] = 0;
		R[1][2] = 0;
		R[0][0] = cang;
		R[1][1] = cang;
		R[1][0] = -sang;
		R[0][1] = sang;
	}
}

void rot3d(double R[3][3], const double reflat, const double reflon) {

	double R1[3][3], R2[3][3];

	rot(R1, 90 + reflon, 3);
	rot(R2, 90 - reflat, 1);

	matrixMultiply(R, R2, R1);
}

bool xyz2enu(double enu[3], const double xyz[3], const double ref_lla[3]) {

	double ref_xyz[3], diff_xyz[3], R[3][3];

	// First, calculate the xyz of reflat, reflon, refalt
	if (!lla2xyz(ref_xyz, ref_lla))
		return 0;

	//Difference xyz from reference point
	diff_xyz[0] = xyz[0] - ref_xyz[0];
	diff_xyz[1] = xyz[1] - ref_xyz[1];
	diff_xyz[2] = xyz[2] - ref_xyz[2];

	rot3d(R, ref_lla[0], ref_lla[1]);

	matrixMultiply(enu, R, diff_xyz);

	return 1;
}

bool lla2enu(double enu[3], const double lla[3], const double ref_lla[3]) {

	double xyz[3];

	if (!lla2xyz(xyz, lla))
		return 0;

	if (!xyz2enu(enu, xyz, ref_lla))
		return 0;

	return 1;
}

void Example1()
{
	double distance, tau_range;
	double x, y, z, x_i, y_i, z_i;
	double vx, vy, vz, vx_i, vy_i, vz_i;
	float lat, lon, alt, lat_i, lon_i, alt_i;
	double x_enu, y_enu, z_enu;
	double bearing2, d_alt, bearing, v_speed;
	int th_tau_TA, th_tau_RA, th_alt_TA, th_alt_RA;

	vx = -22.5373;
	vy = -48.2408;
	vz = 206.0114;

	lat = 19.1095;
	lon = -111.2722;
	alt = 10450.3152;

	double lla[3];
	double lla_i[3];
	double xyz[3];
	double xyz_i[3];
	double enu[3];
	//double ref_lla[3];
	//ref_lla[0] = 0;
	//ref_lla[1] = 0;
	//ref_lla[2] = 0;

	lla[0] = lat;
	lla[1] = lon;
	lla[2] = alt;

	//lla2enu(enu, lla, ref_lla);

	lla2xyz(xyz, lla);
	x = xyz[0];
	y = xyz[1];
	z = xyz[2];

	//geo_to_wgs(lat, lon, alt, &x, &y, &z);

	get_thresold(&th_tau_TA, &th_tau_RA, &th_alt_TA, &th_alt_RA, alt * m_to_ft);
	printf("Tau: %i,%i alt: %i,%i\n", th_tau_TA, th_tau_RA, th_alt_TA, th_alt_RA);

	vx_i = -48.2587;
	vy_i = -0.0913;
	vz_i = -172.8109;

	lat_i = 19.0819;
	lon_i = -111.2680;
	alt_i = 10428.0056;

	//geo_to_wgs(lat_i, lon_i, alt_i, &x_i, &y_i, &z_i);
	lla_i[0] = lat_i;
	lla_i[1] = lon_i;
	lla_i[2] = alt_i;

	lla2xyz(xyz_i, lla_i);
	x_i = xyz_i[0];
	y_i = xyz_i[1];
	z_i = xyz_i[2];

	distance = Distance(x, y, z, x_i, y_i, z_i);
	printf("Distance: %f nm, %f meter\n", distance * m_to_nm, distance);

	tau_range = calculate_tau_r(vx, vy, vz, vx_i, vy_i, vz_i, distance);
	printf("Tau range: %f\n", tau_range);

	/// bearing 2
	double dx_enu, dy_enu, dz_enu;
	wgs_to_enu((x_i - x), (y_i - y), (z_i - z), &dx_enu, &dy_enu, &dz_enu, lat, lon);
	double normD = sqrt(pow(dx_enu, 2) + pow(dy_enu, 2));

	double vx_enu, vy_enu, vz_enu;
	wgs_to_enu(vx, vy, vz, &vx_enu, &vy_enu, &vz_enu, lat, lon);
	double normV = sqrt(pow(vx_enu, 2) + pow(vy_enu, 2));

	if (normD * normV != 0)
	{
		bearing2 = acos((dx_enu * vx_enu + dy_enu * vy_enu) / (normD * normV)) * 180.0 / PI;
		int side = (vx_enu * (y_i - y) - vy_enu * (x_i - x));
		if (side > 0) {
			bearing2 = -bearing2;
		}
	}

	double tau_vertical = 0;

	d_alt = alt_i - alt;
	if (abs(d_alt) > th_d_alt) {
		tau_vertical = calculate_tau_vert(vx, vy, vz, vx_i, vy_i, vz_i, lat, lon, lat_i, lon_i, d_alt, &bearing, &v_speed);
	}

	if ((abs(d_alt * m_to_ft) < th_alt_RA) && (tau_vertical < th_tau_RA) && (tau_range < th_tau_RA))
	{
		printf("RA\n");
	}
	else if ((abs(d_alt) * m_to_ft < th_alt_TA) && (tau_vertical < th_tau_TA) && (tau_range < th_tau_TA))
	{
		printf("TA\n");
	}
	else if ((abs(d_alt) * m_to_ft < 1200) || distance * m_to_nm < 6)
	{
		printf("PT\n");
	}
	else
	{
		printf("OT\n");
	}
}