#include "Exam2.h"

vect3_t vect3_sub(vect3_t a, vect3_t b)
{
	vect3_t v;
	v.x = a.x - b.x;
	v.y = a.y - b.y;
	v.z = a.z - b.z;
	return v;
}

static inline vect3_t VECT3(double x, double y, double z)
{
	vect3_t v;
	v.x = x;
	v.y = y;
	v.z = z;
	return v;
}

static inline vect2_t VECT2(double x, double y)
{
	vect2_t v;
	v.x = x;
	v.y = y;
	return v;
}

vect3_t vect3_add(vect3_t a, vect3_t b)
{
	return (VECT3(a.x + b.x, a.y + b.y, a.z + b.z));
}

vect3_t vect3_scmul(vect3_t a, double b)
{
	return (VECT3(a.x * b, a.y * b, a.z * b));
}

vect2_t vect2_sub(vect2_t a, vect2_t b)
{
	return (VECT2(a.x - b.x, a.y - b.y));
}

double vect2_abs(vect2_t a)
{
	return (sqrt(pow(2, a.x) + pow(2, a.y)));
}

static void make_cpa(double d_t, vect3_t pos_a, vect3_t pos_b)
{
	double d_h = vect2_abs(vect2_sub(VECT2(pos_a.x, pos_a.y), VECT2(pos_b.x, pos_b.y)));
	double d_v = abs(pos_a.z - pos_b.z);
	printf("Examp2 v/h: %f - %f\n", d_v, d_h);
}

void Example2()
{
	double t_cpa;
	double x, y, z, x_i, y_i, z_i;
	double vx, vy, vz, vx_i, vy_i, vz_i;
	vect3_t vel, my_vel, rel_vel, rel_pos_3d, cpa_pos, my_cpa_pos;

	x = -2274258.3676679074;
	y = -5961969.6889963392;
	z = 974.30783268821199;
	vx = 80;
	vy = -8;
	vz = -103;

	x_i = -2274252.1511302376;
	y_i = -5961953.3923390713;
	z_i = 1024.3630242545648;
	vx_i = 100;
	vy_i = -9;
	vz_i = -102;

	rel_pos_3d = VECT3(x_i, y_i, z_i - z);
	vel = VECT3(vx_i, vy_i, vz_i);
	my_vel = VECT3(vx, vy, vz);;
	rel_vel = vect3_sub(vel, my_vel);

	double calc1 = (-(rel_pos_3d.x * rel_vel.x) - (rel_pos_3d.y * rel_vel.y) - (rel_pos_3d.z * rel_vel.z));
	double calc2 = (pow(2, rel_vel.x) + pow(2, rel_vel.y) + pow(2, rel_vel.z));

	t_cpa = floor(calc1 / calc2);

	printf("t_cpa: %f\n", t_cpa);

	vect3_t hedefPos;
	hedefPos.x = x_i;
	hedefPos.y = y_i;
	hedefPos.z = z_i;

	vect3_t ownPos;
	ownPos.x = x;
	ownPos.y = y;
	ownPos.z = z;

	cpa_pos = vect3_add(hedefPos, vect3_scmul(vel, t_cpa));
	my_cpa_pos = vect3_add(ownPos, vect3_scmul(my_vel, t_cpa));

	//printf("hedef: %.0f %.0f %.0f\n", cpa_pos.x, cpa_pos.y, cpa_pos.z);
	//printf("own: %.0f %.0f %.0f\n", my_cpa_pos.x, my_cpa_pos.y, my_cpa_pos.z);

	make_cpa(t_cpa, my_cpa_pos, cpa_pos);
}