#include "registrar.h" 
#include "geometry_collision.h"
#include "state_collision.h"
#include "geometry_1d.h"
#include "state_1d.h"


namespace 
{
	
	GeometryRegistrar<Ball> r1("ball");
	GeometryRegistrar<Disk> r2("disk");	
	
	StateRegistrar<GaussianPressureState> s1("gauss_pressure");

	// for the 2d collision simulation
	GeometryRegistrar<DiskLeft> r3("disk_left");
	GeometryRegistrar<DiskRight> r4("disk_right");
	StateRegistrar<LeftUniformVelocityState> s2("left_uniform_velocity");
	StateRegistrar<RightUniformVelocityState> s3("right_uniform_velocity");

	StateRegistrar<UniformVelocityState> s4("uniform_velocity");


	// for 1d
	GeometryRegistrar<Line> r5("line");
	StateRegistrar<GaussianPressure1DState> s5("1d_gauss_pressure");
}
