#include "registrar.h" 
#include "geometry_collision.h"
#include "state_collision.h"

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
}
