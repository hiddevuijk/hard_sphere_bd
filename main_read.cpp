
#include "vec3.h"
#include "config_file.h"
#include "density.h"
#include "systemBD.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

class Potential {
 public:
  Potential(double energy_scale, double exponent, double cut_off_radius,
            double system_size_x, double system_size_y,
            double system_size_z)
  : energy_scale_(energy_scale), exponent_(exponent),
    cut_off_radius_(cut_off_radius), system_size_x_(system_size_x),
    system_size_y_(system_size_y), system_size_z_(system_size_z)
  {}

  Vec3 Force(const Vec3& r1, const Vec3& r2) const
  {
    Vec3 dr = r2 - r1;
	if (system_size_x_ > 0) {
      dr.x -= system_size_x_ * round(dr.x / system_size_x_);
	}
    if (system_size_y_ > 0) {
      dr.y -= system_size_y_ * round(dr.y / system_size_y_);
    }
    if (system_size_z_ > 0) {
      dr.z -= system_size_z_ * round(dr.z / system_size_z_);
    }

    double dr_length = dr.Length();
    if ( dr_length < cut_off_radius_ ) {

      double f = pow(1 / dr_length, exponent_ + 2);
	  f *=  energy_scale_ * exponent_ ;

	  return dr * f;

	} else {
      return dr * 0;
    }
  }

 double GetCutOffRadius() const { return cut_off_radius_; } 

 private:
  double energy_scale_;
  double exponent_;
  double cut_off_radius_;
  double system_size_x_;
  double system_size_y_;
  double system_size_z_;
};

int main()
{
  cout << "Start Reading parameters\n" << flush;
  Config params("input.txt");
  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  unsigned int number_of_particles =
		params.get_parameter<unsigned int>("number_of_particles");

  double system_size_x =
		params.get_parameter<double>("system_size_x");
  double system_size_y =
		params.get_parameter<double>("system_size_y");
  double system_size_z =
		params.get_parameter<double>("system_size_z");

  double dt = params.get_parameter<double>("dt");

  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");

  double A1 = params.get_parameter<double>("A1");
  //double A2 = params.get_parameter<double>("A2");


  double zlim = params.get_parameter<double>("zlim");
  unsigned int number_of_bins =
		params.get_parameter<double>("number_of_bins");

  double equilibration_time = 
         params.get_parameter<double>("equilibration_time");
  double time_between_samples = 
         params.get_parameter<double>("time_between_samples");
  //unsigned int number_of_samples = 
         //params.get_parameter<unsigned int>("number_of_samples");

  double energy_scale = params.get_parameter<double>("energy_scale");
  double exponent = params.get_parameter<double>("exponent");
  double cut_off_radius = params.get_parameter<double>("cut_off_radius");

  Potential potential(energy_scale, exponent, cut_off_radius,
            system_size_x, system_size_y, system_size_z);

  SystemBD<Potential> system(seed, number_of_particles, system_size_x,
					system_size_y, system_size_z, dt,
					verlet_list_radius, A1, potential);

  double area = system_size_x * system_size_y;
  Density rho_z(-zlim, zlim, number_of_bins, 'z', area);

   
  string positions_name = "equilibrium_positions.dat";
  system.ReadPositions(positions_name);

  cout << "Reading done\n" << flush;
  cout << "Start Equilibration\n" << flush;

  double time = 0;
  while (time < equilibration_time) {
    system.Integrate(time_between_samples);
    time += time_between_samples;
    //cout << equilibration_time << '\t' << time << '\n' << flush;
  }

  positions_name = "equilibrium2_positions.dat";
  system.SavePositions(positions_name);

  cout << "Equilibration done\n" << flush;

	return 0;
}