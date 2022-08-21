
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
            double system_size_x, double system_size_y)
  : energy_scale_(energy_scale), exponent_(exponent),
    cut_off_radius_(cut_off_radius), system_size_x_(system_size_x),
    system_size_y_(system_size_y)
  {}

  Vec3 Force(const Vec3& r1, const Vec3& r2) const
  {
    Vec3 dr = r2 - r1;
    double dr_length = dr.Length();
    if ( dr_length < cut_off_radius_ ) {
        double f = pow(dr_length, - (exponent_ - 1));
		f *= energy_scale_ * exponent_ / dr_length;
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
};

int main()
{
  Config params("input.txt");
  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  unsigned int number_of_particles =
		params.get_parameter<unsigned int>("number_of_particles");

  double system_size_xy =
		params.get_parameter<double>("system_size_xy");

  double dt = params.get_parameter<double>("max_mc_step_size");

  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");

  double A = params.get_parameter<double>("A");


  double zlim = params.get_parameter<double>("zlim");
  unsigned int number_of_bins =
		params.get_parameter<double>("number_of_bins");

  double energy_scale = 100;
  double exponent = 12;
  double cut_off_radius = 1.01;
  double system_size_x = system_size_xy;
  double system_size_y = system_size_xy;

  Potential potential(energy_scale, exponent, cut_off_radius,
            system_size_x, system_size_y);

  SystemBD<Potential> system(seed, number_of_particles, system_size_xy,
					dt, verlet_list_radius, A, potential);

  double area = system_size_xy * system_size_xy;
  Density rho_z(-zlim, zlim, number_of_bins, 'z', area);


	return 0;
}
