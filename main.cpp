
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
    // apply p.b.c in x and y direction
    dr.x -= system_size_x_ * round(dr.x / system_size_x_);
    dr.y -= system_size_y_ * round(dr.y / system_size_y_);
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

  double dt = params.get_parameter<double>("dt");

  double verlet_list_radius =
		params.get_parameter<double>("verlet_list_radius");

  double A1 = params.get_parameter<double>("A1");
  double A2 = params.get_parameter<double>("A2");


  double zlim = params.get_parameter<double>("zlim");
  unsigned int number_of_bins =
		params.get_parameter<double>("number_of_bins");

  double equilibration_time = 
         params.get_parameter<double>("equilibration_time");
  double time_between_samples = 
         params.get_parameter<double>("time_between_samples");
  unsigned int number_of_samples = 
         params.get_parameter<unsigned int>("number_of_samples");

  double energy_scale = params.get_parameter<double>("energy_scale");
  double exponent = params.get_parameter<double>("exponent");
  double cut_off_radius = params.get_parameter<double>("cut_off_radius");
  double system_size_x = system_size_xy;
  double system_size_y = system_size_xy;

  Potential potential(energy_scale, exponent, cut_off_radius,
            system_size_x, system_size_y);

  SystemBD<Potential> system(seed, number_of_particles, system_size_xy,
					dt, verlet_list_radius, A1, potential);

  double area = system_size_xy * system_size_xy;
  Density rho_z(-zlim, zlim, number_of_bins, 'z', area);


  double time = 0;
  while (time < equilibration_time) {
    system.Integrate(time_between_samples);
    time += time_between_samples;
    cout << equilibration_time << '\t' << time << '\n';
  }
  cout << "Equilibration done\n";
  system.SetPotential(A2);
  system.ResetTime();

  rho_z.Sample(system.GetPositions());
  string density_name = "rhoz0.dat";
  rho_z.Save(density_name);
  rho_z.Reset();


  ofstream out_time; 
  out_time.open("time.dat");
  out_time << 0 << '\t' << 0 << '\n';

  for (unsigned int isample = 1; isample < number_of_samples; ++isample) {
    system.Integrate(time_between_samples);

    rho_z.Sample(system.GetPositions());
    string density_name = "rhoz" + to_string(isample) + ".dat";
    rho_z.Save(density_name);
    rho_z.Reset();

    cout << number_of_samples << '\t' << isample << endl;
    out_time << isample << '\t' << system.GetTime() << '\n'; 
  } 
  out_time.close();

	return 0;
}
