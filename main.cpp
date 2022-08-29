
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
    Vec3 dr = r1 - r2;
    // periodic boundary conditions if system_size_..._ > 0
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

vector<Vec3> ReadPositions(string data_name)
{
  vector<Vec3> positions;

  std::ifstream in;
  in.open(data_name);

  string line;

  double x,y,z;
  while (getline(in, line)) {
    std::stringstream ss(line);
    ss >> x;
    ss >> y;
    ss >> z;
    positions.push_back(Vec3(x, y, z));
  }

  return positions;
}

double AverageDistanceMoved(const vector<Vec3>& p1, const vector<Vec3>& p2, double Lx)
{
   double distance = 0;
   for (unsigned int i = 0; i < p1.size(); ++i) {
      Vec3 dr = p1[i] - p2[i];
      distance += dr.Length();
   }
   return distance / p1.size();
}


int main()
{
  cout << "Start Reading parameters\n" << flush;
  Config params("input.txt");
  unsigned long int seed =
		params.get_parameter<unsigned long int>("seed");

  //unsigned int number_of_particles =
  //	params.get_parameter<unsigned int>("number_of_particles");

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

  //unsigned int number_of_particles = 
  //     params.get_parameter<unsigned int>("number_of_particles");

  double energy_scale = params.get_parameter<double>("energy_scale");
  double exponent = params.get_parameter<double>("exponent");
  double cut_off_radius = params.get_parameter<double>("cut_off_radius");

  Potential potential(energy_scale, exponent, cut_off_radius,
            system_size_x, system_size_y, system_size_z);


  vector<Vec3> initial_positions = ReadPositions("eqpos.dat");

  SystemBD<Potential> system(seed, initial_positions, system_size_x,
					system_size_y, system_size_z, dt,
					verlet_list_radius, A1, potential);

  //SystemBD<Potential> system(seed, number_of_particles, system_size_x,
	//				system_size_y, system_size_z, dt,
	//				verlet_list_radius, A1, potential);

  double area = system_size_x * system_size_y;
  Density rho_z(-zlim, zlim, number_of_bins, 'z', area);

  cout << "Start Equilibration\n" << flush;

  double time = 0;
  while (time < equilibration_time) {
    system.Integrate(time_between_samples);
    time += time_between_samples;
    cout << equilibration_time << '\t' << time << '\n' << flush;
  }
  string positions_name = "equilibrium_positions.dat";
  system.SavePositions(positions_name);

  vector<Vec3> eq_positions = system.GetPositions();
  cout << "Average Distance between init and eq:"
       << AverageDistanceMoved(initial_positions, eq_positions, system_size_x)
       << endl;


  cout << "Equilibration done\n" << flush;

  system.SetPotential(A2);
  system.ResetTime();

  rho_z.Sample(system.GetPositions());
  string density_name = "rhoz0.dat";
  rho_z.Save(density_name);
  rho_z.Reset();


  ofstream out_time; 
  out_time.open("time.dat");
  out_time << 0 << '\t' << 0 << '\n';

  cout << "Start Sampling\n" << flush;
  for (unsigned int isample = 1; isample < number_of_samples + 1; ++isample) {
    system.Integrate(time_between_samples);

    rho_z.Sample(system.GetPositions());

    string density_name = "rhoz" + to_string(isample) + ".dat";
    rho_z.Save(density_name);

    rho_z.Reset();

    cout << number_of_samples << '\t' << isample << endl << flush;
    out_time << isample << '\t' << system.GetTime() << '\n'; 
  } 
  out_time.close();

  positions_name = "final_positions.dat";
  system.SavePositions(positions_name);

  	return 0;
}
