/*
 To do:
	- fix initialization
    -- check GetPositions
  -- increase read precision
*/


#ifndef GUARD_SYSTEM_BD_H
#define GUARD_SYSTEM_BD_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <boost/random.hpp>

#include "vec3.h"

namespace systemBD_helper {

// calculate square distance between r1 and r2,
// apply periodic boundaires if L > 0
double distance_squared(Vec3 r1, Vec3 r2, double Lx, double Ly, double Lz)
{
     
	r1 -= r2;
  if (Lx > 0) r1.x -= Lx * round(r1.x/Lx);
	if (Ly > 0) r1.y -= Ly * round(r1.y/Ly);
	if (Lz > 0) r1.z -= Lz * round(r1.z/Lz);
	return r1.LengthSquared();
}

double distance(Vec3 r1, Vec3 r2, double Lx, double Ly, double Lz)
{
	r1 -= r2;
  if (Lx > 0) r1.x -= Lx * round(r1.x/Lx);
	if (Ly > 0) r1.y -= Ly * round(r1.y/Ly);
	if (Lz > 0) r1.z -= Lz * round(r1.z/Lz);
	return r1.Length();
}

};

template <class Potential>
class SystemBD {
 public:
	SystemBD(unsigned long int seed,
		   unsigned int number_of_particles,
		   double system_size_x,
		   double system_size_y,
		   double system_size_z,
		   double dt,
		   double verlet_list_radius,
		   double A,
           Potential potential);

  void MakeTimeStep(double dt);
  void Integrate(double delta_t);

  void SavePositions(std::string name) const;		    
  void ReadPositions(std::string name);

  void SetPotential(double newA) { A_ = newA; }
  double getPotential() const { return A_; }

  std::vector<Vec3> GetPositions() const { return positions_; }

  long unsigned int GetNumberOfVerletListUpdates() const
		{ return number_of_verlet_list_updates_; }

  void ResetTime() { time_ = 0; }
  double GetTime() const { return time_;}

 private:
  const boost::random::uniform_int_distribution<int>
                                       random_int_distribution_;
  const boost::normal_distribution<double> normal_distribution_;

  boost::mt19937 random_number_generator_;
  boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<double> > random_normal_distribution_;

	// initialize the particles on a square lattice
	void RandomInit();

	void UpdateVerletList();

  void UpdateForces();

	// private variable

	unsigned int number_of_particles_;

	// system size in the x, y and z directions
	double system_size_x_;
	double system_size_y_;
	double system_size_z_;
	
    // time step size
    double dt_;

	// Radius of for the Verlet list
	double verlet_list_radius_;

	// particle positions
	std::vector<Vec3> positions_;
	// particle positions at when the Verlet list was last updated
	std::vector<Vec3> positions_at_last_update_;

	// Verlet list
	std::vector<std::vector<unsigned int> > verlet_list_;
	// number of neighbors in the Verlet list
	std::vector<unsigned int> number_of_neighbors_;

	// when distance between position_[i] and position_at_last_update_[i]
	// is larger than max_diff_, the Verlet list needs to be updated
	double max_diff_;	

	// externale potential: U(z) = A_ * z * z
	double A_;

    // potential
    Potential potential_;

	// keep track of number of Verlet list updates
	unsigned long int number_of_verlet_list_updates_;

    // current time	
    double time_;

   // forces
    std::vector<Vec3> forces_;
};

template <class Potential>
SystemBD<Potential>::SystemBD(
	unsigned long int seed,
	unsigned int number_of_particles,
	double system_size_x,
	double system_size_y,
	double system_size_z,
  double dt,
	double verlet_list_radius,
	double A,
    Potential potential)
  : random_int_distribution_(0, number_of_particles - 1),
    normal_distribution_(0.0,1.0),
	random_number_generator_(seed),
    random_normal_distribution_(random_number_generator_,
                                normal_distribution_),
	number_of_particles_(number_of_particles),
	system_size_x_(system_size_x),
	system_size_y_(system_size_y),
	system_size_z_(system_size_z),
	dt_(dt),
	verlet_list_radius_(verlet_list_radius),
    positions_(number_of_particles_),
    positions_at_last_update_(number_of_particles_),
	verlet_list_(number_of_particles_,
				std::vector<unsigned int>(number_of_particles)),
	number_of_neighbors_(number_of_particles_),
	A_(A),
    potential_(potential),
	number_of_verlet_list_updates_(0),
    time_(0.0),
    forces_(number_of_particles_)
{
  max_diff_ = (verlet_list_radius - potential.GetCutOffRadius()) / 2;

  RandomInit();
  UpdateVerletList();
}

template<class Potential>
void SystemBD<Potential>::RandomInit()
{

  std::vector<Vec3> lattice_positions;

  double dx = 1.25;
  double dy = dx;
  double dz = 1.1;
  int n_per_xy = floor(system_size_x_ / dx);

  Vec3 temp;
  int iz = 0;
  while (lattice_positions.size() < number_of_particles_) {
    int ix = 0;
    int iy = 0;
    while (iy < n_per_xy) {
      if (iz == 0) {
        temp.x = ix * dx;
        temp.y = iy * dy;
        temp.z = 0.0;
        lattice_positions.push_back(temp);
      } else {
        temp.x = ix * dx;
        temp.y = iy * dy;
        temp.z = iz * dz;
		if ( (iz % 2) == 1) {
		  temp.x += dx * sqrt(3)/2;
		  temp.y += dy * sqrt(3)/2;
        }
        lattice_positions.push_back(temp);

        temp.z = -iz * dz;
        lattice_positions.push_back(temp);
      }

      ix += 1;
      if (ix == n_per_xy) {
        ix = 0;
        iy += 1;
      }
    }
    iz += 1;
  }
  // shuffle positions

  for (unsigned int i = 0; i < lattice_positions.size(); ++i) {
    unsigned int j = random_int_distribution_(random_number_generator_);
	temp = lattice_positions[i];
    lattice_positions[i] = lattice_positions[j];
    lattice_positions[j] = temp;
  }

  for (unsigned int i = 0; i < number_of_particles_; ++i){
    positions_[i] = lattice_positions[i];
  }
}

template <class Potential>
void SystemBD<Potential>::Integrate(double delta_time)
{
  while (delta_time > dt_) {
    MakeTimeStep(dt_);
    delta_time -= dt_;
  }
  MakeTimeStep(delta_time);
}

template <class Potential>
void SystemBD<Potential>::MakeTimeStep(double dt)
{
  UpdateForces();

  //double lf = 0;
  //for (unsigned int i=0; i < number_of_particles_; ++i) {
    //if ( lf < forces_[i].Length()) lf = forces_[i].Length();
  //}
  //std::cout << lf << std::endl;

  double sqrt_2_dt = sqrt(2 * dt);
  bool update_verlet_list = false;
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_[i].z -= 2.0 * positions_[i].z * A_ * dt;
    positions_[i]   += forces_[i] * dt;

    // add Brownian displacement
    positions_[i].x += sqrt_2_dt *
                random_normal_distribution_();
    positions_[i].y += sqrt_2_dt *
                random_normal_distribution_();
    positions_[i].z += sqrt_2_dt *
                random_normal_distribution_();

    double dist = systemBD_helper::distance_squared(positions_[i],
                  positions_at_last_update_[i], system_size_x_,
                  system_size_y_, system_size_z_);
    if (dist > max_diff_ * max_diff_) update_verlet_list = true;
  }

  time_ += dt;
  if (update_verlet_list) UpdateVerletList();
  
}

template <class Potential>
void SystemBD<Potential>::UpdateForces()
{
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    forces_[i] *= 0;
  }

  Vec3 force_ij;
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    for (unsigned int nj = 0; nj < number_of_neighbors_[i]; ++nj) {
      unsigned int j = verlet_list_[i][nj]; // neighbor index
      force_ij = potential_.Force(positions_[i], positions_[j]);
      forces_[i] += force_ij;
      forces_[j] -= force_ij;
    }
  }
}

template <class Potential>
void SystemBD<Potential>::SavePositions(std::string name) const
{
  std::ofstream out;
  out.open(name);
  out <<std::setprecision(16);
  for (unsigned int i = 0; i < number_of_particles_; ++i) {
	out << positions_[i].x << '\t'
        << positions_[i].y << '\t'		
        << positions_[i].z << '\n';
  }

	out.close();
}

template <class Potential>
void SystemBD<Potential>::ReadPositions(std::string name)
{
  std::ifstream in;
  in.open(name);
  std::string line_str;

  positions_ = std::vector<Vec3>();
  double x, y, z;
  while (getline(in, line_str) ){ 
    std::cout << line_str << std::endl;
    std::stringstream ss(line_str);
    ss >> x;
    ss >> y;
    ss >> z;
    std::cout << x << '\t' << y << '\t' << z << '\n';
    //positions_.push_back(Vec3(x, y, z));
  }

  number_of_particles_ = positions_.size();

  in.close();
}



template <class Potential>
void SystemBD<Potential>::UpdateVerletList()
{
  number_of_verlet_list_updates_ += 1;

  std::fill(number_of_neighbors_.begin(),
		    number_of_neighbors_.end(), 0);

  for (unsigned int i = 0; i < number_of_particles_; ++i) {
    positions_at_last_update_[i] = positions_[i];
    for (unsigned int j = i + 1; j < number_of_particles_; ++j) {
      if (systemBD_helper::distance_squared(positions_[i],
				                    positions_[j], system_size_x_,
                                    system_size_y_, system_size_z_)
			< verlet_list_radius_ * verlet_list_radius_ ) {
        verlet_list_[i][ number_of_neighbors_[i] ] = j;
        ++number_of_neighbors_[i];
        //verlet_list_[j][ number_of_neighbors_[j] ] = i;
        //++number_of_neighbors_[j];
      }
    }
  }

}

#endif
