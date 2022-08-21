#ifndef GUARD_DENSITY_H
#define GUARD_DENSITY_H

#include "vec3.h"

#include <vector>
#include <fstream>
#include <math.h>


class Density {
 public:
  Density(double xmin, double xmax, unsigned int number_of_bins, char xyz, double area);
  void Sample(const std::vector<Vec3>& positions);

  // returns a vector with the current density
  std::vector<double> GetDensity() const;
  // returns a vector with the bin centers
  std::vector<double> GetBins() const { return bin_centers_; }

  long unsigned int GetNumberOfExcludedSamples() const
	{ return number_of_excluded_samples_; }

  long unsigned int GetNumberOfSamples() const
	{ return number_of_samples_; }

  // saves to  file called name
  void Save(std::string name);

 private:

  // range of the density
  double xmin_, xmax_;
  std::vector<double> density_;
  std::vector<double> bin_centers_;
  char xyz_;

  // area perpendicular to xyz
  double area_;

  // bin size
  double dx_;

  // total number of particles sampled
  long unsigned int number_of_samples_;

  // if a position is outside xmin_ -- xmax_
  // the samples is excluded
  long unsigned int number_of_excluded_samples_;
};

Density::Density(double xmin, double xmax, unsigned int number_of_bins, char xyz, double area)
  : xmin_(xmin), xmax_(xmax), density_(number_of_bins, 0.0),
	bin_centers_(number_of_bins), xyz_(xyz), area_(area),
    dx_( (xmax_ - xmin_) / number_of_bins ),
    number_of_samples_(0), number_of_excluded_samples_(0)
{
  for (unsigned int i = 0; i < number_of_bins; ++i) {
	bin_centers_[i] = xmin_ + (i + 0.5) * dx_;
  }
}



void Density::Sample(const std::vector<Vec3>& positions)
{


  for(unsigned int i = 0; i < positions.size(); ++i) {
    unsigned int bin_index = -1;
    if (xyz_ == 'x') { 
      bin_index = floor((positions[i].x - xmin_) / dx_);
    } else if (xyz_ == 'y') {
      bin_index = floor((positions[i].y - xmin_) / dx_);
    } else if (xyz_ == 'z') {
      bin_index = floor((positions[i].z - xmin_) / dx_);
    } else {
		// error
	} 

	if (bin_index < density_.size() ) {
		density_[bin_index] += 1.0;	
	} else {
		// particle outside domain xmin -- xmax
		number_of_excluded_samples_ += 1;
	}	
  }

  //number_of_samples_ += positions.size();
  number_of_samples_ += 1;
}

std::vector<double> Density::GetDensity() const
{
  std::vector<double> rho(density_.size());
  double norm = 1.0 / ( number_of_samples_ * dx_ * area_);
  for (unsigned int i = 0; i < rho.size(); ++i) {
	rho[i] = density_[i] * norm;	
  }
  return rho;
}

void Density::Save(std::string name)
{
  std::ofstream out;
  out.open(name);
  double norm = 1.0 / ( number_of_samples_ * dx_ * area_);
  for (unsigned int i = 0; i < density_.size(); ++i) {
    out << bin_centers_[i] << '\t' << density_[i] * norm;
    if (i < density_.size() - 1) out << '\n';
  }
  out.close();
}
#endif
