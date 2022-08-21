/*

  3 dimensional vector class.
  operators (+,-, +=, -=, *, / ...) are defined element wise

*/



#ifndef GUARD_VEC3_H
#define GUARD_VEC3_H

#include <math.h>
#include <iostream>

class Vec3 {
 public:
  Vec3(double x, double y, double z)
    : x(x), y(y), z(z) {}

  Vec3()
    : x(0.0), y(0.0), z(0.0) {}
  
  double x, y, z;

  double Length() const {
    return sqrt(x * x + y * y + z * z);
  }
  double LengthSquared() const {
    return x * x + y * y + z * z;
  }

  // set length to d
  void Normalize(double d=1.0) {
    double l = Length();
    x *= d / l;
    y *= d / l;
    z *= d / l;
  }


  // dot product
  double DotProduct(const Vec3& right) const{
    return x * right.x + y * right.y + z * right.z;
  }

  Vec3 CrossProduct(const Vec3& right) const {
    return Vec3(y * right.z - z * right.y,
                z * right.x - x * right.z,
                x * right.y - y * right.x);
  }

  // operators
  Vec3 operator+=(const Vec3& r) {
    x += r.x; y += r.y; z += r.z;
    return *this;
  }

  Vec3 operator+=(const double& add) {
    x += add; y += add; z += add;
    return *this;
  }

  Vec3 operator-=(const Vec3& r) {
    x -= r.x; y -= r.y; z -= r.z;
    return *this;
  }  

  Vec3 operator-=(const double& minus) {
    x -= minus; y -= minus; z -= minus;
    return *this;
  }

  Vec3 operator*=(const Vec3& r) {
    x *= r.x; y *= r.y; z *= r.z;
    return *this;
  }

  Vec3 operator*=(const double& mult) {
    x *= mult; y *= mult; z *= mult;
    return *this;
  }

  Vec3 operator/=(const Vec3& r) {
    x /= r.x; y /= r.y; z /= r.z;
    return *this;
  }

  Vec3 operator/=(const double& div) {
    x /= div; y /= div; z /= div;
    return *this;
  }

};


////////////////////////////////////////
//
// Nonmember functions
//
////////////////////////////////////////

std::ostream& operator<<(std::ostream& out_stream, const Vec3& r) {
  out_stream << r.x << '\t' << r.y << '\t' << r.z;
  return out_stream;
}


// math operators, defined elementwise
Vec3 operator+(const Vec3& left, const Vec3& right) {
  return Vec3(left.x + right.x, left.y + right.y, left.z + right.z);
}

Vec3 operator+(const Vec3& left, const double& right) {
  return Vec3(left.x + right, left.y + right, left.z + right);
}

Vec3 operator+(const double& left, const Vec3& right) {
  return Vec3(left + right.x, left + right.y, left + right.z);
}

Vec3 operator-(const Vec3& left, const Vec3& right) {
  return Vec3(left.x - right.x, left.y - right.y, left.z - right.z);
}

Vec3 operator-(const Vec3& left, const double& right) {
  return Vec3(left.x - right, left.y - right, left.z - right);
}

Vec3 operator-(const double& left, const Vec3& right) {
  return Vec3(left - right.x,left - right.y, left - right.z);
}

Vec3 operator*(const Vec3& left, const Vec3& right) {
  return Vec3(left.x * right.x , left.y * right.y, left.z * right.z);
}

Vec3 operator*(const Vec3&left, const double& right) {
  return Vec3(left.x * right, left.y * right, left.z * right);
}

Vec3 operator*(const double& left, const Vec3& right) {
  return Vec3(left * right.x, left * right.y, left * right.z);
}

Vec3 operator/(const Vec3& left, const Vec3& right) {
  return Vec3(left.x / right.x , left.y / right.y, left.z / right.z);
}

Vec3 operator/(const Vec3&left, const double& right) {
  return Vec3(left.x / right, left.y / right, left.z / right);
}

Vec3 operator/(const double& left, const Vec3& right) {
  return Vec3(left / right.x, left / right.y, left / right.z);
}

#endif
