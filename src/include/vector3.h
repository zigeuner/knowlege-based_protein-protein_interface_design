#ifndef __VECTOR3_H
#define __VECTOR3_H

#ifdef _USE_OLD_LIBRARIES
#include <math.h>
#else
#include <cmath>
#endif
#include <string>
#include "config.h"
#include "stdio.h"

class CVector3 {

protected:
  double x, y, z;

public:
  CVector3() 
    : x(0.0), y(0.0), z(0.0) { }
  CVector3(double a, double b, double c) 
    : x(a), y(b), z(c) { }
  CVector3(double a) { 
    x = a; y = a; z = a; 
  }
  CVector3(double* a) { 
    x = a[0]; 
    y = a[1]; 
    z = a[2]; 
  }
  CVector3(int* a) { 
    x = (double) a[0]; 
    y = (double) a[1]; 
    z = (double) a[2]; 
  }
  CVector3(const CVector3& v) 
    : x(v.x), y(v.y), z(v.z) { }

  double X() const {return (x);}
  double Y() const {return (y);}
  double Z() const {return (z);}
  double Length() const { 
    return sqrt(x*x + y*y + z*z); 
  }
  double SqLength() const { 
    return x*x + y*y + z*z;
  }
  double Dist(const CVector3& b) const { 
    CVector3  c = CVector3(x - b.x,y - b.y,z - b.z);
    return (c.Length());
  }
  double DistSq(const CVector3& b) const { 
    CVector3  c = CVector3(x - b.x,y - b.y,z - b.z);
    return (c.SqLength());
  }
  double Normalize() {
    double l = Length(); 

    x /= l; 
    y /= l; 
    z /= l;
    return l;
  }
  CVector3 DivideFloor(const CVector3& a, const CVector3& b) {
    x = floor(a.x/b.x);
    y = floor(a.y/b.y);
    z = floor(a.z/b.z);
    return *this;
  }
  CVector3 MultSubtract(const CVector3& a, const CVector3& b) {
    x -= a.x*b.x;
    y -= a.y*b.y;
    z -= a.z*b.z;
    return *this;
  }
  double Min() {
    double minimum = x;
    if (y < minimum) {minimum = y;}
    if (z < minimum) {minimum = z;}
    return minimum;
  }
  void MinImage(const CVector3& half_ell, const CVector3& ell) {
    if (x > half_ell.x) {
      x -= ell.x;
    } else if (x < -half_ell.x) {
      x += ell.x;
    }

    if (y > half_ell.y) {
      y -= ell.y;
    } else if (y < -half_ell.y) {
      y += ell.y;
    }

    if (z > half_ell.z) {
      z -= ell.z;
    } else if (z < -half_ell.z) {
      z += ell.z;
    }
  }

  void Int(int* array) {
    array[0] = (int) x; 
    array[1] = (int) y; 
    array[2] = (int) z; 
  }

  std::string Display() {
    char     temp[1000];
    std::string   line;
      
    sprintf(temp,"%+12.8f %+12.8f %+12.8f",x,y,z);
    line = temp;
    return(line);
  }

  std::string Display(const double scale) {
    char     temp[1000];
    std::string   line;
      
    sprintf(temp,"%+12.8f %+12.8f %+12.8f",x*scale,y*scale,z*scale);
    line = temp;
    return(line);
  }

  /* Functions declared in the .cpp file */
  double Angle(CVector3&, CVector3&);
  double TorsionAngle(CVector3&, CVector3&);

  /* Operator overloads */
  double& operator[](int i);
  CVector3& operator=(const CVector3& v) {
    x = v.x; 
    y = v.y; 
    z = v.z; 
    return *this;
  }
  CVector3& operator=(const double value) {
    x = value; 
    y = value; 
    z = value; 
    return *this;
  }
  void operator+=(const CVector3& v) {
    x += v.x; 
    y += v.y; 
    z += v.z;
  }
  void operator-=(const CVector3& v) {
    x -= v.x; 
    y -= v.y; 
    z -= v.z;
  }
  void operator/=(const double a) {
    x /= a; 
    y /= a; 
    z /= a;
  }
  void operator*=(const double a) {
    x *= a; 
    y *= a; 
    z *= a;
  }
  friend std::ostream& operator<<(std::ostream& os, const CVector3& v);
  friend std::istream& operator>>(std::istream& is, CVector3& v);
  friend CVector3 operator-(const CVector3& a, const CVector3& b) { 
    return CVector3(a.x - b.x, a.y - b.y, a.z - b.z);
  }
  friend CVector3 operator+(const CVector3& a, const CVector3& b) {
    return CVector3(a.x + b.x, a.y + b.y, a.z + b.z);
  }
  friend double operator*(const CVector3& a, const CVector3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }
  friend CVector3 operator^(const CVector3& a, const CVector3& b) {
    return CVector3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
  }
  friend CVector3 operator*(double a, const CVector3& b) {
    return CVector3(a*b.x, a*b.y, a*b.z);
  }
  friend CVector3 operator*(const CVector3& a, double b) {
    return CVector3(a.x*b, a.y*b, a.z*b);
  }
  friend CVector3 operator/(const CVector3& a, double b) {
    return CVector3(a.x/b, a.y/b, a.z/b);
  }
  friend CVector3 operator/(const CVector3& a, const CVector3& b) {
    return CVector3(a.x/b.x, a.y/b.y, a.z/b.z);
  }
};

#endif

