/*
   3-D vector class. 
   -S. Brown 05.22.01, 
   additions from L. Clark starting 01.01.03
*/

#include <iostream>
#include "include/vector3.h"

using std::ostream;
using std::istream;
using std::cout;
using std::endl;

/*------------------------------------------------------------------------*/
/* Calculate an angle between three position vectors a-b-c */
/* NOTE: the object itself is vector a */
/*     b    */
/*    / \   */
/*   /   \  */
/*  a     c */
/* cos(theta) = (a-b)/|a-b| .dot. (c-b)/|c-b| */
/* Requires:  b -- vector 2 */
/*            c -- vector 3 */
/*------------------------------------------------------------------------*/
double CVector3::Angle(CVector3& b, CVector3& c) {
  CVector3  a = CVector3(x,y,z), v1, v2;
  double    costheta,theta;

  v1 = a - b;
  v1 /= v1.Length();
  v2 = c - b;
  v2 /= v2.Length();
  costheta = v1*v2;
  theta = acos(costheta);

  return(theta);
}

/*------------------------------------------------------------------------*/
/* Calculate a torsion angle using three vectors */
/*  Where phi is the angle between the two plane normals    */
/*        i+1       i-1                                     */
/*          o         o       plane1: defined by bonds a, b */
/*         / \       /        plane2: defined by bonds b, c */
/*      c /   \ b   / a                                     */
/*       /     \   /                                        */
/*      /       \ /                                         */
/* i+2 o       i o                                          */
/* */
/* In the 'trans' conformation above, the angle is zero because the two */
/* plane normals are parallel.  In the 'gauche' conformation, it is Pi */
/* NOTE: the object itself is vector a defined as r_i - r_i-1 */
/* Requires:  b -- vector 2 (r_i+1 - r_i) */
/*            c -- vector 3 (r_i+2 - r_i+1) */
/*------------------------------------------------------------------------*/
double CVector3::TorsionAngle(CVector3& b, CVector3& c) {
  CVector3  a = CVector3(x,y,z);
  double    a_a,a_b,a_c,b_b,b_c,c_c;
  double    T1,T2,T3,T2T3_invroot;
  double    costheta,theta;

  /* calculate the dot products */
  a_a =	a*a;
  a_b =	a*b;
  a_c =	a*c;
  b_b =	b*b;
  b_c =	b*c;
  c_c =	c*c;

  /* calculate T1,T2,T3, terms that will have to used later */
  T1 = a_b*b_c - a_c*b_b;
  T2 = a_a*b_b - a_b*a_b;
  T3 = b_b*c_c - b_c*b_c;
  if (fabs(T2*T3) < 1.0e-11) {
    /* rare occurence of colinear placement of 4 atoms? */
    cout << "WARNING: (vector3.TorsionAngle) could not calculate angle" << endl;
    cout << "  rare colinear placement of 4 atoms?" << endl;
    return(0.0);
  }
  T2T3_invroot = 1.0/sqrt(T2*T3);

  /* Theta=0, when atom.1 and atom.4 are farthest from each other. */
  /* cos(theta) is defined so that the above condn on theta is satisfied */
  costheta = -(T1*T2T3_invroot);
  theta = acos(costheta);

  return(theta);
}

/*  good test, should give ~180 degree torsion angle
	calpha1 = CVector3(71.627,17.234,7.813);
	calpha2 = CVector3(71.610,16.157,6.908);
	cbeta1 = CVector3(72.384,16.168,5.735);
	cbeta2 = CVector3(73.199,17.281,5.457);
	a = cbeta1 - calpha1;
	b = calpha2 - cbeta1;
	c = cbeta2 - calpha2;
*/

ostream& operator<<(ostream& os, const CVector3& v)
{
    os << v.x << " " << v.y << " " << v.z ;
    return os;
}

istream& operator>>(istream& is, CVector3& v)
{
    is >> v.x >> v.y >> v.z;
    return is;
}

double& CVector3::operator[](int i)
{
  switch(i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default :
        cout << " Inappropriate vector component: " << i << "." << endl;
        exit(1);
  }
  return x;
}

