/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef PRIMITIVERECOVERY_HPP_
#define PRIMITIVERECOVERY_HPP_

#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

namespace PrimitiveRecovery
{
// Primitive to conservative variables
template <class data_t>
Tensor<1, data_t> PtoC(const double p)
{
    Tensor<1, data_t> q;
    const double eps = p[0];
    const double ux  = p[1];
    const double uy  = p[2];
    const double nn  = p[3];
    double ux2 = ux*ux;
    double uy2 = uy*uy;
    double ut  = sqrt(1.+ux2+uy2);

    // Ttt                                                                                                         
    q[0] = (eps*(3.+4.*(ux2+uy2)))/3.;
    // Ttx                                                                                                         
    q[1] = (4.*ux*ut*eps)/3.;
    // Tty                                                                                                         
    q[2] = (4.*uy*ut*eps)/3.;
    // Jt                                                                                                          
    q[3] = nn*ut;

    return q;
}

// Conservative to primitive variables
  template <class data_t>
Tensor<1, data_t> CtoP(const double q)
{
    Tensor<1, data_t> p;
    const double Ttt = q[0];
    const double Ttx = q[1];
    const double Tty = q[2];
    const double Jt  = q[3];
    double Ttt2 = Ttt*Ttt;
    double Ttx2 = Ttx*Ttx;
    double Tty2 = Tty*Tty;
    double sqrt1 = sqrt(4.*Ttt2-3.*(Ttx2+Tty2));
    double sqrt2 = sqrt(2.*Ttt2-3.*(Ttx2+Tty2)+Ttt*sqrt1);
    double ut;

    // eps
    p[0] = -Ttt + sqrt1;
    // ux
    p[1] = (3.*Ttx)/(2.*sqrt2);
    // uy
    p[2] = (3.*Tty)/(2.*sqrt2);
    // ut
    ut = sqrt(1.+p[1]*p[1]+p[2]*p[2]);
    p[3] = Jt/ut;
    return p;
}
} 
#endif /* PRIMITIVERECOVERY_HPP_ */
