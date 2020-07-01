/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COORDINATETRANSFORMATIONS_HPP_
#define COORDINATETRANSFORMATIONS_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"

namespace CoordinateTransformations
{

  // Jacobian transformation matrix 
  template <class data_t>
  static Tensor<2, data_t> jacobian(const data_t x, const double y,
					 const double z)
  {
    // calculate useful position quantities
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);
    data_t r2 = simd_max(x * x + y * y + z * z, 1e-12);
    data_t r = sqrt(r2);

    // And the sines and cosines of phi and theta
    data_t cos_phi = x / rho;
    data_t sin_phi = y / rho;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac;
    jac[0][0] = x / r;
    jac[1][0] = cos_phi * z / r2;
    jac[2][0] = -y / rho2;
    jac[0][1] = y / r;
    jac[1][1] = sin_phi * z / r2;
    jac[2][1] = x / rho2;
    jac[0][2] = z / r;
    jac[1][2] = -rho / r2;
    jac[2][2] = 0.0;
    return jac;
  }

  // Incerse Jacobian
  template <class data_t>
  static Tensor<2, data_t> inverse_jacobian(const data_t x, const double y,
					    const double z)
  {
    // calculate useful position quantities
    data_t rho2 = simd_max(x * x + y * y, 1e-12);
    data_t rho = sqrt(rho2);
    data_t r2 = simd_max(x * x + y * y + z * z, 1e-12);
    data_t r = sqrt(r2);

    // And the sines and cosines of phi and theta
    data_t sin_theta = rho / r;
    data_t cos_phi = x / rho;
    data_t sin_phi = y / rho;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac;
    inv_jac[0][0] = x / r;
    inv_jac[1][0] = y / r;
    inv_jac[2][0] = z / r;
    inv_jac[0][1] = z * cos_phi;
    inv_jac[1][1] = z * sin_phi;
    inv_jac[2][1] = -r * sin_theta;
    inv_jac[0][2] = -y;
    inv_jac[1][2] = x;
    inv_jac[2][2] = 0.0;
    return inv_jac;
  }

  // Convert a Tensor (with two lower indices) in spherical coords to cartesian
  // coords
  template <class data_t>
  static Tensor<2, data_t>
  spherical_to_cartesian_LL(Tensor<2, data_t> spherical_g, data_t x, double y,
			    double z)
  {
    Tensor<2, data_t> cartesian_g;

    // derivatives for jacobian matrix - drdx etc
    Tensor<2, data_t> jac = jacobian(x, y, z);

    // Convert the Tensor to cartesian coords
    FOR2(i, j)
    {
        cartesian_g[i][j] = 0;
        FOR2(k, m)
        {
            cartesian_g[i][j] += spherical_g[k][m] * jac[k][i] * jac[m][j];
        }
    }
    return cartesian_g;
}
  
  // Convert a Tensor (with two lower indices) in cartesian coords to spherical
  // coords
  template <class data_t>
  static Tensor<2, data_t>
  cartesian_to_spherical_LL(Tensor<2, data_t> cartesian_g, data_t x, double y,
			    double z)
  {
    Tensor<2, data_t> spherical_g;

    // derivatives for inverse jacobian matrix - drdx etc    
    Tensor<2, data_t> inv_jac = inverse_jacobian(x, y, z);

    // Convert the Tensor to spherical coords
    FOR2(i, j)
      {
	spherical_g[i][j] = 0;
	FOR2(k, m)
	  {
	    spherical_g[i][j] += cartesian_g[k][m] * inv_jac[k][i] * inv_jac[m][j];
	  }
      }
    return spherical_g;
  }

  // Convert a vector (with one upper index) in spherical coords to cartesian
  // coords
  template <class data_t>
  Tensor<1, data_t> spherical_to_cartesian_U(Tensor<1, data_t> spherical_v,
                                           data_t x, double y, double z)
  {
    Tensor<1, data_t> cartesian_v;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> inv_jac = inverse_jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR1(i)
    {
        cartesian_v[i] = 0.0;
        FOR1(j) { cartesian_v[i] += inv_jac[i][j] * spherical_v[j]; }
    }
    return cartesian_v;
  }

  // Convert a vector (with one upper index) in cartesian coords to spherical  
  // coords
  template <class data_t>
  Tensor<1, data_t> cartesian_to_spherical_U(Tensor<1, data_t> cartesian_v,
					     data_t x, double y, double z)
  {
    Tensor<1, data_t> spherical_v;

    // derivatives for inverse jacobian matrix - drdx etc
    Tensor<2, data_t> jac = jacobian(x, y, z);

    // transform the vector to cartesian coords
    FOR1(i)
    {
      spherical_v[i] = 0.0;
      FOR1(j) { spherical_v[i] += jac[i][j] * cartesian_v[j]; }
    }
    return spherical_v;
  }

  // The area element of a sphere
  template <class data_t>
  data_t area_element_sphere(Tensor<2, data_t> spherical_g, data_t x, double y,
				double z)
  {

    // Normal vector s^i in spherical coords (spatial)
    Tensor<1, data_t> si_spher;
    si_spher[0] = 1.0/sqrt(spherical_g[0][0]);
    si_spher[1] = 0.0;
    si_spher[2] = 0.0;

    // Projection operator for the surface of a sphere P^i_j
    Tensor<2, data_t> Proj_spher;
    FOR2(i, j)
      {
	Proj_spher[i][j] =  TensorAlgebra::delta(i, j);
	FOR1(k)
	{
	  Proj_spher[i][j] += -spherical_g[i][k] * si_spher[k] * si_spher[j];
	}
      }

    // This is the metric for the spherical surface
    Tensor<2, data_t> Sigma;
    FOR2(i, j)
      {
        Sigma[i][j] = 0.0;
        FOR2(k, l)
	  {
	    Sigma[i][j] +=
	      Proj_spher[i][k] * Proj_spher[j][l] * spherical_g[k][l];
	  }
      }

    const data_t dArea = sqrt(Sigma[1][1] * Sigma[2][2] - Sigma[1][2] * Sigma[2][1]);

    return dArea;
  }

} // namespace CoordinateTransformations
#endif /* COORDINATETRANSFORMATIONS_HPP_ */
