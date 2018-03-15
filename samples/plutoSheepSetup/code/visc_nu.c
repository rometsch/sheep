/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{ 
  double OmegaK, R, H, alpha;
  R = x1*sin(x2);
  H = g_inputParam[AspectRatio]*R;
  OmegaK = 2.0*CONST_PI/(R*sqrt(R));
  alpha = g_inputParam[ViscosityAlpha];
  *nu1 = alpha*H*H*OmegaK*v[RHO];
  *nu2 = 0.0;
}
