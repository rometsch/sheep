extern "C" {

#include "pluto.h"
#include "plutoCuda.cuh"

#define PI 3.141592654

/* ************************************************************************** */
__device__ void cudaUserDefBoundary (int side, uint3 dataIdx, dim3 dataDim)
/*!
 *	Assign user-defined boundary conditions.
 *
 * \param [in] side	  specifies the boundary side where ghost zones need
 *			  to be filled. It can assume the following
 *			  pre-definite values: X1_BEG, X1_END,
 *					   X2_BEG, X2_END,
 *					   X3_BEG, X3_END.
 * \param [in] dataIdx
 * \param [in] dataDim
 *
 * The primitive variables of the ghost cells can be set the following way:
 *	cudaV.rho[gpu1D(dataIdx, dataDim)] = ...
 *	cudaV.vx1[gpu1D(dataIdx, dataDim)] = ...
 *	cudaV.vx2[gpu1D(dataIdx, dataDim)] = ... (only if COMPONENTS > 1)
 *	cudaV.vx3[gpu1D(dataIdx, dataDim)] = ... (only if COMPONENTS > 2)
 *	cudaV.prs[gpu1D(dataIdx, dataDim)] = ...
 *
 * The coordinates of the cell can be accessed the following way:
 *	x1 coordinate: cudaGrid.x1[dataIdx.x];
 *	x2 coordinate: cudaGrid.x2[dataIdx.y];
 *	x3 coordinate: cudaGrid.x3[dataIdx.z];
 *
 **************************************************************************** */
{
	int idx = gpu1D(dataIdx, dataDim);

	/* Density floor */
	if (cudaV.rho[idx] < cudaSmallDensity)
    {
		cudaV.rho[idx] = cudaSmallDensity;
	}

	/* Wave damping in radial direction */
	if (side == 0)
	{
		real r = cudaGrid.xgc1[dataIdx.x];
		real th = cudaGrid.xgc2[dataIdx.y];

		real beta = abs(th - 0.5*CONST_PI);

		real r_wd_in = cudaInputParam[WDRIn];
		real r_wd_out = cudaInputParam[WDROut];

		int bndIdx = 0;
		// access r_min
		real r_min = cudaGrid.xgc1[cudaIBEG];
		// access r_max
		real r_max = cudaGrid.xgc1[cudaIEND];
		// access theta_min
		real beta_max = 0.5*CONST_PI - cudaGrid.xgc2[cudaJBEG];

		real delta = cudaInputParam[WDThetaBegRel]*beta_max;

		if ( r <= r_wd_in || r >= r_wd_out || beta >= delta ) {
			/* Wave damping */


			real tau_r = 1.0;
			real ramp_r = 1.0;

			real tau_th = 1.0;
			real ramp_th = 1.0;

			real tau = 1.0;
			real ramp = 1.0;

			/* Select the highest damping factor.
			   This avoids using an amplitude higher than 1. */
			if ( r <= r_wd_in ) {
				/* inner wave killing zone */
				tau_r = sqrt(r_min*r_min*r_min);
				ramp_r = (r - r_wd_in)*(r - r_wd_in) / ( (r_wd_in - r_min)*(r_wd_in - r_min) );
			}

			if ( r >= r_wd_out ) {
				/* outer wave killing zone */
				tau_r = sqrt(r_max*r_max*r_max);
				ramp_r = (r - r_wd_out)*(r - r_wd_out) / ( (r_max - r_wd_out)*(r_max - r_wd_out) );
			}

			if ( beta >= delta ) {
				tau_th = sqrt(r*r*r);
				ramp_th = (beta-delta)*(beta-delta) / ( (beta_max - delta)*(beta_max -delta) );
			}

			if ( ramp_th >= ramp_r ) {
				ramp = ramp_th;
				tau = tau_th;
			} else {
				ramp = ramp_r;
				tau = tau_r;
			}

			/* inside a wave killing zone */
			real lambda = cuda_dt/tau*ramp;
			/* Calculate equilibrium values */
			real h = cudaInputParam[AspectRatio];
			real rho_c = cudaInputParam[SigmaRef]/(sqrt(2.0*PI)*h);
			real sin_th = sin(th);
			real R = r*sin_th;
			real H = R*h;
			real z = r*cos(th);
			real R_m32 = 1.0/(R*sqrt(R));
			real OmegaK = 2.0*PI*R_m32;

			/* Equilibrium values as in init.c */
			real rho = rho_c*R_m32*exp((sin_th-1)/(h*h));
			real v_phi = R*OmegaK*sqrt(-2.5*h*h+sin_th);

			cudaV.vc[RHO][idx] += - (cudaV.vc[RHO][idx] - rho)*lambda;
			cudaV.vc[VX3][idx] += - (cudaV.vc[VX3][idx] - v_phi)*lambda;
			/* Damp radial and vertical velocity to 0. */
			cudaV.vc[VX2][idx] += - (cudaV.vc[VX2][idx])*lambda;
			cudaV.vc[VX1][idx] += - (cudaV.vc[VX1][idx])*lambda;

			if (cudaV.rho[idx] < cudaSmallDensity)
			{
				cudaV.rho[idx] = cudaSmallDensity;
			}

		}
	}

}

#if BODY_FORCE != NO
/* ********************************************************************* */
__device__ real cudaBodyForceVector(int dir, uint3 dataIdx, dim3 dataDim)
/*!
 * Returns the component of the acceleration vector in a given direction.
 *
 * \param [in] dir	specifies the component of the acceleration vector.
 * \param [in] dataIdx
 * \param [in] dataDim
 *
 * The coordinates of the cell can be accessed the following way:
 *	x1 coordinate: cudaGrid.x1[dataIdx.x];
 *	x2 coordinate: cudaGrid.x2[dataIdx.y];
 *	x3 coordinate: cudaGrid.x3[dataIdx.z];
 *
 * The cell-centered primitive variables of the cell can be accessed
 * the following way:
 *	rho: cudaV.rho[gpu1D(dataIdx, dataDim)]
 *	vx1: cudaV.vx1[gpu1D(dataIdx, dataDim)]
 *	vx2: cudaV.vx2[gpu1D(dataIdx, dataDim)] (only if COMPONENTS > 1)
 *	vx3: cudaV.vx3[gpu1D(dataIdx, dataDim)] (only if COMPONENTS > 2)
 *	prs: cudaV.prs[gpu1D(dataIdx, dataDim)]
 *
 *********************************************************************** */
{
	return 0.0;
}

/* ********************************************************************* */
__device__ real cudaBodyForcePotential(real x1, real x2, real x3)
/*!
 * Returns the gravitational potential as a function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
	real x, y, z;
	real x_p, y_p, z_p;
	real x_s, y_s, z_s;
	real r_hill, r_sm, r_sm_ip4, r_sm_ip3, r_sm_ip1;
	real mp ,mu, phi, phiplanet;
	real s_star, s_p;

	mu = cudaInputParam[Mplanet];

	/* get planet and star position */
	x_p = cpuCalcParams[P_POS_X];
	y_p = cpuCalcParams[P_POS_Y];
	z_p = cpuCalcParams[P_POS_Z];
	x_s = cpuCalcParams[S_POS_X];
	y_s = cpuCalcParams[S_POS_Y];
	z_s = cpuCalcParams[S_POS_Z];

	/* calculate position of cell center */
	x = x1*sin(x2)*cos(x3);
	y = x1*sin(x2)*sin(x3);
	z = x1*cos(x2);

	/* calculate distance of cell center to star and to planet */
	s_star = sqrt( (x_s-x)*(x_s-x) + (y_s-y)*(y_s-y) + (z_s-z)*(z_s-z) );
	s_p =	 sqrt( (x_p-x)*(x_p-x) + (y_p-y)*(y_p-y) + (z_p-z)*(z_p-z) );

	#if USE_PLANET_POTENTIAL == YES
		/* Smoothing length = parameter*R_hill */
		r_hill = pow(mu/3.0, 1.0/3.0);
		r_sm = r_hill*cudaInputParam[Smoothing];
		/* inverse powers of the smoothing length */
		r_sm_ip1 = 1.0/r_sm;
		r_sm_ip3 = 1.0/(r_sm*r_sm*r_sm);
		r_sm_ip4 = 1.0/(r_sm*r_sm*r_sm*r_sm);
		/* calculate planet potential */
		if ( s_p  > r_sm ) phiplanet = mu/s_p;
		else phiplanet = mu*(s_p*s_p*s_p*r_sm_ip4 - 2.*s_p*s_p*r_sm_ip3 + 2.*r_sm_ip1);
	#else
		phiplanet = 0.0;
	#endif

	/* calculate gravitational potential of star and planet */
	phi	 = - 4.0*PI*PI*(1.0/s_star + phiplanet);

	return phi;
}
#endif

#if VISCOSITY != NO
__device__ void cudaViscNu(real *v, double x1, double x2, double x3,
					double *nu1, double *nu2)
{
	double Omega_K, R, H,alpha;
	R = x1*sin(x2);
	H = cudaInputParam[AspectRatio]*R;
	alpha = cudaInputParam[ViscosityAlpha];
	Omega_K = 2.0*PI/(R*sqrt(R));
	/* Calculate alpha viscosity and multiply by rho b/c pluto wants it. */
	*nu1 = alpha*H*H*Omega_K*v[RHO];
	*nu2 = 0.0;
}
#endif

#if USE_CUDA_REDUCTION == YES
__device__ void calcReductionValues(real *values, uint3 dataIdx, dim3 dataDim)
{

}
#endif

#if USE_CUDA_ANALYSIS == YES
__device__ void calcAnalysisValues(real *analysisValues, uint3 dataIdx, dim3 dataDim)
{
	int idx = gpu1D(dataIdx, dataDim);

	real x_p = cpuCalcParams[P_POS_X];
	real y_p = cpuCalcParams[P_POS_Y];
	real z_p = cpuCalcParams[P_POS_Z];

	real r, phi, th;
	real sin_th, cos_th, sin_phi, cos_phi;
	real vr, vth, vphi;
	real dV, dm, rho;
	real r_hill;
	real s, s_x, s_y, s_z;
	real force, cutoff;
	real mu, mp, b;
	real jx, jy, jz;

	/* get planet mass and calculate hill radius */
	mu = cudaInputParam[Mplanet];
	mp = cudaInputParam[Mplanet];
	r_hill = pow(mu/3.0, 1.0/3.0);


	/* get position of cell in spherical coordinates */
	r = cudaGrid.xgc1[dataIdx.x];
	th = cudaGrid.xgc2[dataIdx.y];
	phi = cudaGrid.xgc3[dataIdx.z];
	/* Get sines and cosines of th and phi */
	sin_th = sin(th);
	cos_th = cos(th);
	sin_phi = sin(phi);
	cos_phi = cos(phi);

	/* get velocities of cells in spherical coordinates */
	vr = cudaV.vc[VX1][idx];
	vth = cudaV.vc[VX2][idx];
	vphi = cudaV.vc[VX3][idx];

	/* Get cell colume and calculate cell mass. */
	rho = cudaV.rho[idx];
	dV =  cudaGrid.dVx1[dataIdx.x]
	 *cudaGrid.dVx2[dataIdx.y]
	 *cudaGrid.dVx3[dataIdx.z];
	dm = rho * dV;

	b = cudaInputParam[ForceCutoff];
	/*************** end of gpu specific code ****************/

	/* Calculate vector and distance from cell center to planet. */
	s_x = x_p - r*cos_phi*sin_th;
	s_y = y_p - r*sin_phi*sin_th;
	s_z = z_p - r*cos_th;
	s = sqrt(s_x*s_x + s_y*s_y + s_z*s_z);


	/* Fermi type cutoff as in Kley et. al. (2009).*/
	cutoff = 1.0/ (1.0 + exp( -10.0*(s/r_hill - b)/b ) );

	/* Calculate force and store it to analysis value array. */
	force = 4.0*PI*PI/(s*s*s)*cutoff;

	analysisValues[AN_F] = - dm*mp*s_x*force;
	analysisValues[AN_F + 1] = - dm*mp*s_y*force;
	analysisValues[AN_F + 2] = - dm*mp*s_z*force;

	/* Calculate total mass of disk.
	   Write to value array, remember to just set the variable insted of +=.
	   the code takes care of the reduction. */
	analysisValues[AN_M_DISK] = dm;

	/* Min and max density */
	analysisValues[AN_RHO_MIN] = rho;
	analysisValues[AN_RHO_MAX] = rho;

	/* Calculate kinetic energy of cell and the contribution from each direction. */
	analysisValues[AN_KE_R] = 0.5*dm*vr*vr;
	analysisValues[AN_KE_TH] = 0.5*dm*vth*vth;
	analysisValues[AN_KE_PHI] = 0.5*dm*vphi*vphi;


	/* Calculate the contribution of each cell to the angular momentum of the disk. */
	jx = dm*r*( -vth*sin_phi - vphi*cos_th*cos_phi );
	jy = dm*r*( vth*cos_phi	 - vphi*cos_th*sin_phi );
	jz = dm*r*( vphi*sin_th );

	analysisValues[AN_J_DISK] = jx;
	analysisValues[AN_J_DISK+1] = jy;
	analysisValues[AN_J_DISK+2] = jz;
}

#endif

#if USE_CPU_CALC_PARAMETERS == YES
void calcCpuParams(real *values)
{
	double t, Omega_P, peri, xi, inc, cos_xi, sin_xi, cos_i, sin_i, mu, r;

	inc = g_inputParam[Inclination]/180*PI;
	mu = g_inputParam[Mplanet];
	peri = g_inputParam[Pericenter];


	/* Planet position is fixed to 1 (choice of units) */
	r = 1.0;

	Omega_P = 2.0*CONST_PI;

	/* Time of actual integration step */
	if (g_intStage == 1)
	t = g_time;
	else if (g_intStage == 2)
	t = g_time + g_dt;
	else if (g_intStage == 3)
	t = g_time + 0.5*g_dt;

	xi = Omega_P*t + peri;

	cos_i = cos(inc);
	sin_i = sin(inc);
	cos_xi = cos(xi);
	sin_xi = sin(xi);

	values[P_POS_X] = 1.0/(1.0+mu)*r*cos_xi;
	values[P_POS_Y] = 1.0/(1.0+mu)*r*cos_i*sin_xi;
	values[P_POS_Z] = 1.0/(1.0+mu)*r*sin_i*sin_xi;
	values[S_POS_X] = -mu/(1.0+mu)*r*cos_xi;
	values[S_POS_Y] = -mu/(1.0+mu)*r*cos_i*sin_xi;
	values[S_POS_Z] = -mu/(1.0+mu)*r*sin_i*sin_xi;

}
#endif

} /* extern "C" */
