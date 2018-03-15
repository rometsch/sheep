/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Disk-Planet interaction problem for planets on cirluar inclined orbits.


  -# g_inputParam[Mstar]: controls the star mass (in solar masses)
  -# g_inputParam[Mplanet]: sets the planet mass (in solar masses)
  -# g_inputParam[ViscosityAlpha]: sets alpha value for the alpha viscosity model
  -# g_inputParam[AspectRatio]: sets the aspect ratio h of the disk
  -# g_inputParam[Inclination]: sets the inclination of the planet


   Computaion is carried out in the corotating frame of the planet.

  \author Thomas Rometsch
  \date   Jul 23, 2017

  \b References:
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double inc, h, b, mu, mp, rho_c, Omega_P, r_hill, r_sm, r_sm_ip4, r_sm_ip3, r_sm_ip1, r_planet, r_star;

void planet_position (double *v, double inc);
void SaveAnalysisValues(double* analysisValues);

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
	static int first_run = 1;

	double r, th, R, R_m32, z, H, OmegaK, cs;
	double sin_th;

	/* ----------------------------------------------------------------------
	  	Calculate constant quantities only once.
	   ---------------------------------------------------------------------- */

	if (first_run == 1) {
		double Rmin, Rmax;
		Rmin = g_domBeg[IDIR];
		Rmax = g_domEnd[IDIR];
		/* disk properties */
		h = g_inputParam[AspectRatio];
		rho_c = g_inputParam[SigmaRef]/(sqrt(2.0*CONST_PI)*h);
		/* planet orbit */
		mu = g_inputParam[Mplanet];
		mp = g_inputParam[Mplanet];
		b = g_inputParam[ForceCutoff];
		inc = g_inputParam[Inclination]/180*CONST_PI;
		Omega_P = 2.0*CONST_PI*sqrt(1.0+mu);
		r_planet = 1.0/(1.0+mu);
		r_star = mu/(1.0+mu);
		r_hill = pow(mu/3.0, 1.0/3.0);
		/* smoothing length and inverse powers of it */
		r_sm = g_inputParam[Smoothing]*r_hill;
		r_sm_ip1 = 1.0/r_sm;
		r_sm_ip3 = 1.0/(r_sm*r_sm*r_sm);
		r_sm_ip4 = 1.0/(r_sm*r_sm*r_sm*r_sm);
		#if EOS == IDEAL
			g_gamma = g_inputParam[Gamma];
		#elif EOS == ISOTHERMAL
			g_isoSoundSpeed = 2.0*CONST_PI*h;
		#endif
		/* density floor */
		g_smallDensity = g_inputParam[DensityFloor];
		/* finished, set flag */
		first_run = 0;
	}


	r     = x1;
	th    = x2;
	R     = r*sin(th);
	z     = r*cos(th);
	R_m32 = 1.0/(R*sqrt(R));

	H      = h*R;
	OmegaK = 2.0*CONST_PI*R_m32;
	cs     = H*OmegaK;

	sin_th = sin(th);

	us[VX1] = us[VX2] = 0.0;

	#if EOS == ISOTHERMAL

	/* Equilibrium density distribution for a locally isothermal dis.
	   with p = -1.5, q = -1 (see Nelson2013+). */
	us[RHO] = rho_c*R_m32*exp((sin_th-1)/(h*h));
	 /* Equilibrium angular velocity for a locally isothermal disk
	    with p = -1.5, q = -1 (see Nelson2013+). */
	us[VX3] = R*OmegaK*sqrt(-2.5*h*h+sin_th);
	g_isoSoundSpeed = h*2.0*CONST_PI;
	#else
		printf("Selected eos not supported in this setup!\n");
		QUIT_PLUTO(1);
	#endif

	/* Density Floor */
	if ( us[RHO] < g_smallDensity )
	{
		us[RHO] = g_smallDensity;
	}

}


/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*
 *  Compute forces as decribed in Bitsch & Kley (2011)
 *
 *********************************************************************** */
{
	// Array for output values for use with MPI
	double analysisValues[NUMBER_OF_ANALYSIS_VALUES] = {0.0};
	Operators reductionOperators[NUMBER_OF_ANALYSIS_VALUES] = {ANA_VAL_OPERATORS};
	for (int n=0; n < NUMBER_OF_ANALYSIS_VALUES; n++)
	{
		if (reductionOperators[n] == OP_PROD)
		{
			analysisValues[n] = 1.0;
		}
	}

	double t;
	double r, th, phi;
	double vr, vth, vphi;
	double x_p, y_p, z_p, xi;
	double s, s_x, s_y, s_z;
	double sin_th, cos_th, sin_phi, cos_phi;
	double dV, dm, rho;
	double jx, jy, jz;
	double force, cutoff;
	int i,j,k;

	t = g_time;

	/* Choose to calculate analysis on cpu or gpu */
	#if USE_CUDA_ANALYSIS == YES
		callAnalysisReductionKernel(analysisValues, reductionOperators);
	#else

		double v[3] = {0.0};
		planet_position(v, inc);
		x_p = r_planet*v[0];
		y_p = r_planet*v[1];
		z_p = r_planet*v[2];

		// Iterate over domain to calculate force acting on planet.
	KDOM_LOOP(k) {
		phi = grid[KDIR].xgc[k];
		vphi = d->Vc[VX3][k][j][i];
		JDOM_LOOP(j) {
			th = grid[JDIR].xgc[j];
			vth = d->Vc[VX2][k][j][i];
			IDOM_LOOP(i) {
			r = grid[IDIR].xgc[i];
			vr = d->Vc[VX1][k][j][i];

			/* Get cell colume and calculate cell mass. */
			rho = d->Vc[RHO][k][j][i];
			dV = grid[IDIR].dV[i] * grid[JDIR].dV[j] * grid[KDIR].dV[k];
			dm =  rho * dV;

			/* sins and coss */
			sin_th = sin(th);
			cos_th = cos(th);
			sin_phi = sin(phi);
			cos_phi = cos(phi);

    		b = g_inputParam[ForceCutoff];
			/*************** end of cpu specific code ****************/


			/* Calculate vector and distance from cell center to planet. */
		    s_x = x_p - r*cos_phi*sin_th;
		    s_y = y_p - r*sin_phi*sin_th;
		    s_z = z_p - r*cos_th;
		    s = sqrt(s_x*s_x + s_y*s_y + s_z*s_z);


		    /* Fermi type cutoff as in Kley et. al. (2009).*/
		    cutoff = 1.0/ (1.0 + exp( -10.0*(s/r_hill - b)/b ) );

			/* Calculate force and store it to analysis value array. */
		    force = 4.0*CONST_PI*CONST_PI/(s*s*s)*cutoff;

		    analysisValues[AN_F] += - dm*mp*s_x*force;
		    analysisValues[AN_F + 1] += - dm*mp*s_y*force;
		    analysisValues[AN_F + 2] += - dm*mp*s_z*force;

			/* Calculate total mass of disk.
		       Write to value array, remember to just set the variable insted of +=.
		       the code takes care of the reduction. */
		    analysisValues[AN_M_DISK] += dm;

			/* Min and max density */
			analysisValues[AN_RHO_MIN] = MIN(rho, analysisValues[AN_RHO_MIN]);
			analysisValues[AN_RHO_MAX] = MAX(rho, analysisValues[AN_RHO_MAX]);

			/* Calculate kinetic energy of cell and the contribution from each direction. */
			analysisValues[AN_KE_R] += 0.5*dm*vr*vr;
			analysisValues[AN_KE_TH] += 0.5*dm*vth*vth;
			analysisValues[AN_KE_PHI] += 0.5*dm*vphi*vphi;


			/* Calculate the contribution of each cell to the angular momentum of the disk. */
			jx = dm*r*( -vth*sin_phi - vphi*cos_th*cos_phi );
			jy = dm*r*( vth*cos_phi  - vphi*cos_th*sin_phi );
			jz = dm*r*( vphi*sin_th );

			analysisValues[AN_J_DISK] += jx;
			analysisValues[AN_J_DISK+1] += jy;
			analysisValues[AN_J_DISK+2] += jz;

			}
		}
	}
	#endif

	/* Reduce arrays for parallel computing. */
	#ifdef PARALLEL
	double tmp;
	for (int n = 0; n < NUMBER_OF_ANALYSIS_VALUES; n++)
	{
		switch(reductionOperators[n]) {
			case OP_SUM: MPI_Allreduce(analysisValues+n, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				break;
			case OP_PROD: MPI_Allreduce(analysisValues+n, &tmp, 1, MPI_DOUBLE, MPI_PROD, MPI_COMM_WORLD);
				break;
			case OP_MIN: MPI_Allreduce(analysisValues+n, &tmp, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
				break;
			case OP_MAX: MPI_Allreduce(analysisValues+n, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				break;
		}
		analysisValues[n] = tmp;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	#endif

	/* Write ascii file "analysis_values.dat" to disk */
	if (prank == 0)
	{
		SaveAnalysisValues(analysisValues);
	}
}

/* ********************************************************************* */
void SaveAnalysisValues(double* analysisValues) {
/*
 *  Save analysis values to a file.
 *
 *********************************************************************** */
	char fname[512];
	static double tpos = -1.0;
	FILE *fp;

	sprintf(fname, "%s/analysis_values.dat", RuntimeGet()->output_dir);

	/* Write header */
	if (g_stepNumber == 0)
	{
		fp = fopen(fname, "w");
		fprintf(fp, "# Analysis values of pluto Simulation\n");
		/* Print the variable string defined in the definition.h file. */
		fprintf(fp, "#vars:time\tniter\t%s\n", ANA_VAL_HEADER);
	}
	else
	{
		if (tpos < 0.0)
		{
			char sline[1024];
			fp = fopen(fname, "r");
			while (fgets(sline, 1024, fp));
			sscanf(sline, "%lf*", &tpos);
			fclose(fp);
		}
		fp = fopen(fname, "a");
	}

	/* Write data. */
	if (g_time > tpos)
	{
		fprintf(fp, "%12.6e\t%12ld" , g_time ,g_stepNumber);
		for (int n=0; n<NUMBER_OF_ANALYSIS_VALUES; n++) {
			fprintf(fp, "\t%12.12e", analysisValues[n]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
}


#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
	double s_p, s_star, R, r, z, th, x, y, phiplanet, rsm;
	double x_p, y_p, z_p, x_s, y_s, z_s, t, phi;

	/* get planet and star position */
	double v[3] = {0.0};
	planet_position(v, inc);
	x_p = r_planet*v[0];
	y_p = r_planet*v[1];
	z_p = r_planet*v[2];
	x_s = -r_star*v[0];
	y_s = -r_star*v[1];
	z_s = -r_star*v[2];

	/* calculate position of cell center */
	x = x1*sin(x2)*cos(x3);
	y = x1*sin(x2)*sin(x3);
	z = x1*cos(x2);

	/* calculate distance of cell center to star and planet */
	s_star =   sqrt( (x-x_s)*(x-x_s) + (y-y_s)*(y-y_s) + (z-z_s)*(z-z_s) );
	s_p =      sqrt( (x-x_p)*(x-x_p) + (y-y_p)*(y-y_p) + (z-z_p)*(z-z_p) );

	#if USE_PLANET_POTENTIAL == YES
		/* calculate planet potential
		   smoothing length r_sm and planet-star mass ratio mu are constant and a global variables */
		if ( s_p  > r_sm ) phiplanet = mu/s_p;
		else phiplanet = mu*(s_p*s_p*s_p*r_sm_ip4 - 2.*s_p*s_p*r_sm_ip3 + 2.*r_sm_ip1);
	#else
		phiplanet = 0.0;
	#endif

	/* calculate gravitational potential of star and planet */
	phi  = - 4.0*CONST_PI*CONST_PI*(1.0/s_star + phiplanet);

	return phi;
}
#endif


/* ************************************************************************ */
void planet_position (double *v, double inc)
/*
 *  Calculate the coordinates of a planet on a circular inclined orbit.
 *
 * \param [in] 	t	time
 * \param [in,out]	double array of size 3 with position of the planet
 *
 *************************************************************************** */
{
    double t, Omega_P, xi, cos_xi, sin_xi, cos_i, sin_i, mp, ms, mu, r;

	mp = g_inputParam[Mplanet];
	mu = mp;

	/* Planet position is fixed to 1 (choice of units) */
	r = 1.0;

    Omega_P = 2.0*CONST_PI*sqrt(1.0+mu);

    /* Time of actual integration step */
    if (g_intStage == 1)
        t = g_time;
    else if (g_intStage == 2)
        t = g_time + g_dt;
    else if (g_intStage == 3)
        t = g_time + 0.5*g_dt;


	double omega_peri = 0.0;
	xi = Omega_P*t + omega_peri;

	cos_i = cos(inc);
	sin_i = sin(inc);
	cos_xi = cos(xi);
	sin_xi = sin(xi);

	v[0] = 1.0/(1.0+mu)*r*cos_xi;
	v[1] = 1.0/(1.0+mu)*r*cos_i*sin_xi;
	v[2] = 1.0/(1.0+mu)*r*sin_i*sin_xi;
}
