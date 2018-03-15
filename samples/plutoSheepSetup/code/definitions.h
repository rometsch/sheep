#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              POTENTIAL
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     13
#define  NBODY_SYS               NO

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               EXPLICIT
#define  ROTATING_FRAME          NO

/* -- nbody declarations -- */

#define CENTRAL_OBJECT        BINARY
#define CO_FEELS_DISK         NO
#define INDIRECT_TERM         NO
#define NO_OF_PLANETS         0

/* -- user-defined parameters (labels) -- */

#define	SigmaRef  			0
#define	Mplanet				1
#define	AspectRatio			2
#define	ViscosityAlpha		3
#define	Inclination			4
#define	Rplanet				5
#define	Smoothing			6
#define	ForceCutoff			7
#define DensityFloor		8
#define Pericenter          9
#define WDThetaBegRel		10
#define WDRIn				11
#define WDROut				12

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH             (g_inputParam[Rplanet]*CONST_au)
#define  UNIT_DENSITY            (CONST_Msun/(UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH))
#define  UNIT_VELOCITY           (sqrt(CONST_G*CONST_Msun/UNIT_LENGTH)/(2.*CONST_PI))

/* user defines for init functions, accessible in both init.c and initCuda.cu */

/* program choices */

/* planets */
#define	USE_PLANET_POTENTIAL	YES

/* Analysis [usrBegAna] */
#define USE_CUDA_ANALYSIS		YES
#define CALC_KINETIC_ENERGY		YES
#define  NUMBER_OF_ANALYSIS_VALUES   12
#define  ANA_VAL_HEADER   "M_DISK\tKE_R\tKE_TH\tKE_PHI\tRHO_MIN\tRHO_MAX\tJ_DISK_0\tJ_DISK_1\tJ_DISK_2\tF_0\tF_1\tF_2"
#define  ANA_VAL_OPERATORS   OP_SUM,OP_SUM,OP_SUM,OP_SUM,OP_MIN,OP_MAX,OP_SUM,OP_SUM,OP_SUM,OP_SUM,OP_SUM,OP_SUM
#define  AN_M_DISK                0
#define  AN_KE_R                  1
#define  AN_KE_TH                 2
#define  AN_KE_PHI                3
#define  AN_RHO_MIN               4
#define  AN_RHO_MAX               5
#define  AN_J_DISK                6
#define  AN_J_DISK_LEN            3
#define  AN_F                     9
#define  AN_F_LEN                 3
/* Analysis [usrEndAna] */

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       YES
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             VANLEER_LIM
