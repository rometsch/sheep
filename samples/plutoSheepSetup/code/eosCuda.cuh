/* Define a custom eos for the locally isothermal disk or include the ideal eos. */
#ifndef EOS_CUDA_CUH
#if EOS == IDEAL
	#include "EOS/Ideal/eosCuda.cuh"
#elif EOS == ISOTHERMAL
	#define EOS_CUDA_CUH
	
	__forceinline__ __device__ real calcSoundSpeed(real v[], int pos, int dir, 
	                                               uint3 dataIdx)
	{
	    real x1 = cudaGrid.xgc1[dataIdx.x];
	    real x2 = cudaGrid.xgc2[dataIdx.y];
	
		/* Choose variables to lie on the interface if the loop
			is in its direction. */
		// in r direction
	    if (dir == IDIR)
	    {
	        x1 = (pos == FACE_CENTER ? cudaGrid.x1_if[dataIdx.x+1] 
	                                 : cudaGrid.xgc1[dataIdx.x]);
	    }
		// in theta direction
	    if (dir == JDIR)
	    {
	        x2 = (pos == FACE_CENTER ? cudaGrid.x2_if[dataIdx.y+1] 
	                                 : cudaGrid.xgc2[dataIdx.y]);
	    }
		// in phi direction there is no difference in the result since only R = r*sin(th) matters
	
	    return cudaIsoSoundSpeed * rsqrt(x1*sin(x2));
	}
#else
	printf("Selected EOS %s not found! Check definitions.h\n", EOS);
	QUIT_PLUTO(1);
#endif
#endif /* EOS_CUDA_CUH */
