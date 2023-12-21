#include <stdlib.h>
#include <math.h>                  
#include <stdint.h>

#ifndef NEUT_H
#define NEUT_H

#define internal static 
#define local_persist static 
#define global_variable static 
typedef int32_t bool32;
typedef float real32;
typedef double real64;
#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))
#define randReal32() (((real32)rand())/((real32) RAND_MAX))

#ifdef GRAPH_WIN32
// Win32 perforomance functions
// NOTE: Speed is inheritly compromised by lack of compiler flags for /O etc
#include <windows.h>
// #include <winnt.h> 
// #include <profileapi.h>
global_variable int64_t PerfCounterFreq;

inline void
Win32SetUpPerf(void)
{
    LARGE_INTEGER CounterFreq;
    QueryPerformanceFrequency(&CounterFreq);
    PerfCounterFreq = CounterFreq.QuadPart;
}
inline LARGE_INTEGER
Win32GetPerfCounter(void)
{
    LARGE_INTEGER Result;
    QueryPerformanceCounter(&Result);
    return (Result);
}

inline real32
Win32GetSecondsElapsed(LARGE_INTEGER Start, LARGE_INTEGER End)
{
    real32 Result = (real32) ((real64)(End.QuadPart - Start.QuadPart) / (real64)PerfCounterFreq);
    return Result;
}
#define SetUpTiming() Win32SetUpPerf()
#define GetTimerCount() Win32GetPerfCounter()
#define GetSecondsElapsed(s, e) Win32GetSecondsElapsed(s, e)
typedef LARGE_INTEGER TimerCount;
#endif


void
GSAlloc(real32** TriArray, real32** PhiArrays, real32** XArray, uint32_t N);

void
GSInitVR(real32** TriAnsArray, real32* XArray, uint32_t N, real32** RI, uint32_t* NAr);

real32
GSStep(real32** TriArray, real32** PhiArrays, uint32_t N);


void
GSAlloc(real32** TriArray, real32** PhiArrays, real32** XArray, uint32_t N)
{
    // Allocates (N + 1) * 7 arrays for Gauss Seidel
    // (DOCS) Initialization to all bits zero does not guarantee that a 
    // floating-point would be initialized to 0.0 (although that is true on all common platforms) (DOCS)
    // Here we operate assuming this to be the case
    TriArray[0]  = (real32*)calloc(N+1, sizeof(real32));
    TriArray[1]  = (real32*)calloc(N+1, sizeof(real32));
    TriArray[2]  = (real32*)calloc(N+1, sizeof(real32));
    TriArray[3]  = (real32*)calloc(N+1, sizeof(real32));
    PhiArrays[0] = (real32*)calloc(N+1, sizeof(real32));
    PhiArrays[1] = (real32*)calloc(N+1, sizeof(real32));
    *XArray      = (real32*)calloc(N+1, sizeof(real32));
		
}

void
GSInitVR(real32** TriAnsArray, real32* XArray, uint32_t N, real32** RI, uint32_t* NAr)
{   
    // Init vacuum right
    // All Arrays are assumed zero 
    uint32_t RIndex = 0;
    for (uint32_t n = 0, nInternal = 0; n < N; n++, nInternal++)
    {
        if (nInternal >= NAr[RIndex])
        {
            nInternal = 0;
            RIndex++;
        }
        real32 s  = RI[0][RIndex];
        real32 SA = RI[1][RIndex];
        real32 D  = RI[2][RIndex];
        real32 dX = RI[3][RIndex];
        TriAnsArray[1][n]   += D/dX+SA*dX/2.f;
        TriAnsArray[2][n]   = D/dX;
        TriAnsArray[3][n]   += s*dX/2.f;
        TriAnsArray[0][n+1] = D/dX;
        TriAnsArray[1][n+1] += D/dX+SA*dX/2.f;
        TriAnsArray[3][n+1] += s*dX/2.f;
        XArray[n+1]         = XArray[n] + dX;
    }
    //Right Vacuum boundary
    TriAnsArray[1][N] += 0.5f;
    
    for (uint32_t n = 0; n < N+1; n++) 
    {
        TriAnsArray[1][n]   = 1.0f/TriAnsArray[1][n];
    }
}

real32
GSStep(real32** TriArray, real32** PhiArrays, uint32_t N)
{
	real32 Convergence = 0.f;
	
	
	PhiArrays[1][0]=(TriArray[3][0]+TriArray[2][0]*PhiArrays[0][1])*TriArray[1][0];
	Convergence = (real32) fmax(Convergence, fabs(PhiArrays[1][0] - PhiArrays[0][0])/PhiArrays[1][0]);
	
	for (uint32_t n = 1; n < N; n++) {
		PhiArrays[1][n]=(TriArray[3][n]+TriArray[2][n]*PhiArrays[0][n+1]+TriArray[0][n]*PhiArrays[1][n-1])*TriArray[1][n];
		Convergence = (real32) fmax(Convergence, fabs(PhiArrays[1][n] - PhiArrays[0][n])/PhiArrays[1][n]);
	}
	
	PhiArrays[1][N]=(TriArray[3][N]+TriArray[0][N]*PhiArrays[1][N-1])*TriArray[1][N];
	Convergence = (real32) fmax(Convergence, fabs(PhiArrays[1][N] - PhiArrays[0][N])/PhiArrays[1][N]);
	
	memcpy(PhiArrays[0],PhiArrays[1],sizeof(real32)*(N+1));
	return Convergence;
}

#endif
