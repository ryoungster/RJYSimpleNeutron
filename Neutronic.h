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
GSInitVR(real32** TriArray, real32* XArray, uint32_t N, real32** RI, uint32_t* NAr);

real32
GSStep(real32** TriArray, real32** PhiArrays, uint32_t N);

void
GSRun(real32** TriArray, real32** PhiArrays, uint32_t N, real32 GSEps);

real32*
CopyArray32(real32* Source, uint32_t Count)
{
    // Allocates New array to act as copy
    // I might not end up using this
    // All arrays are real32 for now (except Number of cells)
    real32* Dest = (real32*)malloc(Count*sizeof(real32));
    memcpy(Dest, Source, Count*sizeof(real32));
    return Dest;
}
inline void
MultArray32(real32* A, real32* B, real32* Dest, uint32_t Count)
{
    // All arrays are real32 for now (except Number of cells)
    // Does work if one of the arguments is the destination
    for (uint32_t i = 0; i < Count; i++)
    {
        Dest[i] = A[i] * B[i];
    }
}

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
GSFreeOut(real32** TriArray, real32** PhiArrays)
{
    // Phi[0] and X are considered output arrays
    // Will be responsibility of caller
    free(TriArray[0]);  
    free(TriArray[1]);
    free(TriArray[2]);
    free(TriArray[3]);
    free(PhiArrays[1]);
}

void
GSInitVR(real32** TriArray, real32* XArray, uint32_t N, real32** RInfo, uint32_t* NAr)
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
        real32 s  = RInfo[0][RIndex];
        real32 SA = RInfo[1][RIndex];
        real32 D  = RInfo[2][RIndex];
        real32 dX = RInfo[3][RIndex];
        TriArray[1][n]   += D/dX+SA*dX/2.f;
        TriArray[2][n]   = D/dX;
        TriArray[3][n]   += s*dX/2.f;
        TriArray[0][n+1] = D/dX;
        TriArray[1][n+1] += D/dX+SA*dX/2.f;
        TriArray[3][n+1] += s*dX/2.f;
        XArray[n+1]         = XArray[n] + dX;
    }
    //Right Vacuum boundary
    TriArray[1][N] += 0.5f;
    
    for (uint32_t n = 0; n < N+1; n++) 
    {
        TriArray[1][n]   = 1.0f/TriArray[1][n];
    }
}

// Double functions to eliminate memcpy
void
GSStepR(real32** TriArray, real32** PhiArrays, uint32_t N)
{
	
	
	PhiArrays[1][0]=(TriArray[3][0]+TriArray[2][0]*PhiArrays[0][1])*TriArray[1][0];
	
	for (uint32_t n = 1; n < N; n++) {
		PhiArrays[1][n]=(TriArray[3][n]+TriArray[2][n]*PhiArrays[0][n+1]+TriArray[0][n]*PhiArrays[1][n-1])*TriArray[1][n];
	}
	
	PhiArrays[1][N]=(TriArray[3][N]+TriArray[0][N]*PhiArrays[1][N-1])*TriArray[1][N];
}
real32
GSStepL(real32** TriArray, real32** PhiArrays, uint32_t N)
{
	real32 Convergence = 0.f;
	
	PhiArrays[0][0]=(TriArray[3][0]+TriArray[2][0]*PhiArrays[1][1])*TriArray[1][0];
	Convergence = (real32) fmax(Convergence, fabs(PhiArrays[0][0] - PhiArrays[1][0])/PhiArrays[0][0]);
	
	for (uint32_t n = 1; n < N; n++) {
		PhiArrays[0][n]=(TriArray[3][n]+TriArray[2][n]*PhiArrays[1][n+1]+TriArray[0][n]*PhiArrays[0][n-1])*TriArray[1][n];
		Convergence = (real32) fmax(Convergence, fabs(PhiArrays[0][n] - PhiArrays[1][n])/PhiArrays[0][n]);
	}
	
	PhiArrays[0][N]=(TriArray[3][N]+TriArray[0][N]*PhiArrays[0][N-1])*TriArray[1][N];
	Convergence = (real32) fmax(Convergence, fabs(PhiArrays[0][N] - PhiArrays[1][N])/PhiArrays[0][N]);
	
	return Convergence;
}

void
GSRun(real32** TriArray, real32** PhiArrays, uint32_t N, real32 Epsilon)
{
    real32 Convergence;
    do
    {
        GSStepR(TriArray, PhiArrays, N);
        Convergence = GSStepL(TriArray, PhiArrays, N);
    }
    while (Convergence > Epsilon);

}

struct RegionDesc
{
    real32* RInfo[4];
    uint32_t* NAr;
};
struct GSOutput
{
// Non standard hidden by w4201
    uint32_t N;
    union
    {
        real32* Arr[2];
        struct
        {
            real32* X;
            real32* Phi;
        };
    };
};
struct GSOutput
GaussSeidel(struct RegionDesc Regions, uint32_t RCount, real32 Epsilon)
{
    // Format
    // Source, SigmaA, D, Region Length
    struct GSOutput Out = {0};
    real32* TriArrays[4];
    real32* PhiArrays[2];
    real32* XArray;
    uint32_t N = 0;
    real32 Length = 0;

    
    Regions.RInfo[3] = CopyArray32(Regions.RInfo[3], RCount);

    for (uint32_t i = 0; i < RCount; i++)
    {
        N += Regions.NAr[i];
        Length += Regions.RInfo[3][i];
        Regions.RInfo[3][i] /= (real32)Regions.NAr[i];
    }
    GSAlloc(TriArrays, PhiArrays, &XArray, N);

    GSInitVR(TriArrays, XArray, N, Regions.RInfo, Regions.NAr);
    GSRun(TriArrays, PhiArrays, N, Epsilon);

    GSFreeOut(TriArrays, PhiArrays);
    free(Regions.RInfo[3]);
    Out.N = N;
    Out.X = XArray;
    Out.Phi = PhiArrays[0];
    return Out;
}


#endif
