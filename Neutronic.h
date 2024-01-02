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

inline real32*
CopyArray32(real32* Source, uint32_t Count)
{
    // Allocates New array to act as copy
    real32* Dest = (real32*)malloc(Count*sizeof(real32));
    memcpy(Dest, Source, Count*sizeof(real32));
    return Dest;
}
inline void
MultArray32(real32* A, real32* B, real32* Dest, uint32_t Count)
{
    // Does work if one of the arguments is the destination (and baring order issues)
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
    real32 s  = RInfo[0][RIndex];
    real32 SA = RInfo[1][RIndex];
    real32 D  = RInfo[2][RIndex];
    real32 dX = RInfo[3][RIndex];
    for (uint32_t n = 0, nInternal = 0; n < N; n++, nInternal++)
    {
        if (nInternal >= NAr[RIndex])
        {
            nInternal = 0;
            RIndex++;
            s  = RInfo[0][RIndex];
            SA = RInfo[1][RIndex];
            D  = RInfo[2][RIndex];
            dX = RInfo[3][RIndex];
        }
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
    uint32_t RCount;
    real32* RInfo[4];
    // NAr is arguably a property of the numerical method, 
    // but it is per region and determines size of output
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
GaussSeidel(struct RegionDesc Regions, real32 Epsilon)
{
    // Format
    // Source, SigmaA, D, Region Length
    // Left Reflecting, Right Vacuum
    struct GSOutput Out = {0};
    real32* TriArrays[4];
    real32* PhiArrays[2];
    real32* XArray;
    uint32_t N = 0;
    real32 Length = 0;

    // Copy to avoid modifying parameter
    Regions.RInfo[3] = CopyArray32(Regions.RInfo[3], Regions.RCount);

    for (uint32_t i = 0; i < Regions.RCount; i++)
    {
        N += Regions.NAr[i];
        Length += Regions.RInfo[3][i];
        Regions.RInfo[3][i] /= (real32)Regions.NAr[i];
    }

    GSAlloc(TriArrays, PhiArrays, &XArray, N);

    GSInitVR(TriArrays, XArray, N, Regions.RInfo, Regions.NAr);

    GSRun(TriArrays, PhiArrays, N, Epsilon);

    GSFreeOut(TriArrays, PhiArrays);
    free(Regions.RInfo[3]); // Free Allocated Edited Array

    Out.N = N;
    Out.X = XArray;
    Out.Phi = PhiArrays[0];
    return Out;
}

inline real32
IsoMu()
{
    return (randReal32()*2.f-1.f)
}

struct MCOutput
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
struct MCOutput
MonteCarlo(struct RegionDesc Regions, uint32_t Histories)
{
    // MonteCarlo using collision estimator
    // Format Source, SS/ST, SigmaT, Length
    struct MCOutput Out = {0};

    uint32_t N = 0;
    real32 Length = 0;

    real32 SourceSpace = 0;
    real32* WeightedSAr = malloc(Regions.RCount * sizeof(real32));
    real32* DXAr = malloc(Regions.RCount * sizeof(real32));
    // TODO: Look at predividing frequently used inverses

    for (uint32_t i = 0; i < Regions.RCount; i++)
    {
        N += Regions.NAr[i];
        Length += Regions.RInfo[3][i];
        WeightedSAr[i] =  Regions.RInfo[0][i]*Regions.RInfo[3][i];
        SourceSpace += WeightedSAr[i];
        DXAr[i] = Regions.RInfo[3][i] / (real32)Regions.NAr[i];
    }
    uint32_t* BinArray = calloc(N, sizeof(uint32_t));
    real32* XArray = calloc(N, sizeof(real32));
    real32* PhiArray = malloc(N * sizeof(real32));

    { //TODO:Pull out midpoint xgen into a function
        uint32_t RIndex = 0 , nInternal = 0;
        real32 dX = DXAr[RIndex];
        for (uint32_t n = 0; n < N-1; n++, nInternal++)
        {
            if (nInternal >= Regions.NAr[RIndex])
            {
                nInternal = 0;
                RIndex++;
                dX = DXAr[RIndex];
            }
            XArray[n] += dX/2;
            XArray[n+1] += XArray[n] + dX/2;
        }
        if (nInternal >= Regions.NAr[RIndex])
        {
            nInternal = 0;
            RIndex++;
            dX = DXAr[RIndex];
        }
        XArray[N-1] += dX/2;
    }
    
    for (uint32_t n = 0; n < Histories; n++)
    {
        real32 x = randReal32() * SourceSpace;
        real32 mu = IsoMu();
        uint32_t RID;
        // Find spawn location using weighted locations
        for (RID = 0; RID < Regions.RCount; RID++)
        {
            if (x < WeightedSAr[RID])
            {
                // Divide by Source
                //TODO:Predivide this?
                x = x / Regions.RInfo[0][RID];
                break;
            }
            else
            {
                x -= WeightedSAr[RID];
            }
        }
        // Distance to collision in region
        //TODO:Predivide this? Wrap it in function?
        real32 DColl = -1.f*log(randReal32()) / Regions.RInfo[2][RID]; 
        x += mu * DColl;
        if (;;) // Infinite loop until neutron death
        {
            if (x<0)
            {
                // Left Overflow
                RID--;
                if (RID < 0)
                {
                    //Reflected Boundary
                    mu = -mu;
                    RID++;
                    x = -x;
                }
                else
                {
                    // Enters new region
                    // Translate (Negative) Distance to new SigT
                    x *= Regions.RInfo[2][RID+1]/ Regions.RInfo[2][RID];
                    // Offset (Negative) Distance by Region Length
                    x += Regions.RInfo[3][RID];
                }

            }
            else if (x > Regions.RInfo[3][RID])
            {
                // Right Overflow
                RID++;
                if (RID >= Regions.RCount)
                {
                    // Vacuum Boundary
                    break
                }
                else
                {
                    // Enters new region
                    // Offset Distance by Old Region Length
                    x -= Regions.RInfo[3][RID-1];
                    // Translate  Distance to new SigT
                    x *= Regions.RInfo[2][RID-1]/ Regions.RInfo[2][RID];
                }
            }
            else
            {
                // Collides in region
                uint32_t BinID = 0;
                for (uint32_t i = 0; i < RID; i++)
                {
                    // Previous regions
                    BinID += Regions.NAr[i];
                }
                BinID += floor((x * (real32) Regions.NAr[RID]) / Regions.RInfo[3][RID]);
                BinArray[BinID]++;

                if (randReal32() > Regions.RInfo[1])
                {
                    // Absorbed
                    break;
                }
                else
                {
                    // Scattered
                    mu = IsoMu();
                    DColl = -1.f*log(randReal32()) / Regions.RInfo[2][RID]; 
                    x += mu * DColl;
                }
            }
        }
    }

    {
    // PhiArray from Bins
    uint32_t RIndex = 0;
    real32 InvST = 1/RInfo[2][RIndex];
    real32 InvDX = 1/DXAr[RIndex];
    real32 HistScale = SourceSpace / (real32) Histories;
    for (uint32_t n = 0, nInternal = 0; n < N; n++, nInternal++)
    {
        if (nInternal >= NAr[RIndex])
        {
            nInternal = 0;
            RIndex++;
            InvST = 1/RInfo[2][RIndex];
            InvDX = 1/DXAr[RIndex];
        }
        PhiArray[n] = BinArray[n] * InvST * InvDX * HistScale;
    }
    }

    free(BinArray);
    Out.N = N;
    Out.X = XArray;
    Out.Phi = PhiArray;
    return Out;
}


#endif
