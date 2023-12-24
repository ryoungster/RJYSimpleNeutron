#include <stdio.h>
#include <stdlib.h>
#include <math.h>                  
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "dislin.h"
#include "Neutronic.h"

//real32 GSErrorArray[GSLEN];
//real32 TTimeArray[TLEN];

#define MaxRows 20000
real32 GraphArray[2][MaxRows];



// TODO: Centralize naming
// TODO: Take arguments, debug level
int main (int argc, char* argv[])
{
    SetUpTiming();
    // TODO: Improve RNG?
    uint32_t Seed;  
    Seed = (uint32_t)time(NULL);
    //Seed = 1699659409;
    srand(Seed);
    int DEBUG = 2;

    if (argc < 8)
    {
        printf("Usage: Provide the following values in order (in cm)\nRegion Count (N),\nRegion 1: Source , SigmaA, Diffusion Coefficient, Length, Cells\nRegion 2: ...\n...\nRegion N: ...\nG-S Epsilon.\n");
        exit(1);
    }
    uint32_t RCount = atoi(argv[1]);
    if ((uint32_t) argc < 3+RCount*5)
    {
        printf("Insufficent Arguments\n");
        exit(1);
    }
    
    real32* SAr = malloc(RCount*sizeof(real32));
    real32* AAr = malloc(RCount*sizeof(real32));
    real32* DAr = malloc(RCount*sizeof(real32));
    real32* aAr = malloc(RCount*sizeof(real32)); 
    uint32_t* NAr = malloc(RCount*sizeof(uint32_t)); 
    for (uint32_t i = 0; i < RCount; i++)
    {
        SAr[i] = strtof(argv[i*5+2],0);
        AAr[i] = strtof(argv[i*5+3],0);
        DAr[i] = strtof(argv[i*5+4],0);
        aAr[i] = strtof(argv[i*5+5],0);
        NAr[i] = atoi(  argv[i*5+6]);
    }
    real32 EP    = strtof(argv[2+RCount*5], 0);
    struct RegionDesc RD = {RCount, {SAr, AAr, DAr, aAr}, NAr};

    if (DEBUG) printf("Generating Graph for Multi-Region Example\n");
    if (DEBUG) printf("\tNumerical Values\n");

    struct GSOutput Output = GaussSeidel(RD, EP);

    if (DEBUG) printf("Starting Graphing\n");

    // Graphing parameters
    real32 Length = Output.X[Output.N];
    real32 MaxPhi = 0.f;
    for (uint32_t i = 0; i < Output.N+1; i++)
    {
        MaxPhi = (real32) fmax(MaxPhi, Output.Phi[i]);
    }

    metafl("XWIN");
    winsiz(1500,1060);
    setpag("da4l");
    disini();
    setvlt("RRAIN");
	clrcyc(4,210); // Replace yellow with a blue
    pagera();
    hwfont();
    //axspos (450, 1800);
    //axslen (2200, 1200);
    texmod ("on");
    name   ("Location (cm)", "x");
    name   ("$\\phi(x)$", "y");
    ticks(10,"xy");
    logtic("FULL");
    ticpos("LABELS","xy");

    titlin ("Gauss-Seidel Finite Difference Method", 1);

    setgrf("NAME","NAME","TICKS","TICKS");
    graf   (0.f, Length, 0.f, Length/10.f, 0.f, MaxPhi*1.05f, 0.f, MaxPhi/10.f);

    title();
	chncrv("COLOR");
	//curve(GraphArray[0], GraphArray[1], MaxRows);
    
	incmrk(-1);
	curve(Output.X, Output.Phi, Output.N+1);

    endgrf();
    disfin();

    if (DEBUG) printf("Finished\n");
    if (DEBUG==1) fgetc(stdin);
    return 0;
}

