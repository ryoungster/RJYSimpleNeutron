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
    if (DEBUG) printf("Generating Graph for Multi-Region Example\n");

    uint32_t RCount = 2;
    real32 SAr[] = {1000.0f, 0.0f};
    real32 cAr[] = {1.f - .1f*1.7f*3.f, 1.f - .01f*3.f};
    real32 TAr[] = {1.f/1.7f/3.f, 1.f/1.0f/3.f}; // Conversions atcompile time for example
    real32 aAr[] = {30.0f, 15.0f}; 
    uint32_t NAr[] = {30, 10}; 
    real32 Hist    = 1e7f;
    struct RegionDesc RD = {RCount, {SAr, cAr, TAr, aAr}, NAr};

    if (DEBUG) printf("\tNumerical Values\n");
    struct MCOutput Output = MonteCarlo(RD, Hist);

    if (DEBUG) printf("\tAnalytical Values\n");

    real32 AAr[] = {.1f, .01f};
    real32 DAr[] = {1.7f, 1.0f};
	const real32 FL = sqrt(DAr[0]/AAr[0]);
	const real32 RL = sqrt(DAr[1]/AAr[1]);
    const real32 RC = (1.0f+2.0f*DAr[1]/RL*tanh(aAr[1]/RL))/(2.0f*DAr[1]/RL+tanh(aAr[1]/RL));
    const real32 FC = 1/(DAr[0]/FL*RL/DAr[1]/RC*sinh(aAr[0]/FL)+cosh(aAr[0]/FL));
    const real32 FA = -SAr[0]/AAr[0]*FC;
    const real32 RA = FA*cosh(aAr[0]/FL)+SAr[0]/AAr[0];

	for (uint32_t  i = 0; i < MaxRows; i++)
    {
        GraphArray[0][i] = ((real32)i)/((real32)(MaxRows-1))*45.f;
        if (GraphArray[0][i] < aAr[0])
        {
            GraphArray[1][i]=SAr[0]/AAr[0]*(1-cosh(GraphArray[0][i]/FL)*FC);
        }
        else
        {
            GraphArray[1][i]=RA * (cosh((GraphArray[0][i]-aAr[0])/RL) - RC * sinh((GraphArray[0][i]-aAr[0])/RL));
        }
    }
    /*
	if (DEBUG) printf("\tError Values\n");
	for (uint32_t i = 0; i < GSLEN; i++)
	{
		uint32_t N = GSVALUES[i];
		GSErrorArray[i]=0; 
		for (uint32_t n = 0; n < N+1; n++) {
			real32 Base = VacConst*(1.f-cosh(GSXArrays[i][n]/L)*VacDiv);
			real32 Error = fabs(Base-GSPhiArrays[i][0][n])/Base;
			GSErrorArray[i] = fmax(GSErrorArray[i], Error);
		}
	}
	
	printf("Displaying Data Tables\n");
    printf("N\tError\n");

    for (uint32_t i = 0; i < GSLEN; i++)
    {
        
        printf("%i\t&%.3e\\\\\n", GSVALUES[i], GSErrorArray[i]);
        //printf("%i&\t%.3f\\\\\n", N, TrapArrays[i][1][N]);
    }
	*/

    
    if (DEBUG) printf("Starting Graphing\n");

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

    titlin ("MonteCarlo Method", 1);

    setgrf("NAME","NAME","TICKS","TICKS");
    graf   (0.f, 45.f, 0.f, 45.f/10.f, 0.f, 11000.f, 0.f, 2000.f);

    title();
	chncrv("COLOR");
	curve(GraphArray[0], GraphArray[1], MaxRows);
    
	incmrk(-1);
	curve(Output.X, Output.Phi, Output.N+1);

    endgrf();
    disfin();

    if (DEBUG) printf("Finished\n");
    if (DEBUG==1) fgetc(stdin);
    return 0;
}

