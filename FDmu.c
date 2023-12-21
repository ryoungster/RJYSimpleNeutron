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
    if (DEBUG) printf("Generating Graph for HW3-3\n");

    uint32_t RCount = 2;
    uint32_t NAr[] = {30, 10}; 
    real32 SAr[] = {1000.0f, 0.0f};
    real32 AAr[] = {.1f, .01f};
    real32 DAr[] = {1.7f, 1.0f};
    real32 aAr[] = {30.0f, 15.0f}; 
    real32 EP    = 1e-7f;

    real32* GSTriAnsArrays[4];
    real32* GSPhiArrays[2];
    real32* GSXArrays;

    uint32_t N = 0;
    real32 Length = 0;
    for (uint32_t i = 0; i < RCount; i++)
    {
        N += NAr[i];
        Length += aAr[i];
        aAr[i] /= (real32)NAr[i];
    }

    
    if (DEBUG) printf("Allocing Arrays\n");
    GSAlloc(GSTriAnsArrays, GSPhiArrays, &GSXArrays, N);

    if (DEBUG) printf("Generating Data\n");

    // Format:
    // Source, SigmaA, D, deltaX
    real32* RegionInfo[] = {SAr, AAr, DAr, aAr};


	if (DEBUG) printf("\tSetting-Up Arrays\n");

    GSInitVR(GSTriAnsArrays, GSXArrays, N, RegionInfo, NAr, RCount);
	
	if (DEBUG) printf("\tPerforming Iteration\n");
    real32 Convergence;
    do
    {
        Convergence = GSStep(GSTriAnsArrays, GSPhiArrays, N);
    }
    while (Convergence > EP);

    /*
    if (DEBUG) printf("\tAnalytical Values\n");
	const real32 L = sqrt(D/SA);
	const real32 VacConst = s0/SA;
	const real32 VacDiv = 1.0f/(cosh(a/L)+2.f*D/L*sinh(a/L));
	for (uint32_t  i = 0; i < MaxRows; i++)
    {
        GraphArray[0][i] = ((real32)i)/((real32)(MaxRows-1))*a;
        GraphArray[1][i]=VacConst*(1-cosh(GraphArray[0][i]/L)*VacDiv);
    }
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

    if (1)//(DEBUG > 1) 
    {
        metafl("XWIN");
    } 
    else
    {
#define OutFile "../texImage/NSE451-HW3-3.EPS"
        metafl("EPS");
        remove(OutFile);
    }
    winsiz(1500,1060);
    setpag("da4l");
    setfil(OutFile);
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
    graf   (0.f, Length, 0.f, 20.f, 0.f, 11000.f, 0.f, 2000.f);

    title();
	chncrv("COLOR");
	//curve(GraphArray[0], GraphArray[1], MaxRows);
    
	incmrk(-1);
	curve(GSXArrays,GSPhiArrays[0],N+1);

	color("FORE");
	
    char LegendBuffer[300];
    legini(LegendBuffer, 4, 30);

    leglin(LegendBuffer, "Analytical Solution", 1);
    leglin(LegendBuffer, "N$=10$", 2);
    leglin(LegendBuffer, "N$=40$", 3);
	leglin(LegendBuffer, "N$=100$", 4);

	legtit("Legend");
    //legbgd(0);
    //legend(LegendBuffer, 3);
	
    endgrf();
    disfin();

    if (DEBUG) printf("Finished\n");
    if (DEBUG==1) fgetc(stdin);
    return 0;
}

