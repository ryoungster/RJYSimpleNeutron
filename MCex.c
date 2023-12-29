#include <stdio.h>
#include <stdlib.h>
#include <math.h>                  
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "dislin.h"
#include "Neutronic.h"



const uint32_t MCVALUES[] = {1e4,1e5,1e6};
#define NBIN 100
#define MCLEN ArrayCount(MCVALUES)
uint32_t MCBinArrays[MCLEN][NBIN];
real32 MCPhiArrays[MCLEN][NBIN];
real32 MCXArray[NBIN];
//real32 GSErrorArray[GSLEN];
//real32 TTimeArray[TLEN];

#define MaxRows 20000
real32 GraphArray[2][MaxRows];


#define a 	10.0f
#define SA 	0.1f
#define SS	0.2f
#define s0	1000.0f
#define EP	1e-7f
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
    if (DEBUG) printf("Generating Graph for HW3-4-2\n");

    if (DEBUG) printf("Generating Data\n");
	const real32 DX = a/(real32)NBIN;
	if (DEBUG) printf("\tSetting-Up Arrays\n");
	for (uint32_t i = 0; i < NBIN; i++)
    {
		MCXArray[i] = (real32)i*DX+DX/2.f;
    }
    
	
	if (DEBUG) printf("\tPerforming Monte Carlo\n");
	const real32 DivST = 1.f/(SA+SS);
	const real32 DivDX = 1.f/DX;
#define SpawnX()	(randReal32()*a)
#define IsoMu()		(randReal32()*2.f-1.f)
#define DColl()		(-1.f*log(randReal32())*DivST)
	for (uint32_t i = 0; i < MCLEN; i++)
    {
		uint32_t N = MCVALUES[i];

		for (uint32_t n = 0; n < N; n++)
		{
			real32 X = SpawnX();		
			real32 Mu = IsoMu();
			for(;;)
			{

				X += Mu*DColl();
				if (X <= 0.f)
				{
					Mu = -Mu;
					X = 0.f;

				} 
				else if (X < a)
				{
					MCBinArrays[i][(uint32_t)floor(X*DivDX)]++;
					if (randReal32()>(SA*DivST))
					{
						Mu = IsoMu();
					}
					else
					{
						break;
					}
				}
				else
				{
					break;
				}
			}

			
		}
		for (uint32_t j = 0; j < NBIN; j++)
		{
			MCPhiArrays[i][j] = (real32)MCBinArrays[i][j]*a*s0*DivDX/(real32)N*DivST;
		}
	}

    if (DEBUG) printf("\tAnalytical Values\n");
	const real32 D = 1/(3*(SA+SS));
	const real32 L = sqrt(D/SA);
	const real32 VacConst = s0/SA;
	const real32 VacDiv = 1.0f/(cosh(a/L)+2.f*D/L*sinh(a/L));
	for (uint32_t  i = 0; i < MaxRows; i++)
    {
        GraphArray[0][i] = ((real32)i)/((real32)(MaxRows-1))*a;
        GraphArray[1][i]=VacConst*(1-cosh(GraphArray[0][i]/L)*VacDiv);
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

    if (DEBUG > 1) 
    {
        metafl("XWIN");
    } 
    else
    {
#define OutFile "../texImage/NSE451-HW3-4-2.EPS"
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

    titlin ("Monte Carlo Method, $a=5$", 1);

    setgrf("NAME","NAME","TICKS","TICKS");
    graf   (0.f, a, 0.f, 1.f, 0.f, 10000.f, 0.f, 2000.f);

    title();
	chncrv("COLOR");
	curve(GraphArray[0], GraphArray[1], MaxRows);
    
	incmrk(-1);
	curve(MCXArray,MCPhiArrays[0],NBIN);
	curve(MCXArray,MCPhiArrays[1],NBIN);
	curve(MCXArray,MCPhiArrays[2],NBIN);
	

	color("FORE");
	
    char LegendBuffer[300];
    legini(LegendBuffer, 4, 30);

    leglin(LegendBuffer, "Analytical Solution", 1);
    leglin(LegendBuffer, "N=1e4", 2);
    leglin(LegendBuffer, "N=1e5", 3);
	leglin(LegendBuffer, "N=1e6", 4);

	legtit("Legend");
    //legbgd(0);
    legend(LegendBuffer, 3);
	
    endgrf();
    disfin();
	
    if (DEBUG) printf("Finished\n");
    if (DEBUG==1) fgetc(stdin);
    return 0;
}

