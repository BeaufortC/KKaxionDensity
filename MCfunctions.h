#ifndef MCFUNCTIONS_H_INCLUDED
#define MCFUNCTIONS_H_INCLUDED

//C++
#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <pthread.h>

//ROOT
#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TRandom3.h"

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

//Physical quantities
#define HBAR 6.58212e-19// keV.s
#define PI TMath::Pi()
#define GING10 1e-32// g_{ayy}^2 = 1e-32 * g10^2 [keV^{-2}]
#define SOLARRADIUS 6.9551e10// cm
#define SPEEDLIGHT 2.99792e10// cm/s
#define AGESUN 1.452e17// s


/*!
    \class Contains all functions and parameters to generated the KK axions
    and to transport them via the Equations of Motion.
*/

class AxionMC{

    public:
        AxionMC();
        ~AxionMC();

        ///Set Parameters
        inline void SetDelta(double x){delta = x;};
        inline void SetRCompactification(double x){Rcompact = x;};
        inline void SetG10(double x){g10 = x;};
        inline void SetNiterMC(int x){nIterMC = x;};
        inline void SetNenergy(int x){nEnergy = x;};
        inline void SetTminPower(double x){tMinPower = x;};
        inline void SetTmaxPower(double x){tMaxPower = x;};
        inline void SetMinVelocity(double x){minVelocity = x;};
        inline void SetMaxVelocity(double x){maxVelocity = x;};
        inline void SetNbinsDensity(int x){nBinsDensity = x;};
        inline void SetMinDist(double x){minDist = x;};
        inline void SetMaxDist(double x){maxDist = x;};
        inline void SetBoxWidthFactor(double x){boxWidthFactor = x;};

        ///Get Parameters
        inline int GetNbinsDensity(){return nBinsDensity;};
        inline double GetMinDist(){return minDist;};
        inline double GetMaxDist(){return maxDist;};
        inline double GetDelta(){return delta;};
        inline double GetRcompact(){return Rcompact;};
        inline double GetG10(){return g10;};
        inline std::vector<double> GetVR(){return vRadius;};

        ///Methods
        int MCthread(AxionMC *axMC, int threadID, double mass, double massWidth);
        void ReadSolarModel();
        int EquationsOfMotion(double t, const double y[], double f[]);
        static int EquationsOfMotion_runner(double t, const double y[], double f[], void * const opaque){
            assert(opaque);
            return(static_cast<AxionMC *>(opaque)->EquationsOfMotion(t, y, f));
        }
        std::vector<double> SolveEOM(double y[], double time, AxionMC *ax);
        int PlotEOM(double y[], int nIter, double tFirst, double tLast, AxionMC *ax);
        TH1F CoalRate(double mass, double energy, int ID);
        void WriteHeader(int nMass, double minMass, double maxMass);


    private:
        ///Structures
        struct Axion{
            double theta0;//angle perpendicular to the plane of motion ; in radian
            double phi0;//angle in the plane of motion (x-y) ; in radian
            double x0, y0;//Initial position at production in the Sun ; in unit of solar radius
            double vx0, vy0;//Initial velocity at production in unit of sqrt{GM/R} = sqrt{2} * v_esc
                    //where v_esc is the escape velocity at the Sun's surface
            double time;//time at which we look at the position of the axion ; in unit of omega^{-1} = sqrt{R^3/(GM)}
            double xt, yt, zt, rt;//positions at t=time ; in unit of solar radius
        };

        ///Attributes
        double delta;//Number of extra dimensions
        double Rcompact;//Radius of compactification, in keV^{-1}
        double g10;//Coupling to photon, g10 = 10^{-10} GeV^{-1} * g_{ayy}
        int nIterMC;//number of iterations per mass and per energy ;
        int nEnergy;//number of energy points -> this is a key parameter to be accurate enough to get trapped axions at large distances ;
        double tMinPower;//10^(tMinPower) = Minimum time of the integration of the EOM
        double tMaxPower;//10^(tMaxPower) = Maximal time of the integration of the EOM
        double minVelocity;//range of initial velocity of interest in EOM units
        double maxVelocity;//range of initial velocity of interest in EOM units
        int nBinsDensity;//output histo binning
        double minDist;//output histo first bin
        double maxDist;//output histo last bin
        double boxWidthFactor;//Scaling of the box in which the axions must be located

        std::vector<double> vRadius;//Solar model, vector of radius inside the Sun ; Unit of solar radius
        std::vector<double> vTemp;//Solar model, vector of temperature inside the Sun ; in keV
};


///Structure for passing the arguments to the multithreads
struct ArgThread{
    AxionMC *ax;
    int threadID;
    double mass;
    double massWidth;
};

///Helper to run the threads
inline void* ThreadRunner(void *voidArg) {
    ArgThread *args = (ArgThread*)voidArg;
    int res = args->ax->MCthread(args->ax, args->threadID, args->mass, args->massWidth);
    return new int(res); // Return an `int` pointer to pass back thread runner's results
};

#endif //MCFUNCTIONS_H_INCLUDED
