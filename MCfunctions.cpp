#include "MCfunctions.h"

using namespace std;


AxionMC::AxionMC(){

    //Initialize the solar model vectors
    ReadSolarModel();

    ///Default values
    nIterMC = 2e1;//number of iterations per mass and per energy ;
    nEnergy = 1e6;//number of energy points -> this is a key parameter to be accurate enough to get trapped axions at large distances ;
    tMinPower = 3.;//10^(tMinPower) = Minimum time of the integration of the EOM
    tMaxPower = 4.;//10^(tMaxPower) = Maximal time of the integration of the EOM
    boxWidthFactor = 2e11;//Scaling of the box in which the axions must be located
    minVelocity = 0.8*sqrt(2);
    maxVelocity = 1.4*sqrt(2);//range of initial velocity of interest in EOM units
        //Reminder: v_escp = sqrt(2GM/R) = sqrt(2) in EOM units (at the Sun's surface)
    nBinsDensity = 299;//Number of bins in the output histo
    minDist = 1.;//First bin
    maxDist= 300.;//Last bin

}

///--------------------------///

AxionMC::~AxionMC(){

}


///--------------------------///

    /*!
        \brief MC loop.
        For a given mass, loop over the energy, generate the axions,
        solve the EOM, and check if the axion is trapped within a box.
        The info of each trapped axion is written in an individual output file
    */


int AxionMC::MCthread(AxionMC* axMC, int threadID, double mass, double massWidth){

    ///Factors
    double highDimVolume = 2 * pow(PI, delta/2.) / TMath::Gamma(delta/2.) * pow(Rcompact, delta);
        //prefactor when counting the number of KK modes
    double prefactorWidth = GING10 * pow(g10, 2) / (64. * PI * HBAR);
        //Width = g^2 m^3/(64 pi hBar). With this prefactor we can express m in keV
    double cInUnitEOM  = 686.32;
        //c = 686.32 * sqrt(GM/R) ; where sqrt(GM/R) is the velocity unit in the EOM solver

    ///Random number generator
    TRandom3 *rand = new TRandom3(time(NULL) * threadID);

    ///Output file for this thread
    char filename[64];
    sprintf(filename, "output/density_%.0fkeV.txt", mass);
    ofstream outFileThread;
    outFileThread.open(filename);

    ///Declaration
    struct Axion ax;
    double r0, v0;//initial radius and velocity at production
    double initial_y[4];//Initial conditions for the EOM solver
    vector<double> vPosAtT;//X and T at a given time

    double energy, energyWidth, velocityWidth;
    double modeMultiplicity;//KK mode multiplicity for a given mass
    double decayWidth;//Decay width of the axion mode in s^{-1}
    double amountTrapped;//Gives the amount of trapped KK axions of mass m accumulated over Sun history
    double integralRate;//in s^{-1} ; For a given mass, integrate the production rate over the distance in the Sun and over dE
    double weight;//Weighting factor to convert one entry in the histogram into the correct physical number density

    double boxWidth = boxWidthFactor/SOLARRADIUS;//we do not set 1cm directly because the statistics would be too low
    int inBox;//Number of trapped axions that are located in a box of dimension 1BoxWidth * 1BoxWidth * 1bin
    double progress = 0.1;//To print the progress

    modeMultiplicity = highDimVolume * pow(mass, delta - 1.) * massWidth;
        //We multiply by the massWidth so that the loop over masses can be seen as an integral over the masses
    decayWidth = prefactorWidth * pow(mass, 3);
    amountTrapped = modeMultiplicity * (1. - exp(-AGESUN * decayWidth)) / decayWidth;

    ///Loop over energy
    for(int iEnergy=0 ; iEnergy < nEnergy ; iEnergy ++){

        //We only consider the energy range for which v ~ v_escp
        v0 = minVelocity + double(iEnergy)*(maxVelocity - minVelocity)/nEnergy;//in EOM units
        velocityWidth = (maxVelocity - minVelocity) / nEnergy;
        energy = mass * sqrt(1 + pow(v0/cInUnitEOM, 2));//Semi-classical approximation
        energyWidth = pow(mass / cInUnitEOM, 2) * v0 / energy * velocityWidth;
            //dE = m^2/c^2 * v/E * dv

        TH1F hRate = CoalRate(mass, energy, threadID);//Production rate histogram
        integralRate = hRate.Integral();//in s^{-1}.keV^{-1}
        integralRate *= energyWidth;//in s^{-1}. So that the loop over E is the Riemannian integral

        inBox = 0;//Reinit

        ///Loop over produced axions
        for(Long64_t iAxions = 0 ; iAxions < nIterMC ; iAxions++){

            r0 = hRate.GetRandom();//radius of production

            ///Generate random axions parameters
            ax.theta0   = rand->Uniform(0., PI);//isotrope production
            ax.phi0     = rand->Uniform(0., 2*PI);//isotrope production
            ax.x0       = r0 * cos(ax.phi0);
            ax.y0       = r0 * sin(ax.phi0);
            ax.vx0      = v0 * cos(ax.phi0);
            ax.vy0      = v0 * sin(ax.phi0);
            ax.time     = pow(10, rand->Uniform(tMinPower, tMaxPower));

            ///Solve EOM
            initial_y[0] = ax.x0; initial_y[1] = ax.vx0;
            initial_y[2] = ax.y0; initial_y[3] = ax.vy0;
            vPosAtT = SolveEOM(initial_y, ax.time, axMC);
            if(vPosAtT.size() !=2)
                cerr<<"<AxionMC::MCthread> Error when solving the EOM"<<endl;

            ax.xt = vPosAtT.front();
            ax.yt = vPosAtT.back();
            ax.rt = sqrt(pow(ax.xt, 2) + pow(ax.yt, 2));
            ax.zt = ax.rt * cos(ax.theta0);

            ///Test is the axion is within a 1Rsun * 1Rsun * 1bin box
                //Due to isotropy we can locate this box anywhere in y and z
            if(abs(ax.yt) < boxWidth/2. && abs(ax.zt) < boxWidth/2. && ax.rt > 1 && ax.xt > 0){

                weight = amountTrapped * integralRate;
                weight /= nIterMC;
                weight /= pow(boxWidthFactor,2) * (maxDist - minDist)/nBinsDensity * SOLARRADIUS;//1Box*1Box*1bin -> cm^3

                inBox ++;

                outFileThread << ax.rt <<"\t"<< weight<<"\n";
            }


        }//End loop over axions

        ///Print progress
        if(threadID == 0 && progress < double(iEnergy)/nEnergy*100){
            cout<<progress<<"% done"<<endl;
            if(abs(progress - 99.8) < 1e-5)
                cout<<"End of the first thread. Waiting for the others to finish."<<endl;
            progress += 0.1;
        }

    }//End loop over energy

    outFileThread.close();

    return 0;
}


///--------------------------///

/*!
    \brief Read the solar model and store the data into vectors.
    We use the Saclay model, cf:
    "Solar Neutrino Emission Deduced from a Seismic Model", Turck-Chieze et al.,"
*/

void AxionMC::ReadSolarModel(){

    ///Read the Saclay solar model
    ifstream file("data/solarModel.dat");
    double rTmp, Ttmp;
    string bufferline;

    while (getline(file, bufferline)){

        if(bufferline.size() ==0 || bufferline.at(0) == '#')
            continue;
        istringstream iss(bufferline);

        if(!(iss >> rTmp >> Ttmp))
            cerr<<"<AxionMC::ReadSolarModel> Cannot properly read the Saclay data."<<endl;

        vRadius.push_back(rTmp);
        vTemp.push_back(Ttmp);

    }//End loop over file lines
    file.close();
}

///--------------------------///

   /*!
        \brief Equations of motion of a point-like particle
            in a classical gravitational field
        - Distance in unit of Solar Radius
        - Time in unit of omega^{-1} = sqrt{R^3/(GM)}

        We separate the cases inside and outside the Sun
    */

int AxionMC::EquationsOfMotion(double t, const double y[], double f[]){

    (void)(t); //avoid unused parameter warning
    double radius = sqrt(y[0]*y[0] + y[2]*y[2]);

    //Case r<R
    if(radius <= 1){
        //We use phi(r) = -0.5*w2*R^2*(3 - (r/R)^2)
        f[0] = y[1]; // dx/dt = u
        f[1] = - y[0]; //du/dt = -w2*x
        f[2] = y[3]; // dy/dt = v
        f[3] = - y[2]; //dv/dt = -w2*y
    }
    //Case r>R
    else{
        f[0] = y[1]; // dx/dt = u
        f[1] = - y[0] / pow(radius, 3); //du/dt = -w2*x/r3
        f[2] = y[3]; // dy/dt = v
        f[3] = - y[2] / pow(radius, 3); //dv/dt = -w2*y/r3
    }

    return GSL_SUCCESS;
}

///--------------------------///
    /*!
        \brief Solve the EOM for a given time
        Several solvers are available (GSL solvers)
    */

vector<double> AxionMC::SolveEOM(double y[], double time, AxionMC *ax){

    gsl_odeiv2_system sys = {&AxionMC::EquationsOfMotion_runner, NULL, 4, reinterpret_cast<void *>(ax)};

    //Solver Runge-Kutta-Fehlberg (4, 5)
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 1e-8);


    vector<double> vOut;
    double initialT = 0.;
    if(gsl_odeiv2_driver_apply (driver, &initialT, time, y) != GSL_SUCCESS){
        //Initial conditions at t=initialT ; evolve up to t=time
        cerr<<"<AxionMC::SolveEOM> Cannot integrate the equations of motion. EXIT."<<endl;
        return vOut;
    }

    vOut.push_back(y[0]);
    vOut.push_back(y[2]);

    gsl_odeiv2_driver_free (driver);

    return vOut;
}


///--------------------------///

    /*!
        \brief Solve and plot the equations of motions
        for a given set of intial conditions.
        Several solvers are available (GSL solvers)

        \note Can be used for debugging but is not called in the main program
    */

int AxionMC::PlotEOM(double y[], int nIter, double tFirst, double tLast, AxionMC *ax){

    gsl_odeiv2_system sys = {&AxionMC::EquationsOfMotion_runner, NULL, 4, reinterpret_cast<void *>(ax)};

    //Solver Runge-Kutta-Fehlberg (4, 5)
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 1e-8);

    vector<double> vX, vY;
    double initialT = tFirst;
    for(int i = 1; i <= nIter; i++){
        double ti = i * (tLast - tFirst) / nIter + initialT;

        if(gsl_odeiv2_driver_apply (driver, &tFirst, ti, y) != GSL_SUCCESS){
            cerr<<"<plotEOM> Cannot integrate the equations of motion. EXIT."<<endl;
            return EXIT_FAILURE;
        }

        vX.push_back(y[0]);
        vY.push_back(y[2]);
    }

    gsl_odeiv2_driver_free (driver);

    TCanvas *cEOM = new TCanvas("cEOM");
    TGraph *g = new TGraph(vX.size(), &(vX[0]), &(vY[0]));
    g->SetTitle("Motion of the axion around the Sun; X [unit of solar radius]; Y [unit of solar radius]");
    g->Draw();

    cEOM->WaitPrimitive();

    return EXIT_SUCCESS;
}


///--------------------------///
    /*!
        \brief Integration of the inverse decay * blackbody spectrum
        over a solar model to obtain the differential rate of axion production from the coalescence of
        2 plasmons in the Sun. Differential rate in s^{-1}.keV^{-1}.

        \arg ID is the ID of the thread

        \return the histo of the differential rate.
    */

TH1F AxionMC::CoalRate(double mass, double energy, int ID){

    //Prepare histo
    char hName[64];
    sprintf(hName, "hRate_%d", ID);
    TH1F hRate(hName, "Axion production rate from coalescence; Distance to the center [in solar radius]; Production rate [s^{-1}.keV^{-1}]", vRadius.size(), vRadius.front(), vRadius.back());
    double rate;
    double step = (vRadius.back() - vRadius.front())/vRadius.size();

    ///Integrate the production over the solar model
    for(int i=0 ; i<vRadius.size() ; i++){

        rate = GING10 * pow(g10, 2) /(128. * pow(PI, 3));//prefactor
        rate /= pow(HBAR, 4) * pow(SPEEDLIGHT, 3);//Unit conversion
        rate *= pow(mass, 4) * sqrt(pow(energy, 2) - pow(mass, 2));//energy and mass
        rate /= (exp(energy / vTemp.at(i)) -1);//Bose-Einstein distribution
        rate *= 4. * PI * pow(vRadius.at(i) * SOLARRADIUS, 3);// integral over dr^3 = 4 pi r^2 dr
        rate *= step;//binWidth

        hRate.SetBinContent(i+1, rate);
    }//end integration solar model

    return hRate;
}

///--------------------------///
    /*!
        \brief Write the simulation parameter in the header of the output file
    */

void AxionMC::WriteHeader(int nMass, double minMass, double maxMass){

    ofstream file;

    file.open("output/density.txt");
    file <<"#PARAMETERS:\n";
    file << "#delta: "<<delta<<"\t R: "<<Rcompact << "\tg10: "<<g10<<"\n";
    file << "#nIterMC: "<<nIterMC<<"\tnEnergy: "<<nEnergy<<"\tnMass: "<<nMass<<"\n";
    file << "#Mass range: ["<<minMass<<","<<maxMass<<"]\tVelocity range: ["<<minVelocity<<","<<maxVelocity<<"]\tTime range: ["<<tMinPower<<","<<tMaxPower<<"]\n";
    file << "#NBinsdensity: "<<nBinsDensity<<"\tminDist: "<<minDist <<"\tmaxDist: "<<maxDist<< "\tboxWidthFactor: "<<boxWidthFactor  <<"\n";
    file << "#\n#Distance to the Sun [solar radius]\t weight\n";

    file.close();
}
