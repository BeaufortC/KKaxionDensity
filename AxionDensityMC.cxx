#include "MCfunctions.h"

using namespace std;

static void usage(const string &name) {
  cerr << "Usage: " << name << " <delta> <R [keV^{-1}]> <g10>\n"
       << "\tExample: "<< name<< " 2 1e3 9.2e-4\n"
       <<endl;
}

int main(int argc, char** argv){

    AxionMC *axionMC = new AxionMC();

//=============================//
//  INITIALIZATION
//=============================//

    ///Start clock
    chrono::high_resolution_clock::time_point startingTime = chrono::high_resolution_clock::now();

    ///Reading input parameters
    if(argc!=4) {
        cerr<<"<main> Something wrong with the arguments\n"<<endl;
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    axionMC->SetDelta(atof(argv[1]));
    axionMC->SetRCompactification(atof(argv[2]));
    axionMC->SetG10(atof(argv[3]));

    ///Range
    int nMass = 25;//One thread per mass will be created
    double minMass = 1., maxMass = 26.; //keV
    double massWidth = (maxMass - minMass)/nMass;

//=============================//
//  PREPARE THE THREADS
//=============================//

    ROOT::EnableThreadSafety();

    ///Prepare threads
    int numThreads=nMass;
    pthread_t threads[numThreads];
    pthread_attr_t attr; //for joinable threads
    void *status;

    // Initialize and set thread joinable
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    ///Multithreading structure
    struct ArgThread arg[numThreads];

    cout<<"\nSTATUS: Initialization completed. Start the MC loop.\n"<<endl;

//=============================//
//  MC LOOP
//=============================//

    ///Loop over mass -> generate a thread for each mass
    for(int iMass=0 ; iMass < nMass ; iMass++){

        //Fill the argument structure
        arg[iMass].ax = axionMC;
        arg[iMass].threadID = iMass;
        arg[iMass].mass = minMass + double(iMass)*(maxMass - minMass)/nMass;;
        arg[iMass].massWidth = massWidth;

        //Create thread
        if(pthread_create(&threads[iMass], &attr, &ThreadRunner, (void *)&arg[iMass])){
            cerr<<"<ERROR> Unable to create the thread number "<<iMass<<". Exit"<<endl;
            return EXIT_FAILURE;
        }
    }

    //Free attribute and wait for the other threads
    pthread_attr_destroy(&attr);
    for(int iThread = 0 ; iThread < numThreads ; iThread ++){
        if(pthread_join(threads[iThread], &status)){
            cerr<<"<ERROR> Unable to join threads. Exit"<<endl;
            return EXIT_FAILURE;
        }
    }

    cout<<"\nSTATUS: end of the MC loops. The threads have been joined."<<endl;

    chrono::high_resolution_clock::time_point endingTime = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedTime = chrono::duration_cast<chrono::duration<double>>(endingTime - startingTime);
    cout << "Elapsed time: "<<elapsedTime.count()/60.<<" minutes"<<endl;

//=============================//
//  MERGE THE FILES AND GENERATE HISTO
//=============================//

    ///Histo declaration
    TH1F *hDensity = new TH1F("hDensity", ";Distance to the Sun [unit of solar radius];KK axion density [10^{14} cm^{-3}]", axionMC->GetNbinsDensity(), axionMC->GetMinDist(), axionMC->GetMaxDist());
    double scalingFactor = 1e14;

    ///Output file
    axionMC->WriteHeader(nMass, minMass, maxMass);
    ofstream outputFile;
    outputFile.open ("output/density.txt", ofstream::app);

    ///Read the generated files
    double radius, weight;
    for(int i=0 ; i<numThreads ; i++){

        char filename[64];
        sprintf(filename, "output/density_%.0fkeV.txt", arg[i].mass);
        ifstream file(filename);
        string bufferline;

        while (getline(file, bufferline)){
            if(bufferline.size() ==0 || bufferline.at(0) == '#')
                continue;
            istringstream iss(bufferline);
            if(!(iss >> radius >> weight)){
                cerr<<"<main> Cannot properly read the mass files."<<endl;
                return EXIT_FAILURE;
            }//End test non empty line

            outputFile << radius <<"\t"<< weight << "\n";
            hDensity->Fill(radius, weight/scalingFactor);
        }//End loop over file lines
    }

    outputFile.close();

    TCanvas *cDensity = new TCanvas("cDensity");
    cDensity->SetLogy();
    cDensity->SetTicks(1,1);
    cDensity->SetGrid(1,1);

    hDensity->Draw("hist");

    char outputCanvas[128];
    sprintf(outputCanvas, "output/density_delta%1.0f_R%3.0f_g14-%.0f.png", axionMC->GetDelta(), axionMC->GetRcompact(), axionMC->GetG10()*1e4);
    cDensity->SaveAs(outputCanvas);
    char outputHisto[128];
    sprintf(outputHisto, "output/density_delta%1.0f_R%3.0f_g14-%.0f.root", axionMC->GetDelta(), axionMC->GetRcompact(), axionMC->GetG10()*1e4);
    hDensity->SaveAs(outputHisto);

    return 0;
}
