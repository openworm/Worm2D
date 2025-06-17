// =============================================================
// Evolving Locomotion
// Dec 2021
// Eduardo Izquierdo
// =============================================================



#include <iostream>
#include <iomanip>
#include <math.h>
#include "TSearch.h"
#include "VectorMatrix.h"
//#include "Worm.h"
//#include <vector>

#include "jsonUtils.h"

#include "Mainvars.h"
#include <stdio.h>
#include <string.h>
//#include <sys/stat.h>
#include "argUtils.h"

SuppliedArgs supArgs1;

#define PRINTTOFILE

int skip_steps = 10;
using namespace std;

// Integration parameters
const int Duration = 24;
const double Transient = 8.0;
const double StepSize = 0.005;
const int N_curvs = 23;

// Fitness traj
const double    AvgSpeed = 0.0001; //0.00022;              // Average speed of the worm in meters per seconds
const double    BBCfit = AvgSpeed*Duration;

// Genotype -> Phenotype Mapping Ranges
const double    BiasRange               = 16.0; //15.0;
const double    SCRange                 = 16.0; //15.0;
const double    CSRange                 = 16.0; //15.0;
const double    ESRange                 = 2.0;
const double    SRmax                   = 200.0;
const double    NMJmax                  = 0.8; //1.2;
const double    NMJmin                  = 0.0;

const int SR_A = 1;
const int SR_B = 2;




// Size of genotype
int	VectSize = 17;

// ------------------------------------
// Genotype-Phenotype Mapping
// ------------------------------------
void GenPhenMapping(TVector<double> &gen, TVector<double> &phen)
{
  // Parameters for the Stretch Receptors
  phen(SR_A) = MapSearchParameter(gen(SR_A), 0.0, SRmax);
  phen(SR_B) = MapSearchParameter(gen(SR_B), 0.0, SRmax);

  // Bias
  int k=3;
  for (int i = 1; i <= 3; i++){
    phen(k) = MapSearchParameter(gen(k), -BiasRange, BiasRange);k++;
  }
  // Self connections
  for (int i = 1; i <= 3; i++){
    phen(k) = MapSearchParameter(gen(k), -SCRange, SCRange);k++;
  }
  // DA, DB, VA, VB Chemical synapses (excitatory)
  for (int i = 1; i <= 2; i++){
    phen(k) = MapSearchParameter(gen(k), 0.0, CSRange);k++;
  }
  // VD Chemical synapses (Inhibitory)
  for (int i = 1; i <= 2; i++){
    phen(k) = MapSearchParameter(gen(k), -CSRange, 0.0);k++;
  }
  // Interunits Gap junctions
  for (int i = 1; i <= 2; i++){
    phen(k) = MapSearchParameter(gen(k), 0.0, ESRange);k++;
  }
  // Excitatory NMJ Weight
  for (int i = 1; i <= 2; i++){
    phen(k) = MapSearchParameter(gen(k), NMJmin, NMJmax);k++;
  }
  // Inhibitory NMJ Weight
  for (int i = 1; i <= 1; i++){
    phen(k) = MapSearchParameter(gen(k), -NMJmax, -NMJmin);k++;
  }
}

// ------------------------------------
// Fitness function
// ------------------------------------
double Evaluation(TVector<double> &v, RandomState &rs, int direction){
    double fitA,fitB;
    double bodyorientation, anglediff;
    double movementorientation, distancetravelled = 0, temp;
    double distance;
    double xt, xtp, oxt, fxt;
    double yt, ytp, oyt, fyt;

    // Genotype-Phenotype Mapping
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(v, phenotype);
    Worm w(phenotype, 1);
    w.InitializeState(rs);

    if (direction == 1){
        w.AVA_output =  0.0;
        w.AVB_output =  1.0;
    }
    else{
        w.AVA_output =  1.0;
        w.AVB_output =  0.0; // Command Interneuron Activation Backward
    }

    // Transient
    for (double t = 0.0; t <= Transient; t += StepSize){
        w.Step(StepSize, 1);
    }
    xt = w.CoMx(); yt = w.CoMy();
    oxt = w.CoMx(); oyt = w.CoMy();
    // Run
    for (double t = 0.0; t <= Duration; t += StepSize) {
        w.Step(StepSize, 1);
        // Current and past centroid position
        xtp = xt; ytp = yt;
        xt = w.CoMx(); yt = w.CoMy();
        // Integration error check
        if (isnan(xt) || isnan(yt) || sqrt(pow(xt-xtp,2)+pow(yt-ytp,2)) > 10*AvgSpeed*StepSize) {return 0.0;}
        // Velocity Fitness
        bodyorientation = w.Orientation();                  // Orientation of the body position
        movementorientation = atan2(yt-ytp,xt-xtp);         // Orientation of the movement
        anglediff = movementorientation - bodyorientation;  // Check how orientations align
        if (direction == 1){
            temp = cos(anglediff) > 0.0 ? 1.0 : -1.0;           // Add to fitness only movement forward
        }
        else{
            temp = cos(anglediff) > 0.0 ? -1.0 : 1.0;           // Add to fitness only movement backward
        }
        distancetravelled += temp * sqrt(pow(xt-xtp,2)+pow(yt-ytp,2));
    }
    fxt = w.CoMx(); fyt = w.CoMy();
    distance = sqrt(pow(oxt-fxt,2)+pow(oyt-fyt,2));
    fitA = 1 - (fabs(BBCfit - distance)/BBCfit);
    fitA = (fitA > 0)? fitA : 0.0;

    fitB = 1 - (fabs(BBCfit-distancetravelled)/BBCfit);
    fitB = (fitB > 0)? fitB : 0.0;
    return fitB;
}

// ------------------------------------
// Fitness Function
// ------------------------------------
double EvaluationFunction(TVector<double> &v, RandomState &rs){
  double sra = v(SR_A);
  double srb = v(SR_B);
  double fitnessForward, fitnessBackward;
  v(SR_A)= -1.0;
  v(SR_B)= srb;
  fitnessForward = Evaluation(v, rs, 1);
  //  v(SR_A)= sra;
  //  v(SR_B)= -1.0;
  //  fitnessBackward = Evaluation(v, rs, -1);
  //  return (fitnessForward + fitnessBackward)/2;
  return fitnessForward;
  // return fitnessBackward;
}





// ------------------------------------
// Plotting
// ------------------------------------
double save_traces(TVector<double> &v, RandomState &rs){
    ofstream curvfile(supArgs1.rename_file("curv.dat"));
    ofstream bodyfile(supArgs1.rename_file("body.dat"));
    ofstream actfile(supArgs1.rename_file("act.dat"));
    // Genotype-Phenotype Mapping
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(v, phenotype);
    double sra = phenotype(SR_A);
    double srb = phenotype(SR_B);
    
    
    Worm w(phenotype, 1);
   /*  {
    ofstream phenfile(rename_file("phenotype.dat"));
    w.DumpParams(phenfile);
    phenfile.close();
    } */
    
    w.InitializeState(rs);
    w.sr.SR_A_gain = 0.0;
    w.sr.SR_B_gain = srb;
    w.AVA_output =  w.AVA_inact;
    w.AVB_output =  w.AVB_act;

  

    for (double t = 0.0; t <= Transient + supArgs1.traceDuration; t += StepSize){
        w.Step(StepSize, 1);
        w.DumpBodyState(bodyfile, skip_steps);
        w.DumpCurvature(curvfile, skip_steps);
        w.DumpActState(actfile, skip_steps);
    }

     w.sr.SR_A_gain = 0.0;
     w.sr.SR_B_gain = 0.0;

     for (double t = 0.0; t <= (12); t += StepSize){
         w.Step(StepSize, 1);
         w.DumpBodyState(bodyfile, skip_steps);
         w.DumpCurvature(curvfile, skip_steps);
         w.DumpActState(actfile, skip_steps);
     }

     w.sr.SR_A_gain = sra;
     w.sr.SR_B_gain = 0.0;
     w.AVA_output =  w.AVA_act;
     w.AVB_output =  w.AVB_inact;

     for (double t = 0.0; t <= (20); t += StepSize){
         w.Step(StepSize, 1);
         w.DumpBodyState(bodyfile, skip_steps);
         w.DumpCurvature(curvfile, skip_steps);
         w.DumpActState(actfile, skip_steps);
     }

     w.sr.SR_A_gain = 0.0;
     w.sr.SR_B_gain = 0.0;

     for (double t = 0.0; t <= (12); t += StepSize){
         w.Step(StepSize, 1);
         w.DumpBodyState(bodyfile, skip_steps);
         w.DumpCurvature(curvfile, skip_steps);
         w.DumpActState(actfile, skip_steps);
     }

    
    bodyfile.close();
    curvfile.close();
    actfile.close();
    return 0;
}

// ------------------------------------
// Display functions
// ------------------------------------
void EvolutionaryRunDisplay(int Generation, double BestPerf, double AvgPerf, double PerfVar)
{
    cout << BestPerf << " " << AvgPerf << " " << PerfVar << endl;
}

void ResultsDisplay(TSearch &s)
{
    TVector<double> bestVector;
    ofstream BestIndividualFile;
    bestVector = s.BestIndividual();
    BestIndividualFile.open(supArgs1.rename_file("best.gen.dat"));
    BestIndividualFile << setprecision(32);
    BestIndividualFile << bestVector << endl;
    BestIndividualFile.close();
}




// ------------------------------------
// The main program
// ------------------------------------
int main (int argc, const char* argv[])
{
    std::cout << std::setprecision(10);
    long randomseed = static_cast<long>(time(NULL));
   
    if (argc==2) randomseed += atoi(argv[1]);

    if (argc>2) if (!supArgs1.setArgs(argc,argv,randomseed)) return 0;
 
    InitializeBodyConstants();

    if (supArgs1.do_evol){

    supArgs1.doOrigNS = 1;

    TSearch s(VectSize);

    // save the seed to a file

    ofstream seedfile;
    seedfile.open(supArgs1.rename_file("seed.dat"));
    seedfile << supArgs1.randomseed << endl;
    seedfile.close();

    supArgs1.writeMessage();


    // configure the search
    s.SetRandomSeed(supArgs1.randomseed);
    s.SetPopulationStatisticsDisplayFunction(EvolutionaryRunDisplay);
    s.SetSearchResultsDisplayFunction(ResultsDisplay);
    s.SetSelectionMode(RANK_BASED);               //{FITNESS_PROPORTIONATE,RANK_BASED}
    s.SetReproductionMode(GENETIC_ALGORITHM);	    // {HILL_CLIMBING, GENETIC_ALGORITHM}
    s.SetPopulationSize(supArgs1.pop_size);
    s.SetMaxGenerations(supArgs1.max_gens);
    s.SetMutationVariance(0.05);                   // For 71 parameters, an estimated avg change of 0.25 for weights (mapped to 15).
    s.SetCrossoverProbability(0.5);
    s.SetCrossoverMode(UNIFORM);                  //{UNIFORM, TWO_POINT}
    s.SetMaxExpectedOffspring(1.1);
    s.SetElitistFraction(0.02);
    s.SetSearchConstraint(1);
    s.SetReEvaluationFlag(0);
  // redirect standard output to a file

     
  #ifdef PRINTTOFILE
      ofstream evolfile;
      evolfile.open (supArgs1.rename_file("fitness.dat"));
      

      std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
      cout.rdbuf(evolfile.rdbuf());
  #endif
    // Code to run simulation:
    
    s.SetEvaluationFunction(EvaluationFunction);
    s.ExecuteSearch();

    #ifdef PRINTTOFILE
        std::cout.rdbuf(coutbuf); //reset to standard output again
        evolfile.close();
    #endif

    cout << "Finished, now saving the best fit...\n";
    
    
    RandomState rs;
    long seed = static_cast<long>(time(NULL));
    rs.SetRandomSeed(supArgs1.randomseed);
    ifstream Best;
    Best.open(supArgs1.rename_file("best.gen.dat"));
    TVector<double> best(1, VectSize);
    Best >> best;
    
    TVector<double> phenotype(1, VectSize);
    GenPhenMapping(best, phenotype);
    double sra = phenotype(SR_A);
    double srb = phenotype(SR_B);
    
    Worm w(phenotype, 1);
    {
    ofstream phenfile(supArgs1.rename_file("phenotype.dat"));
    w.DumpParams(phenfile);
    phenfile.close();
    }
    
    w.InitializeState(rs);
    w.sr.SR_A_gain = 0.0;
    w.sr.SR_B_gain = srb;
    w.AVA_output =  w.AVA_inact;
    w.AVB_output =  w.AVB_act;

    writeParsToJson(w);

    }

    supArgs1.setSimRandomInit();
    RandomState rs;
    //long seed = static_cast<long>(time(NULL));
    rs.SetRandomSeed(supArgs1.randomseed);
    ifstream Best;
    Best.open(supArgs1.rename_file("best.gen.dat"));
    TVector<double> best(1, VectSize);
    Best >> best;
    
    if (supArgs1.do_nml){
      supArgs1.doOrigNS = 0;
    cout << "Performing nml run and saving data\n" << endl;
    }
    else{
      supArgs1.doOrigNS = 1;
    cout << "Performing C++ run and saving data\n" << endl;
    }
    save_traces(best, rs);

    
    return 0;
}
