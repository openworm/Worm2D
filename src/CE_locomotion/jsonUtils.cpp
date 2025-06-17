#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <functional>
#include <vector>
#include <nlohmann/json.hpp>
//#include "utils.h"
#include "jsonUtils.h"
#include "argUtils.h"
#include "Mainvars.h"

extern SuppliedArgs supArgs1;
using json = nlohmann::json;


using std::istringstream;
using std::string;
using std::cout;
using std::function;
using std::vector;



//// write all worm paramters to json and txt file 

void writeWormParams(wormForJson & w)
{
    writeParsToJson(w);

    {
    ofstream nv_file(supArgs1.rename_file("w_verb.dat"));
    writeNSysToFile(nv_file, static_cast<NervousSystem&>(*w.n_ptr));
    nv_file << endl;
    writeWSysToFile(nv_file, w);
    nv_file << endl;
    writeMuscSysToFile(nv_file, w.m);
    nv_file << endl;
    writeGlobalParsToFile(nv_file);
    nv_file << endl;
    writeStretchSysToFile(nv_file, w.sr);
    nv_file << endl;
    writeBodySysToFile(nv_file, w.b);
    nv_file.close();
    }
   
   
}


//// Params structure




Params<int> getNervousSysParamsIntNH(NervousSystem& c)
{
Params<int> par;    
par.names = {"size", "maxchemcons", "maxelecconns"};
par.vals = {c.size, c.maxchemconns, c.maxelecconns};
return par;
}

ParamsHead<int> getNervousSysParamsInt(NervousSystem& c)
{
Params<int> par = getNervousSysParamsIntNH(c);    
ParamsHead<int> parH("Nervous system", par);
return parH;
}
      

Params< vector<int> > getNervousSysCellGroups(NervousSystem& c)
{
Params< vector<int> > par;
par.names = {"Cell group"};
vector<int> group_names; 
for (int i=0;i<10;i++){
vector<int> vec1;
for (int j=0;j<6;j++) vec1.push_back(i);
group_names.insert(group_names.end(),vec1.begin(),vec1.end());
}
par.vals = {group_names};
return par;
}




Params< vector<double> > getNervousSysParamsDoubleNH(NervousSystem& c)
{
Params< vector<double> > par;
par.names = {"taus", "biases", "gains", "outputs", "states", "paststates", "Rtaus", "externalinputs"};
par.vals = {
getVector<double>(c.taus, c.size), 
getVector<double>(c.biases, c.size), 
getVector<double>(c.gains, c.size),
getVector<double>(c.outputs, c.size),
getVector<double>(c.states, c.size),
getVector<double>(c.paststates, c.size),
getVector<double>(c.Rtaus, c.size),
getVector<double>(c.externalinputs, c.size),
};
return par;
}

Params< vector<int> > getNervousSysVecInt(NervousSystem& c)
{
Params< vector<int> > par;
par.names = {"NumChemicalConns", "NumElectricalConns"};
par.vals = {
getVector<int>(c.NumChemicalConns, c.size), 
getVector<int>(c.NumElectricalConns, c.size), 
};
return par;
}


Params< vector<string> > getNervousSysCellNames(NervousSystem& c)
{
Params< vector<string> > par;
par.names = {"Cell name"};
vector<string> cell_names_all;
vector<string> cell_names = {"DA", "DB", "DD", "VD", "VA", "VB"};
for (int i=0;i<10;i++) cell_names_all.insert(cell_names_all.end(),cell_names.begin(),cell_names.end());
par.vals = {cell_names_all};
return par;
}


ParamsHead< vector<double> > getNervousSysParamsDouble(NervousSystem& c)
{
Params< vector<double> > par = getNervousSysParamsDoubleNH(c);
ParamsHead< vector<double> > parH("Nervous system",par);
return parH;
}


Params<double> getMusclesParamsDouble(Muscles & m)
{
Params<double> par;
par.names = {"T_muscle"};
par.vals = {m.T_muscle};
return par;
}

Params<int> getMusclesParamsInt(Muscles & m)
{
Params<int> par;
par.names = {"Nmuscles"};
par.vals = {m.Nmuscles};
return par;
}



Params<double> getWormParams(wormForJson & w)
{

Params<double> par;
par.names = {"NMJ_DA", "NMJ_DB", "NMJ_VD", "NMJ_VB", "NMJ_VA", "NMJ_DD"};
par.vals = {w.NMJ_DA, w.NMJ_DB, w.NMJ_VD, w.NMJ_VB, w.NMJ_VA, w.NMJ_DD};
append<string>(par.names,{"AVA_act", "AVA_inact", "AVB_act", "AVB_inact"});
append<string>(par.names,{"AVA_output", "AVB_output"});
append<double>(par.vals,{w.AVA_act, w.AVA_inact, w.AVB_act, w.AVB_inact});
append<double>(par.vals,{w.AVA_output, w.AVB_output});

return par;

}


vector<ParamsHead<double> > getGlobalParamsDouble()
{
vector<ParamsHead<double> > parvec;
{ParamsHead<double> par;
par.head = "Worm";
par.names = {"T_muscle"};
par.vals = {T_muscle};
par.messages ={"Muscle time constant"};
par.messages_inds = {0};
parvec.push_back(par);}
{ParamsHead<double> par;
par.head =  "Integration parameters";
par.names =  {"Transient", "StepSize"};
par.vals = {Transient, StepSize};
parvec.push_back(par);}
{ParamsHead<double> par;
par.head =  "Fitness traj";
par.names =  {"AvgSpeed", "BBCfit"};
par.vals = {AvgSpeed, BBCfit};
par.messages ={"Average speed of the worm in meters per second", "AvgSpeed*Duration"};
par.messages_inds = {0,1};
parvec.push_back(par);}
{ParamsHead<double> par;
par.head =  "Genotype -> Phenotype Mapping Ranges";
par.names =  {"BiasRange", 
     "SCRange", 
     "CSRange", 
     "ESRange",  
      "SRmax",  
      "NMJmax",  
      "NMJmin"};  
par.vals = {BiasRange, SCRange, CSRange, ESRange, SRmax, NMJmax, NMJmin}; 
parvec.push_back(par);}

return parvec;
}

vector<ParamsHead<int> > getGlobalParamsInt()
{
vector<ParamsHead<int> > parvec;
{ParamsHead<int> par;
par.head = "Worm";
par.names = { "N_muscles", "N_units", "N_neuronsperunit", "N_stretchrec", "NmusclePerNU"};
par.vals = {N_muscles, N_units, N_neuronsperunit, N_stretchrec, NmusclePerNU};
par.messages = {"Number of muscles on each side of the body (ventral and dorsal)", 
"Number of neural units", "Number of neurons in a neural unit",  
"Number of stretch receptors in each of dorsal A, ventral A, dorsal B and ventral B", 
"Number of muscles per neural unit on each side (dorsal and ventral)"};
par.messages_inds = {0,1,2,3,4};
parvec.push_back(par);}
{ParamsHead<int> par;
par.head = "Name conventions";
par.names =  {"DA","DB","DD","VD","VA","VB","Head","Tail"};
par.vals = {DA,DB,DD,VD,VA,VB,Head,Tail};
parvec.push_back(par);}
{ParamsHead<int> par;
par.head =  "Integration parameters";
par.names =  {"skip_steps", "N_curvs"};
par.vals = {skip_steps, N_curvs};
parvec.push_back(par);}
{ParamsHead<int> par;
par.head =  "Stretch receptor parameters";
par.names =  {"SR_A", "SR_B"};
par.vals = {SR_A, SR_B};
parvec.push_back(par);}
{ParamsHead<int> par;
par.head =  "Worm";
par.names =  {"VectSize"};
par.vals = {VectSize};
par.messages = {"Size of genotype"};
par.messages_inds = {0}; 
parvec.push_back(par);}

return parvec;
}


Params<double> getStretchReceptorParams(StretchReceptor& s)
{
Params<double> par;
par.names = {"NSR", "NSEGS", "NSEGSSR", "SR_A_gain", "SR_B_gain"};
par.vals = {s.NSR, s.NSEGS, s.NSEGSSR, s.SR_A_gain, s.SR_B_gain};
par.messages = {"Number of stretch receptor in DA, DB, VA and VB, equal to number of units",
                                "Number of segments in the body",
                                "Number of segments sensed by each stretch receptor"};
par.messages_inds = {0,1,2}; //must be ordered
return par;
}

struct stretchReceptorProjection{

    vector<toFromWeight> dorsalA, ventralA, dorsalB, ventralB;
};


struct VDProjection{

    vector<toFromWeight> dorsal, ventral;
};

VDProjection getMuscleToBodyProjection()
{

VDProjection bp;
const double weight = 0.5;

  {const int to = 1; 
  bp.dorsal.push_back(toFromWeight(weightentry{1,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{1,weight},to));}
  {const int to = 2; 
  bp.dorsal.push_back(toFromWeight(weightentry{1,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{1,weight},to));}

  //  All other segments receive force from two muscles
  for (int i = 3; i <= N_segments-2; i++)
  {
    const int to = i;
    const int mi = (int) ((i-1)/2);
    bp.dorsal.push_back(toFromWeight(weightentry{mi,weight},to));
    bp.ventral.push_back(toFromWeight(weightentry{mi,weight},to));
    bp.dorsal.push_back(toFromWeight(weightentry{mi+1,weight},to));
    bp.ventral.push_back(toFromWeight(weightentry{mi+1,weight},to));
  }

  //  Last two segments receive special treatment because they are only affected by a single muscle
  {const int to = N_segments-1; 
  bp.dorsal.push_back(toFromWeight(weightentry{N_muscles,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{N_muscles,weight},to));}
  {const int to = N_segments;
  bp.dorsal.push_back(toFromWeight(weightentry{N_muscles,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{N_muscles,weight},to));}

return bp;
}

VDProjection getHubToMuscleProjection()
{
VDProjection bp;

const double weight = 1.0;
// muscles 1 to 3
for (int mi=1; mi<=3; mi++){ 
    const int to = mi;
    const int from = 1;
    bp.dorsal.push_back(toFromWeight(weightentry{from,weight},to));
    bp.ventral.push_back(toFromWeight(weightentry{from,weight},to));
  }

  {const int to = 4; // 4th muscle
  bp.dorsal.push_back(toFromWeight(weightentry{1,weight},to));
  bp.dorsal.push_back(toFromWeight(weightentry{2,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{1,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{2,weight},to));}

  {const int to = 5; // 5th muscle
  bp.dorsal.push_back(toFromWeight(weightentry{2,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{2,weight},to));}
 
  int mt = 2; // Muscles 6-19
  for (int mi=6; mi<=19; mi++){
    const int to = mi;
    bp.ventral.push_back(toFromWeight(weightentry{mt,weight},to));
    bp.ventral.push_back(toFromWeight(weightentry{mt+1,weight},to));
    bp.dorsal.push_back(toFromWeight(weightentry{mt,weight},to));
    bp.dorsal.push_back(toFromWeight(weightentry{mt+1,weight},to));
    mt += mi%2; // increment the index for the innervating unit each two muscles, starting from mi = 7
  }

  {const int to = 20; // 20th muscle
  bp.dorsal.push_back(toFromWeight(weightentry{9,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{9,weight},to));}

  {const int to = 21; // 21th muscle
  bp.dorsal.push_back(toFromWeight(weightentry{9,weight},to));
  bp.dorsal.push_back(toFromWeight(weightentry{10,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{9,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{10,weight},to));}

  // Muscles 22-24
  for (int mi=22; mi<=24; mi++){
  const int to = mi; 
  bp.dorsal.push_back(toFromWeight(weightentry{10,weight},to));
  bp.ventral.push_back(toFromWeight(weightentry{10,weight},to));}

return bp;
}


//projection from 60 NS neurons to 10 muscle hubs, here included NMJs as weights

VDProjection getMuscleInputHubProjection(wormForJson & w){

VDProjection bp;

for (int i=1; i<=N_units; i++){

    const int to = i;
    const vector<int> dorsal_list = {DA,DB,DD};
    const vector<double> dorsal_list_NMJ = {w.NMJ_DA,w.NMJ_DB,w.NMJ_DD};
    for (size_t j=0;j<dorsal_list.size();j++) 
    bp.dorsal.push_back(toFromWeight(weightentry{nn(dorsal_list[j],i),dorsal_list_NMJ[j]},to));

    const vector<int> ventral_list = {VA,VB,VD};
    const vector<double> ventral_list_NMJ = {w.NMJ_VA,w.NMJ_VB,w.NMJ_VD};
    for (size_t j=0;j<ventral_list.size();j++) 
    bp.ventral.push_back(toFromWeight(weightentry{nn(ventral_list[j],i),ventral_list_NMJ[j]},to));    

    }
return bp;
}


//projection from 10 stretch receptors to 60 neurons

stretchReceptorProjection getStretchToNSProjection()
{
stretchReceptorProjection bp;
const double weight = 1.0;

for (int i = 1; i <= N_units; i++){
   
    const int from = i;
    {const int to = nn(DA,i);
    bp.dorsalA.push_back(toFromWeight(weightentry{from,weight},to));}
    {const int to = nn(VA,i);
    bp.ventralA.push_back(toFromWeight(weightentry{from,weight},to));}
    {const int to = nn(DB,i);
    bp.dorsalB.push_back(toFromWeight(weightentry{from,weight},to));}
    {const int to = nn(VB,i);
    bp.ventralB.push_back(toFromWeight(weightentry{from,weight},to));}

}
return bp;
}

//projection from 50 body segments to 10 stretch receptors

stretchReceptorProjection getBodyToStretchReceptorProjection(StretchReceptor& s)
{
stretchReceptorProjection bp;    

const double weight = 1.0;
{const int to = 1;
for (int j = 1; j <= s.NSEGSSR; j++){
const int from = j;
bp.dorsalA.push_back(toFromWeight(weightentry{from,weight},to));
bp.ventralA.push_back(toFromWeight(weightentry{from,weight},to));
}}

for (int to = 2; to <= 10; to++){
for (int j = 1; j <= s.NSEGSSR; j++)
        {
            const int from = j+(to-2)*4;
            bp.dorsalA.push_back(toFromWeight(weightentry{from,weight},to));
            bp.ventralA.push_back(toFromWeight(weightentry{from,weight},to));
        }
}

for (int to = 1; to <= 9; to++){
        for (int j = 1; j <= s.NSEGSSR; j++)
        {
            const int from = 12+j+(to-1)*4;
            bp.dorsalB.push_back(toFromWeight(weightentry{from,weight},to));
            bp.ventralB.push_back(toFromWeight(weightentry{from,weight},to));
        }
}

{const int to = 10;
for (int j = 1; j <= s.NSEGSSR; j++){
const int from = j+44;
bp.dorsalB.push_back(toFromWeight(weightentry{from,weight},to));
bp.ventralB.push_back(toFromWeight(weightentry{from,weight},to));
}}

return bp;
}



Params<int> getBodyParamsInts(WormBody& b)
{
Params<int> par;
par.names = {"N_segments"};
par.vals = {N_segments};
par.messages = {"Number of body segments on each side, dorsal and ventral"};
par.messages_inds  = {0};
return par;
}

Params<double> getBodyParams(WormBody& b)
{

Params<double> par;

par.names = {"Medium", "L_worm", "R_min", "C_agar_par_total", 
"C_agar_perp_total", "C_water_par_total", "C_water_perp_total", "kappa_L", "kappa_D", 
"kappa_M0", "beta_L", "beta_D", "beta_M0", "delta_M"};

par.vals = {Medium, L_worm, R_min, C_agar_par_total, C_agar_perp_total, 
C_water_par_total, C_water_perp_total, kappa_L, kappa_D, kappa_M0, beta_L, 
beta_D, beta_M0, delta_M};

par.messages_inds.resize(par.vals.size());

for (size_t i=0;i<par.messages_inds.size();i++) par.messages_inds[i]=i;

par.messages = {    
"Normalized medium drag coefficient (0 = water, 1 = agar)",
"Length of worm in m",
"Minor radius of prolate ellipse body in m",
"Total tangential drag coefficient for agar in kg/s",
"Total rod normal drag coefficient in agar in kg/s",
"Total rod tangential drag coefficient for water in kg/s",
"Total rod normal drag coefficient for water in kg/s",
"Lateral spring constant in kg/s",
"Diagonal spring constant in kg/s",
"Baseline active muscle spring constant in kg/s",
"Lateral passive damping constant in s",
"Diagonal passive damping constant in s",
"Baseline active damping constant in s",
"Rest muscle length scaling constant"

};

return par;

}

//////// functions for reading and writing json



void to_json(json & j, const weightentry & w)
{
  j = json{{"from", w.from}, {"weight", w.weight}};
}

void to_json(json & j, const toFromWeight & w)
{
  j = json{{"to", w.to}, {"from", w.w.from}, {"weight", w.w.weight}};
}

void from_json(const json& j, toFromWeight & w) 
{
        j.at("to").get_to(w.to);
        j.at("from").get_to(w.w.from);
        j.at("weight").get_to(w.w.weight);
}

template<class T>
void appendToJson(json & j, const Params<T> & par)
{
    size_t mess_ind = 0;
    for (size_t i=0;i<par.names.size(); i++) {
        if (par.messages_inds.size()>mess_ind && par.messages_inds[mess_ind]==static_cast<int>(i)) 
        {j[par.names[i]]["message"] = par.messages[i];mess_ind++;}
        j[par.names[i]]["value"] = par.vals[i];
        }
               
}

void appendMatrixToJson(json & j, TMatrix<weightentry> & vec, TVector<int> & sizes, int tot_size)
{    
    vector<toFromWeight> newvec;
    for (int i=1; i<=tot_size; i++){    
        for (int j=1; j<=sizes[i]; j++) { 
            toFromWeight tv(vec[i][j], i);
            newvec.push_back(tv);}        
    }
    j["value"] = newvec;

}

void appendMuscleToBodyProjection(json & j){

VDProjection bp = getMuscleToBodyProjection();
j["dorsal"]["value"] = bp.dorsal;
j["ventral"]["value"] = bp.ventral;
j["dorsal"]["message"] = "Projection from dorsal muscles to dorsal body segments";
j["ventral"]["message"] = "Projection from ventral muscles to ventral body segments";

}


void appendHubToMuscleProjection(json & j){

VDProjection bp = getHubToMuscleProjection();
j["dorsal"]["value"] = bp.dorsal;
j["ventral"]["value"] = bp.ventral;
j["dorsal"]["message"] = "Projection from dorsal hubs (10) to dorsal muscles (24)";
j["ventral"]["message"] = "Projection from ventral hubs (10) to ventral muscles (24)";

}

void appendMuscleInputHubProjection(json & j, wormForJson & w){

VDProjection bp = getMuscleInputHubProjection(w);
j["dorsal"]["value"] = bp.dorsal;
j["ventral"]["value"] = bp.ventral;
j["dorsal"]["message"] = "Projection from dorsal motorneurons to dorsal hubs with NMJ weights";
j["ventral"]["message"] = "Projection from ventral motorneurons to ventral hubs with NMJ weights";
}


void appendBodyStretchProjToJson(json & j, StretchReceptor& s)
{
stretchReceptorProjection bp = getBodyToStretchReceptorProjection(s);
j["dorsalA"]["value"] = bp.dorsalA;
j["dorsalB"]["value"] = bp.dorsalB;
j["ventralA"]["value"] = bp.ventralA;
j["ventralB"]["value"] = bp.ventralB;
j["dorsalA"]["message"] = "Projection from body segs to dorsalA (DA) stretch receptors";
j["dorsalB"]["message"] = "Projection from body to dorsalB (DB) stretch receptors";
j["ventralA"]["message"] = "Projection from body to ventralA (VA) stretch receptors";
j["ventralB"]["message"] = "Projection from body to ventralB (VB) stretch receptors";
}

void appendStretchToNSProjToJson(json & j)
{
stretchReceptorProjection bp = getStretchToNSProjection();
j["dorsalA"]["value"] = bp.dorsalA;
j["dorsalB"]["value"] = bp.dorsalB;
j["ventralA"]["value"] = bp.ventralA;
j["ventralB"]["value"] = bp.ventralB;
j["dorsalA"]["message"] = "Projection from SR to DA neurons";
j["dorsalB"]["message"] = "Projection from SR to DB neurons";
j["ventralA"]["message"] = "Projection from SR to VA neurons";
j["ventralB"]["message"] = "Projection from SR to VB neurons";
}


void appendNSToJson(json & j, NervousSystem& c)
{
    j["Chemical weights"]["message"] = "chemical weights in sparse format";
    appendMatrixToJson(j["Chemical weights"], c.chemicalweights, c.NumChemicalConns, c.size);
    appendMatrixToJson(j["Electrical weights"], c.electricalweights, c.NumElectricalConns, c.size);
    j["Electrical weights"]["message"] = "electrical weights in sparse format";
}

json getNSJson(NervousSystem & n)
{
json j;
string nsHead = "Nervous system";

{Params<vector<double> > parvec = getNervousSysParamsDoubleNH(n);
appendToJson<vector<double> >(j[nsHead],parvec);}

{Params<int> parvec = getNervousSysParamsIntNH(n);
appendToJson<int>(j[nsHead],parvec);}

{Params< vector<string> > parvec = getNervousSysCellNames(n);
appendToJson<vector<string> >(j[nsHead],parvec);}

appendNSToJson(j[nsHead], n);

return j;
}

ofstream & writeNSJson(ofstream & ofs, NervousSystem & n)
{
json j = getNSJson(n);
ofs << std::setw(4) << j << std::endl;
return ofs;
}

json getJsonFromFile(ifstream & ifs)
{
    //json j; 
    return json::parse(ifs);
    //ifs >> j;
    //return j;

}


void setNSFromJsonFile(ifstream & ifs, NervousSystem & n)
{

json j = getJsonFromFile(ifs);
json jns = j["Nervous system"];
auto outputs = jns["outputs"]["value"].template get< vector<double> >();
auto biases = jns["biases"]["value"].template get< vector<double> >();
auto chem_weights = jns["Chemical weights"]["value"].template get< vector<toFromWeight> >();
auto elec_weights = jns["Electrical weights"]["value"].template get< vector<toFromWeight> >();
auto gains = jns["gains"]["value"].template get< vector<double> >();
auto time_consts = jns["taus"]["value"].template get< vector<double> >();
auto maxchemcons = jns["maxchemcons"]["value"].template get<int>();
auto maxelecconns = jns["maxelecconns"]["value"].template get<int>();
auto size = jns["size"]["value"].template get<int>(); 
json jw = j["Worm"];

auto N_units_val = jw["N_units"]["value"].template get<int>();
auto N_neuronsperunit_val = jw["N_neuronsperunit"]["value"].template get<int>();

n.SetCircuitSize(N_units_val*N_neuronsperunit_val, 3, 2);
for (size_t i=0;i<biases.size();i++) n.SetNeuronBias(i+1, biases[i]);
for (size_t i=0;i<time_consts.size();i++) n.SetNeuronTimeConstant(i+1, time_consts[i]);
for (size_t i=0;i<outputs.size();i++) n.SetNeuronOutput(i+1, outputs[i]);


for (size_t i=0;i<chem_weights.size();i++)
n.SetChemicalSynapseWeight(chem_weights[i].w.from, chem_weights[i].to, chem_weights[i].w.weight);
for (size_t i=0;i<elec_weights.size();i++)
n.InternalSetElectricalSynapseWeight(elec_weights[i].w.from, elec_weights[i].to, elec_weights[i].w.weight);

//return n;
}



void writeParsToJson(json & j, wormForJson & w, string file_name)
{

{Params<double> par = getBodyParams(w.b);
appendToJson<double>(j["Body"],par);}

{Params<int> par = getBodyParamsInts(w.b);
appendToJson<int>(j["Body"],par);}

{Params<double> par = getStretchReceptorParams(w.sr);
appendToJson<double>(j["Stretch receptor"],par);
}
{
Params<double> par = getWormParams(w);
appendToJson<double>(j["Worm"],par);
}
{
Params<double> par =getMusclesParamsDouble(w.m);
appendToJson<double>(j["Muscle"],par);
}
{
Params<int> par =getMusclesParamsInt(w.m);
appendToJson<int>(j["Muscle"],par);
}

{vector<ParamsHead<int> > parvec = getGlobalParamsInt();
for (size_t i=0;i<parvec.size(); i++) {
appendToJson<int>(j[parvec[i].head],parvec[i]);
}}
{vector<ParamsHead<double> > parvec = getGlobalParamsDouble();
for (size_t i=0;i<parvec.size(); i++) {
appendToJson<double>(j[parvec[i].head],parvec[i]);
}}

/* {ParamsHead<vector<double> > parvec = getNervousSysParamsDouble(static_cast<NervousSystem&>(*w.n_ptr));
appendToJson<vector<double> >(j[parvec.head],parvec);}

ParamsHead<int> parvec = getNervousSysParamsInt(static_cast<NervousSystem&>(*w.n_ptr));
appendToJson<int>(j[parvec.head],parvec); */

string nsHead = "Nervous system";

{Params<vector<double> > parvec = getNervousSysParamsDoubleNH(static_cast<NervousSystem&>(*w.n_ptr));
appendToJson<vector<double> >(j[nsHead],parvec);}

{Params<int> parvec = getNervousSysParamsIntNH(static_cast<NervousSystem&>(*w.n_ptr));
appendToJson<int>(j[nsHead],parvec);}

{Params< vector<string> > parvec = getNervousSysCellNames(static_cast<NervousSystem&>(*w.n_ptr));
appendToJson<vector<string> >(j[nsHead],parvec);} 

{Params< vector<int> > parvec = getNervousSysVecInt(static_cast<NervousSystem&>(*w.n_ptr));
appendToJson<vector<int> >(j[nsHead],parvec);}

//{Params< vector<int> > parvec = getNervousSysCellGroups(static_cast<NervousSystem&>(*w.n_ptr));
//appendToJson<vector<int> >(j[nsHead],parvec);} 

appendNSToJson(j[nsHead], static_cast<NervousSystem&>(*w.n_ptr));

appendStretchToNSProjToJson(j[nsHead]);

appendBodyStretchProjToJson(j["Stretch receptor"], w.sr);

appendMuscleInputHubProjection(j["NSToMuscleHub"],w);
appendHubToMuscleProjection(j["HubToMuscle"]);
appendMuscleToBodyProjection(j["MuscleToBody"]);

ofstream json_out(supArgs1.rename_file(file_name));
json_out << std::setw(4) << j << std::endl;
json_out.close();
}

/* void writeParsToJson(wormForJson & w)
{
writeParsToJson(w, "worm_data.json");
} */


void writeParsToJson(wormForJson & w)
{
    if (supArgs1.doOrigNS){
  cout << "making json" << endl;
  vector<doubIntParamsHead> evolutionParams;
  doubIntParamsHead var1 = supArgs1.getParams();
  evolutionParams.push_back(var1);

  writeParsToJson(w, "worm_data.json", evolutionParams);
  //writeParsToJson(w, "worm_data.json");
  testNervousSystemJson("worm_data.json", static_cast<NervousSystem &>(*w.n_ptr)); 
    }
}

void writeParsToJson(wormForJson & w, string file_name)
{
json j;
writeParsToJson(j,w,file_name);
}

void writeParsToJson(wormForJson & w, string file_name, vector<doubIntParamsHead> & parvec)
{
json j;

for (size_t i=0;i<parvec.size(); i++) {
if (strcmp(parvec[i].parDoub.head.c_str(),"NULL")!=0)
appendToJson<double>(j[parvec[i].parDoub.head],parvec[i].parDoub);
if (strcmp(parvec[i].parInt.head.c_str(),"NULL")!=0)
appendToJson<long>(j[parvec[i].parInt.head],parvec[i].parInt);
}
writeParsToJson(j,w,file_name);
}




void readJson(json j, ifstream & ifs)
{
ifs >> j;
}

void testNervousSystemJson(string fileName, NervousSystem & n)
{
ifstream NS_ifs(supArgs1.rename_file(fileName));
setNSFromJsonFile(NS_ifs, n);
NS_ifs.close();
}

//////////// functions for writing text file



template<class T>
ostream& writeVectorFormat(ostream& os, 
const vector<string> & names, const vector<T> & vals, 
const vector<int> & messages_inds, const vector<string> & messages)
{
    size_t mess_ind = 0;
    os << setprecision(32);
    for (size_t i=0;i<names.size(); i++) {
        if (messages_inds.size()>mess_ind && messages_inds[mess_ind]==static_cast<int>(i)) {os << messages[mess_ind] << endl;mess_ind++;}
        os << names[i] + ": " << vals[i] << endl;
        }

return os;        
        
}


template<class T>
ostream& writeVectorFormat(ostream& os, 
const vector<string> & names, const vector<T> & vals)
{
    os << setprecision(32);
    for (size_t i=0;i<names.size(); i++) {
        os << names[i] + ": " << vals[i] << endl;
        }

return os;        
        
}    


ostream& writeBodySysToFile(ostream& os, WormBody& b)
{   

os << "Worm Body parameters" << endl;    
{Params<double> par = getBodyParams(b);
writeVectorFormat<double>(os,par.names,par.vals,par.messages_inds,par.messages);
}

{Params<int> par = getBodyParamsInts(b);
writeVectorFormat<int>(os,par.names,par.vals,par.messages_inds,par.messages);
}

return os;
}


ostream& writeStretchSysToFile(ostream& os, StretchReceptor& s)
{   
    
    os << setprecision(32);
    Params<double> par = getStretchReceptorParams(s);
    
 
    os << "StretchReceptor parameters" << endl;
    writeVectorFormat<double>(os,par.names,par.vals,par.messages_inds,par.messages);
    
return os;
}


ostream& writeGlobalParsToFile(ostream& os)
{   
    os << setprecision(32);
    os << "Worm" << endl;

    os << "N_muscles: " << N_muscles << endl;
    os << "N_units: " << N_units << endl;
    os << "N_neuronsperunit: " << N_neuronsperunit << endl;
    os << "N_stretchrec: " << N_stretchrec << endl;
    os << "T_muscle: " << T_muscle << endl;
    os << "NmusclePerNU: " << NmusclePerNU << endl << endl;

    os << "Name conventions" << endl;
    vector<string> names = {"DA","DB","DD","VD","VA","VB","Head","Tail"};
    vector<int> vals = {DA,DB,DD,VD,VA,VB,Head,Tail};
    for (size_t i=0; i<names.size(); i++)
    {
    os << names[i] + ": " << vals[i] << endl;
    }

    os << endl;
    os << "Integration parameters" << endl;
    os << "skip_steps: " << skip_steps << endl;
    os <<  "Transient: " << Transient << endl; 
    os <<  "StepSize: " << StepSize << endl;
    os <<  "N_curvs: " << N_curvs << endl;

    os << endl;
    os << "Fitness traj" << endl;
    os <<  "AvgSpeed: " << AvgSpeed << endl;
    os <<  "BBCfit: " <<  BBCfit << endl;

    os << endl;
    os <<  "Genotype -> Phenotype Mapping Ranges" << endl;
    os <<  "BiasRange: " << BiasRange << endl;
    os <<  "SCRange: " << SCRange << endl;
    os <<  "CSRange: " << CSRange  << endl;
    os <<  "ESRange: "  << ESRange << endl;
    os <<  "SRmax: "  << SRmax << endl;
    os <<  "NMJmax: "  << NMJmax << endl;
    os <<  "NMJmin: "  << NMJmin << endl;

    os << endl;
    os << "Stretch receptor parameters" << endl;
    os << "SR_A :" << SR_A << endl;
    os << "SR_B :" << SR_B << endl;
 
    os << endl;
    os << "Size of genotype" << endl;
    os << "Vectsize :" << VectSize << endl;

    
    return os;
}


ostream& writeWSysToFile(ostream& os, wormForJson & w)
{   
    
    os << setprecision(32);

    os << "Worm parameters" << endl;
    {
    vector<string> names = {"NMJ_DA", "NMJ_DB", "NMJ_VD", "NMJ_VB", "NMJ_VA", "NMJ_DD"};
    vector<double> vals = {w.NMJ_DA, w.NMJ_DB, w.NMJ_VD, w.NMJ_VB, w.NMJ_VA, w.NMJ_DD};

    for (size_t i=0; i<names.size(); i++)
    {
    os << names[i] + ": " << vals[i] << endl;
    }}
    os << endl;
    {
    vector<string> names = {"AVA_act", "AVA_inact", "AVB_act", "AVB_inact"};
    vector<double> vals = {w.AVA_act, w.AVA_inact, w.AVB_act, w.AVB_inact};

    for (size_t i=0; i<names.size(); i++)
    {
    os << names[i] + ": " << vals[i] << endl;
    }}
    os << endl;
    {
    vector<string> names = {"AVA_output", "AVB_output"};
    vector<double> vals = {w.AVA_output, w.AVB_output};

    for (size_t i=0; i<names.size(); i++)
    {
    os << names[i] + ": " << vals[i] << endl;
    }}
    return os;

}    

ostream& writeMuscSysToFile(ostream& os, Muscles& m)
{   
    os << setprecision(32);
    os << "T_muscle: " << m.T_muscle << endl;
    os << "Nmuscles: " << m.Nmuscles << endl;
    return os;
}


ostream& writeNSysToFile(ostream& os, NervousSystem& c)
{
    // Set the precision
    os << setprecision(32);
    os << "size: " << c.size << endl;
    os << "maxchemcons: " << c.maxchemconns << endl;
    os << "maxelecconns: " <<  c.maxelecconns << endl << endl;
    //os << "size, maxchemconns and maxelecconns" << endl;
    // Write the size, maxchemconns and maxelecconns
    //os << c.size << " " << c.maxchemconns << " " << c.maxelecconns << endl << endl;
    // Write the time constants
    os << "time constants: " << endl;
    for (int i = 1; i <= c.size; i++)
        os << c.taus[i] << " ";
    os << endl << endl;
    os << "biases: " << endl;
    // Write the biases
    for (int i = 1; i <= c.size; i++)
        os << c.biases[i] << " ";
    os << endl << endl;
    // Write the gains
    os << "gains: " << endl;
    for (int i = 1; i <= c.size; i++)
        os << c.gains[i] << " ";
    os << endl << endl;
    os << "chemical weights: " << endl;
    // Write the chemical weights in sparse format (N from1 weight1 ... fromN weightN)
    for (int i = 1; i <= c.size; i++) {
        //cout << c.NumChemicalConns[i] << "  ";
        os << c.NumChemicalConns[i] << "  "; //fixed
        for (int j = 1; j <= c.NumChemicalConns[i]; j++)
            os << c.chemicalweights[i][j].from << " " << c.chemicalweights[i][j].weight << "  ";
        os << endl;
    }
    os << endl;
    os << "electrical weights: " << endl;
    // Write the electrical weights in sparse format (N from1 weight1 ... fromN weightN)
    for (int i = 1; i <= c.size; i++) {
        //cout << c.NumElectricalConns[i] << "  ";
        os << c.NumElectricalConns[i] << "  "; //fixed 
        for (int j = 1; j <= c.NumElectricalConns[i]; j++)
            os << c.electricalweights[i][j].from << " " << c.electricalweights[i][j].weight << "  ";
        os << endl;
    }
    // Return the ostream
    return os;
}




istream& readNSysFromFile(istream& is, NervousSystem& c)
{
    // Read the sizes
    int size;
    is >> size;
    int maxchemconns;
    is >> maxchemconns;
    int maxelecconns;
    is >> maxelecconns;
    c.SetCircuitSize(size, maxchemconns, maxelecconns);
    // Read the time constants
    for (int i = 1; i <= size; i++) {
        is >> c.taus[i];
        c.Rtaus[i] = 1/c.taus[i];
    }
    // Read the biases
    for (int i = 1; i <= size; i++)
        is >> c.biases[i];
    // Read the gains
    for (int i = 1; i <= size; i++)
        is >> c.gains[i];
    // Read the chemical weights
    int n;
    for (int i = 1; i <= size; i++) {
        is >> n;
        for (int j = 1; j <= n; j++) {
            is >> c.chemicalweights[i][j].from;
            is >> c.chemicalweights[i][j].weight;
            c.NumChemicalConns[i]++;
        }
    }
    // Read the electrical weights
    for (int i = 1; i <= size; i++) {
        is >> n;
        for (int j = 1; j <= n; j++) {
            is >> c.electricalweights[i][j].from;
            is >> c.electricalweights[i][j].weight;
            c.NumElectricalConns[i]++;
        }
    }
    // Return the istream
    return is;
}