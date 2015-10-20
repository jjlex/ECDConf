/**********************************************************************
2015-10-16
ECDConf : Conformational search based OpenBabel and GALib libraries
***********************************************************************/

// ---- For Customized code -------
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>

// ---- For OpenBabel Lib -----
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/parsmart.h>
#include <openbabel/rotamer.h>
//#include <unistd.h>
#include <openbabel/obconversion.h>
#include <openbabel/base.h>
#include <openbabel/forcefield.h>

using namespace std;
using namespace OpenBabel;
//
// ----  For GALib -----------
#include <stdio.h>
#include <math.h>

#include <ga/GASimpleGA.h>
#include <ga/GABin2DecGenome.h>
#include <ga/std_stream.h>


//#define cout STD_COUT
// ------------------------------------------------------------
//  customized class 
// ------------------------------------------------------------
class ECDtor {
  public:
    int a1,a2,a3,a4 ;  // the atomic id of a torsional bond. a2-a3 is intersection.
    float torsion_min, torsion_max, torsion_step;
    double torsion_value;
    ECDtor () ;
    ECDtor ( int,int,int,int,float,float,float,double );
    void ini ();
};
// constructor:
  ECDtor::ECDtor () {
    a1=0 ; a2=0; a3=0; a4=0;
    torsion_min=0; torsion_max=360; torsion_step=0;
    torsion_value=0;
  }
  ECDtor::ECDtor (int b1, int b2, int b3, int b4, float t1, float t2, float t3, double t4) {
    a1=b1; a2=b2; a3=b3; a4=b4;
    torsion_min=t1; torsion_max=t2; torsion_step=t3;
    torsion_value=t4;
  }
  void ECDtor::ini () {
    a1=0 ; a2=0; a3=0; a4=0;
    torsion_min=0; torsion_max=360; torsion_step=0;
    torsion_value=0;
  }  
//
class ctrlvars {
  public:
    string molIn, molOut;
    unsigned int seed ;
    int popsize, ngen, nBestGenomes ;
    float pmut,pcross ;
    ctrlvars ();
};
// constructor:
  ctrlvars::ctrlvars () {
    molIn=" ";
    molOut=" ";
    seed=0;
    popsize= 50;
    ngen = 1000 ;
    nBestGenomes = 1;
    pmut   = 0.01;
    pcross = 0.6;
  }
    

// ------------------------------------------------------------
//  Gloable objects 
// ------------------------------------------------------------
float Objective(GAGenome &);
int gasimple();
double molEnergy();
void dumpGeom() ;
// ----------------------------
  OBMol mol; 
  OBForceField* pFF;
  OBConversion conv; //NF...
  ctrlvars eInp;
  ofstream ofs;
  vector<ECDtor> whiteTorsion(1);
  vector<ECDtor> blockTorsion(1);

// ------------------------------------------------------------
//  Custom Function
// ------------------------------------------------------------
 // read input: - Molecular filename
 //             - torsion list with user-defined angle range
 //             - blocked torsion list for dihedrals that should be defined as unrotable
 void readInput (ifstream & ifs) {
  string line,stmp;
  bool read_whitelist=false, read_blocklist=false;
  int i=0 ;
  ECDtor d1;
  whiteTorsion.clear();
  blockTorsion.clear();
  while( getline(ifs,line) ) {
     // 1). A: remove comment inputs
     if (line.find("#") < line.size()) line.erase(line.find("#"));
     //     B: remove all "," in line
     while ( line.find(",") < line.size() ) line.erase(line.find(","),1);
     //     C: replace all "=" with space
     while ( line.find("=") < line.size() ) line.replace(line.find("="),1," ");
     //     C: remove empty lines
     if (line.size() == 0 ) continue;
     // 2). setup section flag
     i++ ;
     if (line.find("whitelist") < line.size() ) {read_whitelist=true; read_blocklist=false ; i=0;}
     if (line.find("blocklist") < line.size() ) {read_whitelist=false; read_blocklist=true ; i=0;}

     // 3). read control parameters
     if ( line.find("molIn") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.molIn ;
     } else if ( line.find("molOut") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.molOut ;
     } else if ( line.find("seed") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.seed ;
     } else if ( line.find("popsize") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.popsize ;
     } else if ( line.find("ngen") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.ngen ;
     } else if ( line.find("nbestgenomes") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.nBestGenomes ;
     } else if ( line.find("pmut") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.pmut ;
     } else if ( line.find("pcross") < line.size() ) {
         istringstream iss(line);
         iss >> stmp >> eInp.pcross ;
     }
  
     // 4) A: read whitelist
     if ( i>0 && read_whitelist ) {
         d1.ini() ;
         istringstream iss(line);
         iss >> d1.a1 >> d1.a2 >> d1.a3 >> d1.a4;
         if ( iss.good() ) iss >> d1.torsion_min >> d1.torsion_max;
         if ( iss.good() ) iss >> d1.torsion_step ;
         whiteTorsion.push_back(d1);
     }
     //    B: read blocklist
     if ( i>0 && read_blocklist ) {
         d1.ini() ;
         istringstream iss(line);
         iss >> d1.a1 >> d1.a2 >> d1.a3 >> d1.a4;
         if ( iss.good() ) iss >> d1.torsion_min >> d1.torsion_max;
         if ( iss.good() ) iss >> d1.torsion_step ;
         blockTorsion.push_back(d1);
     }
  } // end of reading file (while loop)
 };

// ---------------------------------------------------------------
//  dumpGeom: output geometry based on torsionals
// ---------------------------------------------------------------
void dumpGeom() {
  OBAtom *a1, *a2, *a3, *a4;
  int n;
  double angle, energy;
  n=whiteTorsion.size() ;
  vector<ECDtor>::iterator it = whiteTorsion.begin();
  for(it=whiteTorsion.begin();it!=whiteTorsion.end();it++) {
    a1 = mol.GetAtom((*it).a1) ;
    a2 = mol.GetAtom((*it).a2) ;
    a3 = mol.GetAtom((*it).a3) ;
    a4 = mol.GetAtom((*it).a4) ;
    angle = (*it).torsion_value ;
    mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);
  }
  conv.Write(&mol,&ofs); //NF
}

///////////////////////////////////////////////////////////////////////////////
int main(int argc,char **argv)
{
  char *FileIn =NULL;
  char FileOut[200];

// 1). read command line for control parameters
   if (argc != 2) {
       cerr << "Usage: ecdconf <filename> " << endl;
       exit(-1);
   }
   FileIn = argv[1];
//OBPlugin::List("forcefields");
// 2).
   ifstream ifs;
   ifs.open(FileIn);
   readInput(ifs); // stmp carries filename for molecule
   ifs.close();
   strcpy(FileIn,eInp.molIn.c_str());
   strcpy(FileOut,eInp.molOut.c_str());
cout << "here" << endl;

  OBFormat* format = conv.FormatFromExt(FileIn);
cout << "here 1" << endl;
  if(!(format && conv.SetInAndOutFormats(format, format))) { //in and out formats same
    cerr << "ecdconf: cannot read and/or write this file format!" << endl;
    exit (-1);
  } //...NF

  // Read the file
  ifs.open(FileIn);
  ofs.open(FileOut);
  if (!ifs) {
    cerr << "ecdconf: cannot read input file!" << endl;
    exit (-1);
  }

  for (;;) {
    mol.Clear();
    //NF      ifs >> mol;                   // Read molecule
    conv.Read(&mol,&ifs); //NF
    if (mol.Empty())
      break;

// -- newly added
//    string ff="mmff94.ff" ;
    string ff="" ;
    pFF = OBForceField::FindForceField(ff);
    if (!pFF) {
      cerr << " could not find forcefield '" << ff << "'." <<endl;
      exit (-1);
    }
    pFF->SetLogFile(&cout);
//    pFF->SetLogLevel(OBFF_LOGLVL_MEDIUM);

  cout << "initial dihedral values:\n" ;
  OBAtom *a1, *a2, *a3, *a4;
  vector<ECDtor>::iterator it ;
  double angle;
  for(it=whiteTorsion.begin();it!=whiteTorsion.end();it++) {
    a1 = mol.GetAtom((*it).a1) ;
    a2 = mol.GetAtom((*it).a2) ;
    a3 = mol.GetAtom((*it).a3) ;
    a4 = mol.GetAtom((*it).a4) ;
    (*it).torsion_value = mol.GetTorsion(a1, a2, a3, a4);
    angle =  (*it).torsion_value ;
    if ( angle < 0 ) angle += 360;
    cout.width(10); 
    cout << angle << " " ;
  }
  cout << "\n" ;

    gasimple();

//----
    //NF      cout << mol;
//    conv.Write(&mol,&cout); //NF
  }

  return(0);
}
/* ----------------------------------------------------------------------------
   GALib code
---------------------------------------------------------------------------- */
int  gasimple()  {
// Declare variables for the GA parameters and set them to some default values.
  unsigned int seed = eInp.seed;
/*
  int popsize  = 50;
  int ngen     = 1000;
  int nBestGenomes = 5;
  float pmut   = 0.01;
  float pcross = 0.6;
*/

// Create a phenotype then fill it with the torsionals which will be mapped to binary
// n : number of torsionals
// The arguments to the add() method of a  Bin2Dec phenotype are :
//   - (1) number of bits,
//   - (2) min value, 
//   - (3) max value.

  GARandomSeed(seed);
  unsigned int n=whiteTorsion.size();
  float min[n], max[n];
  unsigned int nbits=8; // number of bits for each gene
  unsigned int i;

  vector<ECDtor>::iterator it = whiteTorsion.begin() ;
  for(it=whiteTorsion.begin(), i=0;i<n;i++, it++) {
      min[i] = (*it).torsion_min ;
      max[i] = (*it).torsion_max ;
  }

//  unsigned int nbits=16; 
  GABin2DecPhenotype map;
  for(i=0; i<n; i++)
    map.add(nbits, min[i], max[i]);

// Create the template genome using the phenotype map we just made.  The
// GA will use this genome to clone the population that it uses to do the
// evolution.  We pass the objective function to create the genome. 

//  GABin2DecGenome genome(map, Objective, (void *)target);
  GABin2DecGenome genome(map, Objective );

// Now create the GA using the genome, set the parameters, and run it.

  GASimpleGA ga(genome);
  ga.populationSize(eInp.popsize);
  ga.nGenerations(eInp.ngen);
  ga.pMutation(eInp.pmut);
  ga.pCrossover(eInp.pcross);
  ga.nBestGenomes(eInp.nBestGenomes);
//  ga.scoreFilename("bog.dat");
//  ga.flushFrequency(50);	// dump scores to disk every 50th generation
  ga.evolve(seed);

// Dump the results of the GA to the screen.  We print out first what a random
// genome looks like (so we get a bit of a feel for how hard it is for the
// GA to find the right values) then we print out the best genome that the
// GA was able to find.

  genome.initialize();
/*
  cout << "random values in the genome:\n";;
  for(i=0; i<map.nPhenotypes(); i++){
    cout.width(10); cout << genome.phenotype(i) << " ";
  }
  cout << "\n";
*/
//------
  GAPopulation gapop = ga.statistics().bestPopulation();
  n=gapop.size() ;
  gapop.sort();
  int j;
  for(i=0;i<n;i++) {
     genome = gapop.individual(i) ;
     for(j=0, it=whiteTorsion.begin(); j<genome.nPhenotypes(); j++, it++) {
       (*it).torsion_value = genome.phenotype(j) ;
     }
//     dumpGeom() ;
     cout << "The energy for " << i << " is: " << molEnergy() << endl;
     conv.Write(&mol,&ofs); //NF
  }

//------

  genome = ga.statistics().bestIndividual();
  cout << "the ga generated:\n";
  for(i=0; i<map.nPhenotypes(); i++){
    cout.width(10); cout << genome.phenotype(i) << " ";
  }
  cout << "\n\n"; cout.flush();

// We could print out the genome directly, like this:
// cout << genome << "\n";

  return 0;
}
 

// For this objective function we try to match the values in the array of float
// that is passed to us as userData.  If the values in the genome map to 
// values that are close, we return a better score.  We are limited to positive
// values for the objective value (because we're using linear scaling - the
// default scaling method for SimpleGA), so we take the reciprocal of the 
// absolute value of the difference between the value from the phenotype and 
// the value in the sequence.
double molEnergy() {
  OBAtom *a1, *a2, *a3, *a4;
  OBMol moltmp;
  int n;
  double angle, energy;
  n=whiteTorsion.size() ;
  vector<ECDtor>::iterator it = whiteTorsion.begin();
  for(it=whiteTorsion.begin();it!=whiteTorsion.end();it++) {
    a1 = mol.GetAtom((*it).a1) ;
    a2 = mol.GetAtom((*it).a2) ;
    a3 = mol.GetAtom((*it).a3) ;
    a4 = mol.GetAtom((*it).a4) ;
    angle = (*it).torsion_value ;
    mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);
  }

  moltmp = mol ;
  moltmp.DeleteHydrogens();
  if (!pFF->Setup(moltmp)) {
    cerr << "could not setup force field." << endl;
    exit (-1);
  }
  energy = pFF->Energy(false);
/*
  energy = pFF->E_Torsion(false) + pFF->E_VDW(false);
  cout << "energy is:" << energy << endl ;
  cout << "energy is:" << pFF->E_Bond(false) << endl ;
  cout << "energy is:" << pFF->E_Angle(false) << endl ;
*/
   return energy;
}
//
///////////////////////////////////////////////
// Objective function for GA
float
Objective(GAGenome& g)
{
  GABin2DecGenome & genome = (GABin2DecGenome &)g;
  int i;
//  float *sequence = (float *)g.userData();

  float value=genome.nPhenotypes();
  double energy;
  vector<ECDtor>::iterator it = whiteTorsion.begin();
  for( i=0, it=whiteTorsion.begin(); i<genome.nPhenotypes(); i++, it++) {
    (*it).torsion_value = genome.phenotype(i) ;
  }
  value = 500.0 - molEnergy() ;
/*
  energy = molEnergy() ;
  value = -energy/fabs(energy) * log(1+fabs(energy)) + 10 ;
*/
  if ( value < 0 ) value=0;

  return value;
}
