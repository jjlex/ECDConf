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
    torsion_min=-180; torsion_max=180; torsion_step=0;
    torsion_value=0;
  }
  ECDtor::ECDtor (int b1, int b2, int b3, int b4, float t1, float t2, float t3, double t4) {
    a1=b1; a2=b2; a3=b3; a4=b4;
    torsion_min=t1; torsion_max=t2; torsion_step=t3;
    torsion_value=t4;
  }
  void ECDtor::ini () {
    a1=0 ; a2=0; a3=0; a4=0;
    torsion_min=-180; torsion_max=180; torsion_step=0;
    torsion_value=0;
  }  
// ------------------------------------------------------------------------
class ECDdelBond {
  public:
    int a1,a2;
    ECDdelBond ();
    ECDdelBond ( int,int );
    void ini () ;
};
// constructor:
    ECDdelBond::ECDdelBond () {
      a1=0 ; a2=0 ;
    }
    ECDdelBond::ECDdelBond( int b1, int b2 ) {
      a1=b1; a2=b2 ;
    }
    void ECDdelBond::ini () {
      a1=0 ; a2=0;
    }
// ------------------------------------------------------------------------
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
void setGeom();
void DelBond ();
void RestoreBond() ;
void fixNonHydrogens();
// ----------------------------
  OBMol mol; 
  OBForceField* pFF;
  OBConversion conv; //NF...
  ctrlvars eInp;
  ofstream ofs;
  vector<ECDtor> RotorTorsion(1);
  vector<ECDtor> AddRotorTorsion(1);
  vector<ECDtor> RemoveRotorTorsion(1);
  vector<ECDdelBond> rDelBond(1);
  vector<int> VecDelbondID;
  bool doDeleteBond = false ;

// ------------------------------------------------------------
//  Custom Function
// ------------------------------------------------------------
 // read input: - Molecular filename
 //             - torsion list with user-defined angle range
 //             - blocked torsion list for dihedrals that should be defined as unrotable
 void readInput (ifstream & ifs) {
  string line,stmp;
  bool read_whitelist=false, read_blocklist=false, read_delBond=false;
  int i=0 ;
  ECDtor d1;
  AddRotorTorsion.clear();
  RemoveRotorTorsion.clear();
  ECDdelBond bond1;
  rDelBond.clear();
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
     if (line.find("whitelist") < line.size() ) { 
        read_whitelist=true; 
        read_blocklist=false ; 
        read_delBond=false;
        i=0;
      } else if (line.find("blocklist") < line.size() ) {
        read_whitelist=false; 
        read_blocklist=true ; 
        read_delBond=false;
        i=0;
      } else if (line.find("deletebond") < line.size() ) {
        read_whitelist=false; 
        read_blocklist=false ; 
        read_delBond=true;
        i=0;
      }

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
         AddRotorTorsion.push_back(d1);
     }
     //    B: read blocklist
     if ( i>0 && read_blocklist ) {
         d1.ini() ;
         istringstream iss(line);
         iss >> d1.a1 >> d1.a2 >> d1.a3 >> d1.a4;
         if ( iss.good() ) iss >> d1.torsion_min >> d1.torsion_max;
         if ( iss.good() ) iss >> d1.torsion_step ;
         RemoveRotorTorsion.push_back(d1);
     }
     //    C: read delBond list
     if ( i>0 && read_delBond ) {
         bond1.ini() ;
         istringstream iss(line);
         iss >> bond1.a1 >> bond1.a2 ;
         rDelBond.push_back(bond1);
     }
  } // end of reading file (while loop)
 };

// ---------------------------------------------------------------
// ---------------------------------------------------------------
void getAllTorsion () {
     int datoms[4];
     unsigned int i ;
     OBRotorList rl;
     OBRotorIterator ri;
     ECDtor d1;
     RotorTorsion.clear();

     rl.Setup(mol);
     OBRotor *rotor = rl.BeginRotor(ri);
     cout << " Number of rotatable bonds: " << rl.Size() << endl;
     for ( i = 1; i < rl.Size() + 1; ++i, rotor = rl.NextRotor(ri)) {
        d1.ini();
        rotor->GetDihedralAtoms(datoms) ;
        printf ("%5d %5d %5d %5d\n", datoms[0] ,datoms[1],datoms[2],datoms[3]);
        d1.a1 = datoms[0] ;
        d1.a2 = datoms[1] ;
        d1.a3 = datoms[2] ;
        d1.a4 = datoms[3] ;
        RotorTorsion.push_back(d1);
     }
}

// ---------------------------------------------------------------
void RestoreBond() {
  vector<ECDdelBond>::iterator it1;
  for(it1=rDelBond.begin(); it1!=rDelBond.end(); it1++) {
      mol.AddBond((*it1).a1,(*it1).a2,1);
  }
}

//;
void DelBond () {
  unsigned int a1,a2, b1,b2, bID;
  vector<ECDdelBond>::iterator it;
  for(it=rDelBond.begin(); it!=rDelBond.end(); it++) {
     a1 = (*it).a1 ; 
     a2 = (*it).a2;
     for( OBMolBondIter blist(&mol); blist; ++blist ) {
        b1 = blist->GetBeginAtomIdx() ;
        b2 = blist->GetEndAtomIdx() ;
        if ( (a1==b1&&a2==b2) || (a1==b2&&a2==b1) ) {
          bID = blist->GetIdx() ;
          mol.DeleteBond(mol.GetBond(bID));
        }
     }
  } // end of big for loop
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc,char **argv)
{
  char *FileIn =NULL;
  char FileOut[200];
  int i;

// 1). read command line for control parameters
   if (argc != 2) {
       cerr << "Usage: ecdconf <filename> " << endl;
       exit(-1);
   }
   FileIn = argv[1];

// 2).
   ifstream ifs;
   ifs.open(FileIn);
   readInput(ifs); // stmp carries filename for molecule
   ifs.close();
   strcpy(FileIn,eInp.molIn.c_str());
   strcpy(FileOut,eInp.molOut.c_str());
   if (rDelBond.size() > 0 ) doDeleteBond=true;

  OBFormat* format = conv.FormatFromExt(FileIn);
  if(!(format && conv.SetInAndOutFormats(format, format))) { //in and out formats same
    cerr << "obrotate: cannot read and/or write this file format!" << endl;
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
    mol.Center();
//--- this for comparison with initial 
// conv.Write(&mol,&ofs); //NF
    // Delte bonds in Ring(user defined);
    if ( doDeleteBond ) DelBond() ;
    getAllTorsion ();

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
  for(it=RotorTorsion.begin();it!=RotorTorsion.end();it++) {
cout << RotorTorsion[0].a1 << endl;
        printf ("%5d %5d %5d %5d\n", (*it).a1,(*it).a2,(*it).a3,(*it).a4);
    a1 = mol.GetAtom((*it).a1) ;
    a2 = mol.GetAtom((*it).a2) ;
    a3 = mol.GetAtom((*it).a3) ;
    a4 = mol.GetAtom((*it).a4) ;
    (*it).torsion_value = mol.GetTorsion(a1, a2, a3, a4);
    angle =  (*it).torsion_value ;
//    if ( angle < 0 ) angle += 360;
    cout.width(10); 
    cout << angle << " " ;
  }
  cout << "\n" ;

    gasimple();

//----
    //NF      cout << mol;
//    conv.Write(&mol,&cout); //NF
  } // end of "for loop" for molecules

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
  unsigned int n=RotorTorsion.size();
  float min[n], max[n];
  unsigned int nbits=8; // number of bits for each gene
  unsigned int i;

  vector<ECDtor>::iterator it = RotorTorsion.begin() ;
  for(it=RotorTorsion.begin(), i=0;i<n;i++, it++) {
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
  int steps = 100;
  double crit = 1e-3;
  bool done=true;
  double energy, energyLast=-1000, energyWindow;
  energyWindow = 1.0;
  // optimize hydrogens for the best population
  for(i=0;i<n;i++) {
     genome = gapop.individual(i) ;
     for(j=0, it=RotorTorsion.begin(); j<genome.nPhenotypes(); j++, it++) {
       (*it).torsion_value = genome.phenotype(j) ;
     }
     if (doDeleteBond) DelBond() ;
     setGeom() ;
     if (doDeleteBond) RestoreBond() ;
     // optimization
     pFF->Setup(mol);
     fixNonHydrogens();
     pFF->SteepestDescent(steps, crit);

     energy =  pFF->Energy(false);
     if ( energy > energyLast + energyWindow ) { 
        cout << "The energy for " << i << " is: " << energy << endl;
        energyLast = energy;
        pFF->GetCoordinates(mol);
        if (doDeleteBond) RestoreBond() ;
        mol.Center();
        conv.Write(&mol,&ofs); //NF
     }

  }

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
void fixNonHydrogens () {
  unsigned int i,NAtoms;
  OBForceField:OBFFConstraints Constraints;
  NAtoms = mol.NumAtoms();
  for (i=1;i<=NAtoms;i++) {
    if( ! mol.GetAtom(i)->IsHydrogen() ) Constraints.AddAtomConstraint(i);
  }
  pFF->SetConstraints(Constraints);
}
// ---
void setGeom () {
  OBAtom *a1, *a2, *a3, *a4;
  int n;
  double angle;
  n=RotorTorsion.size() ;
  vector<ECDtor>::iterator it;
  for(it=RotorTorsion.begin();it!=RotorTorsion.end();it++) {
    a1 = mol.GetAtom((*it).a1) ;
    a2 = mol.GetAtom((*it).a2) ;
    a3 = mol.GetAtom((*it).a3) ;
    a4 = mol.GetAtom((*it).a4) ;
    angle = (*it).torsion_value ;
    mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);
  }
//  if (doDeleteBond) RestoreBond() ;
}
// ---
double molEnergy() {
  OBMol moltmp;
  vector<ECDdelBond>::iterator it1;
  double energy;

  setGeom() ;
  moltmp = mol ;
  if ( doDeleteBond ) {
    for(it1=rDelBond.begin(); it1!=rDelBond.end(); it1++) {
        moltmp.AddBond((*it1).a1,(*it1).a2,1);
    }
  }

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
  vector<ECDtor>::iterator it = RotorTorsion.begin();
  for( i=0, it=RotorTorsion.begin(); i<genome.nPhenotypes(); i++, it++) {
    (*it).torsion_value = genome.phenotype(i) ;
  }
 // value = 5000.0 - molEnergy() ;
  energy = molEnergy() ;
  value = -energy/fabs(energy) * log(1+fabs(energy)) + 10000 ;
  if ( value < 0 ) value=0;

  return value;
}
