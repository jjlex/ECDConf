/* ----------------------------------------------------------------------------
   GALib code
---------------------------------------------------------------------------- */
// ----  For GALib -----------
#include <stdio.h>
#include <math.h>

#include <ga/GASimpleGA.h>
#include <ga/GABin2DecGenome.h>
#include <ga/std_stream.h>


//#define cout STD_COUT
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
