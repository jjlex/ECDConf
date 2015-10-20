/**********************************************************************
obrotate = rotate a tortional bond matched by a SMART pattern
Copyright (C) 2003 Fabien Fontaine
Some portions Copyright (C) 2004-2005 Geoffrey R. Hutchison
Some portions Copyright (C) 2008 Tim Vandermeersch
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/


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

using namespace std;
using namespace OpenBabel;

///////////////////////////////////////////////////////////////////////////////
//! \brief Set a tortional bond to a given angle
int main(int argc,char **argv)
{
  OBAtom *a1, *a2, *a3, *a4;
//  unsigned int smartor[4]= {0,0,0,0};// atoms of the tortional in the SMART
  unsigned long atomid[4]= {0,0,0,0};// atoms of the tortional
  float angle =   0;      // tortional angle value to set in degree
  char *FileIn =NULL, *Pattern=NULL;
  unsigned int i, t, errflg = 0;
  int c;
  string err;
  //bool changeAll = false; // default to only change the last matching torsion

  // parse the command line -- optional -a flag to change all matching torsions
  if (argc < 7 || argc > 8) {
    errflg++;
  } else {
    FileIn = argv[1];
    // Read the atom position
      for(i=2, t=0; i<6; ++i, ++t) {
      c = sscanf(argv[i], "%lu", &atomid[t]);
      if (c != 1) {
        break;
        errflg++;
      }
    }
    c = sscanf(argv[6], "%f", &angle);
  }

  if (errflg) {
//    cerr << "Usage: obrotate \"PATTERN\" <filename> <atom1> <atom2> <atom3> <atom4> <angle>" << endl;
    cerr << "Usage: obrotate <filename> <atom1> <atom2> <atom3> <atom4> <angle>" << endl;
    exit(-1);
  }

  // create pattern
/*
  OBSmartsPattern sp;
  sp.Init(Pattern);
  if (sp.NumAtoms() < 4) {
    cerr << "obrotate: The number of atoms in the SMART pattern must be higher than 3." << endl;
    exit(-1);
  }

  for (i=0; i<4; ++i) {
    if ( smartor[i] < 1 || smartor[i] > sp.NumAtoms()) {
      cerr << "obrotate: The torsional atom values must be between 1 and "
           <<  sp.NumAtoms()
           << ", which is the number of atoms in the SMART pattern." << endl;
      exit(-1);
    }
  }

*/
  OBConversion conv; //NF...
  OBFormat* format = conv.FormatFromExt(FileIn);
  if(!(format && conv.SetInAndOutFormats(format, format))) { //in and out formats same
    cerr << "obrotate: cannot read and/or write this file format!" << endl;
    exit (-1);
  } //...NF

  //Open the molecule file
  ifstream ifs;

  // Read the file
  ifs.open(FileIn);
  if (!ifs) {
    cerr << "obrotate: cannot read input file!" << endl;
    exit (-1);
  }

  OBMol mol; 
  vector< vector <int> > maplist;      // list of matched atoms
  vector< vector <int> >::iterator m;  // and its iterators
  //   int tindex;

  // Set the angles
  for (;;) {
    mol.Clear();
    //NF      ifs >> mol;                   // Read molecule
    conv.Read(&mol,&ifs); //NF
    if (mol.Empty())
      break;

/*
    if (sp.Match(mol)) {          
      // if match perform rotation
      maplist = sp.GetUMapList(); // get unique matches
      
      if (maplist.size() > 1)
        cerr << "obrotate: Found " << maplist.size() << " matches. Only last one will be rotated." << endl;

      // look at all the mapping atom but save only the last one.
      for (m = maplist.begin(); m != maplist.end(); ++m) {
        a1 = mol.GetAtom( (*m)[ smartor[0] - 1] );
        a2 = mol.GetAtom( (*m)[ smartor[1] - 1] );
        a3 = mol.GetAtom( (*m)[ smartor[2] - 1] );
        a4 = mol.GetAtom( (*m)[ smartor[3] - 1] );
        //if (changeAll)
        //  mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);                
      }

      if ( !a2->IsConnected(a3) ) {
        cerr << "obrotate: The atoms of the rotating bond must be bonded." << endl;
        exit(-1);
      }

      mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);
    } else {
      cerr << "obrotate: Found 0 matches for the SMARTS pattern." << endl;
      exit(-1);
    }
*/
// -- newly added
    a1 = mol.GetAtomById( atomid[0] - 1 );
    a2 = mol.GetAtomById( atomid[1] - 1 );
    a3 = mol.GetAtomById( atomid[2] - 1 );
    a4 = mol.GetAtomById( atomid[3] - 1 );
    mol.SetTorsion(a1, a2, a3, a4, angle * DEG_TO_RAD);

    cerr << a1->GetIdx() << endl;
    cerr << a2->GetIdx() << endl;
    cerr << a3->GetIdx() << endl;
    cerr << a4->GetIdx() << endl;
    cerr << "" << endl;
    cerr << a1->GetAtomicNum() << endl;
    cerr << a2->GetAtomicNum() << endl;
    cerr << a3->GetAtomicNum() << endl;
    cerr << a4->GetAtomicNum() << endl;
    exit(1);
//--
    //NF      cout << mol;
    conv.Write(&mol,&cout); //NF
  }

  return(0);
}


