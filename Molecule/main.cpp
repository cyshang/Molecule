#include <iostream>
#include <fstream>
#include "Molecule.h"

using namespace std;

ofstream debug;


int main()
{
#ifdef DEBUG_MOLECULE
	debug.open("debug.dat", ofstream::out);
#endif // DEBUG_MOLECULE

	ifstream fin;
	fin.open("OH3", ifstream::in);

	Molecule::usingBond();
	Molecule::InputInfo(fin);

	fin.close();


#ifdef DEBUG_MOLECULE
	debug.close();
#endif // DEBUG_MOLECULE

    return 0;
}