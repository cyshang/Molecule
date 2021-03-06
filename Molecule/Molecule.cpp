#include "Molecule.h"
#include <sstream>

using namespace Eigen;
using std::vector;
using std::string;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::endl;

// ---------- define static member variable ----------
string Molecule::name;

int Molecule::nElem;
vector<string> Molecule::elem_list;
Molecule::str2int Molecule::elem2num;
Molecule::int2str Molecule::num2elem;

int Molecule::totAtom;
vector<int> Molecule::nAtom;
vector<int> Molecule::atom_list;
vector<vector<int>> Molecule::atomTravlist;

int Molecule::totBond;
int Molecule::nBondtype;
vector<Molecule::BondType> Molecule::bondTypelist;
vector<int> Molecule::nBond;
vector<vector<Molecule::Array2>> Molecule::bondTravlist;

bool Molecule::ifVectorR = false;
bool Molecule::ifBond = false;

// ---------- construct function ---------- 
Molecule::Molecule()
	:energy_str(), atom_str(),
	Energy(0.0),
	X(3, totAtom),
	vectorR(),
	bond()
{
	if (ifVectorR)
		vectorR.resize(totAtom * (totAtom - 1) / 2);

	if (ifBond) {
		bond.resize(nBondtype);
		for (int iBondtype = 0; iBondtype < nBondtype; ++iBondtype) {
			bond[iBondtype].resize(nBond[iBondtype]);
		}
	}
}

void Molecule::InputInfo(istream & fin)
{
	string line;
	while (std::getline(fin, line)) {
		if (line[0] == '#' || line[0] == '\n' || line[0] == ' ')
			continue;
		
		for (size_t i = 0; i < line.size(); ++i) {
			if (line[i] == '=' || line[i] == '(' || line[i] == ')' || line[i]=='\'' || line[i]=='"')
				line[i] = ' ';
		}
		
		istringstream sin(line);
		string key_word;

		sin >> key_word;
		if (key_word == "molecule") {
			sin >> name;
		}
		else if (key_word == "elem_num") {
			sin >> nElem;
		}
		else if (key_word == "elem_list") {
			elem_list.resize(nElem);

			for (size_t i = 0; i < elem_list.size(); ++i) {
				sin >> elem_list[i];
				elem2num[elem_list[i]] = i;
				num2elem[i] = elem_list[i];
			}
		}
		else if (key_word == "atom_num") {
			sin >> totAtom;
		}
		else if (key_word == "atom_list") {
			atom_list.resize(totAtom);

			nAtom.resize(nElem);
			for (size_t i = 0; i < nAtom.size(); ++i) {
				nAtom[i] = 0;
			}

			string elem;
			for (size_t i = 0; i < atom_list.size(); ++i) {
				sin >> elem;
				atom_list[i] = elem2num[elem];
				nAtom[atom_list[i]]++;
			}
		}
	}

	atomTravlist.resize(nElem, vector<int>());
	for (int iAtom = 0; iAtom < totAtom; ++iAtom) {
		atomTravlist[atom_list[iAtom]].push_back(iAtom);
	}

#ifdef DEBUG_MOLECULE

	debug << "name: " << name << endl;
	debug << "nElem: " << nElem << endl;
	debug << "elem_list: ";
	for (auto & i : elem_list) debug << i << ' ';
	debug << endl;
	for (auto & i : elem_list) {
		debug << elem2num[i] << "<->" << num2elem[elem2num[i]] << endl;
	}
	debug << "totAtom: " << totAtom << endl;
	debug << "nAtom: ";
	for (auto & i : nAtom) debug << i << ' ';
	debug << endl;
	debug << "atom_list: ";
	for (auto & i : atom_list) debug << i << ' ';
	debug << endl;
	debug << "atomTravlist:" << endl;
	for (size_t i = 0; i < atomTravlist.size(); ++i) {
		debug << num2elem[i] << ": ";
		for (const auto & j : atomTravlist[i])
			debug << j << ' ';
		debug << endl;
	}
#endif // DEBUG_MOLECULE


	if (ifBond)
		BondInfo();
}

void Molecule::BondInfo()
{
	totBond = (totAtom * (totAtom - 1)) / 2;

	nBondtype = 0;
	bondTypelist.clear();
	for (int iE = 0; iE < nElem; ++iE) {
		if (nAtom[iE] > 1) {
			nBondtype++;
			bondTypelist.push_back(BondType(iE, iE));
		}
		for (int jE = iE + 1; jE < nElem; ++jE) {
			nBondtype++;
			bondTypelist.push_back(BondType(iE, jE));
		}
	}

	bondTravlist.clear();
	nBond.clear();
	int iBondtype = 0;
	for (int iE = 0; iE < nElem; ++iE) {
		if (nAtom[iE] > 1) {
			nBond.push_back((nAtom[iE] * (nAtom[iE] - 1)) / 2);
			bondTravlist.push_back(vector<Array2>());
			for (int i = 0; i < nAtom[iE] - 1; ++i) {
				for (int j = i + 1; j < nAtom[iE]; ++j) {
					bondTravlist[iBondtype].push_back(Array2(atomTravlist[iE][i], atomTravlist[iE][j]));
				}
			}
			iBondtype++;
		}
		for (int jE = iE + 1; jE < nElem; ++jE) {
			nBond.push_back(nAtom[iE] * nAtom[jE]);
			bondTravlist.push_back(vector<Array2>());
			for (int i = 0; i < nAtom[iE]; ++i) {
				for (int j = 0; j < nAtom[jE]; ++j) {
					bondTravlist[iBondtype].push_back(Array2(atomTravlist[iE][i], atomTravlist[jE][j]));
				}
			}
			iBondtype++;
		}
	}

#ifdef DEBUG_MOLECULE

	debug << "totBond: " << totBond << endl;
	debug << "nBondtype: " << nBondtype << endl;
	debug << "bondTypelist: ";
	for (const auto & i : bondTypelist) debug << i << ' ';
	debug << endl;
	debug << "nBond: ";
	for (const auto & i : nBond) debug << i << ' ';
	debug << endl;
	debug << "bondTravlist:" << endl;
	for (const auto & i : bondTravlist) {
		for (const auto & j : i) {
			debug << '(' << j.iAtom << ", " << j.jAtom << ") ";
		}
		debug << endl;
	}

#endif // DEBUG_MOLECULE

}

void Molecule::InputX(std::istream & fin)
{
	for (int i = 0; i < totAtom; ++i) {
		string elem;
		fin >> elem;

		for (int j = 0; j < 3; ++j) {
			fin >> X(j, i);
		}
	}
}

std::istream & operator >> (std::istream & fin, Molecule & m)
{
	using std::getline;
	
	string line;
	if (getline(fin, line)) {
		m.atom_str.resize(m.totAtom);
		m.X.resize(3, m.totAtom);

		getline(fin, m.energy_str);
		for (int i = 0; i < m.totAtom; ++i) {
			getline(fin, m.atom_str[i]);

			istringstream sin(m.atom_str[i]);
			string elem;
			sin >> elem;
			for (int j = 0; j < 3; ++j) {
				sin >> m.X(j, i);
			}
		}
	}

	return fin;
}

void Molecule::CalcVectorR()
{
	int pos = 0;
	for (int i = 0; i < totAtom - 1; ++i) {
		for (int j = i + 1; j < totAtom; ++j) {
			vectorR(pos++) = (X.col(i) - X.col(j)).norm();
		}
	}
}

ostream & operator << (ostream & fout, const Molecule & m)
{
	fout << m.totAtom << endl;
	fout << m.energy_str << endl;
	for (size_t i = 0; i < m.atom_str.size(); ++i) {
		fout << m.atom_str[i] << endl;
	}

	return fout;
}