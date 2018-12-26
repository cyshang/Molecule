#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>

extern std::ofstream debug;

#define DEBUG_MOLECULE

class Molecule
{
	// ===========================================================
	// ========================= public ==========================
	// ===========================================================
public:
	// ---------- typedef ----------
	class Bond;
	class BondType;
	struct Array2;
	typedef std::map<std::string, int> str2int;
	typedef std::map<int, std::string> int2str;

	// ---------- construct ----------
	Molecule();

	// ---------- variable usage ---------- 
	inline static void usingVectorR() { ifVectorR = true; }
	inline static void usingBond() { ifBond = true; }

	// ---------- input ----------
	static void InputInfo(std::istream &);

	inline bool InputEnergy(std::istream & fin) {
		return (fin >> Energy) ? true : false;
	}
	// input X
	void InputX(std::istream &);
	// input energy_string and atom_string, then calculate X
	friend std::istream & operator >> (std::istream &, Molecule &);

	// ---------- calculate ----------
	void CalcVectorR();
	template<typename Derived>
	void CalcVectorR(const Eigen::MatrixBase<Derived> & R);
	inline double operator - (const Molecule & m) const { return (vectorR - m.vectorR).norm(); }

	// ---------- data pointer ----------
	inline double* X_ptr() { return X.data(); }
	inline double* vectorR_ptr() { return vectorR.data(); }
	
	// ---------- data reference ----------
	inline double & getEnergy() { return Energy; }	
	inline Eigen::MatrixXd & getX() { return X; }
	inline Eigen::VectorXd & getVectorR() { return vectorR; }
	inline std::vector<std::vector<Bond>> & getBond() { return bond; }

	// ---------- output ----------
	friend std::ostream & operator << (std::ostream &, const Molecule &);

	// ---------- molecule description ----------
	static std::string name;
	
	static int nElem;
	static std::vector<std::string> elem_list;
	static str2int elem2num;
	static int2str num2elem;

	static int totAtom;
	static std::vector<int> nAtom;	// number of atoms for each element
	static std::vector<int> atom_list;
	static std::vector<std::vector<int>> atomTravlist;

	// ---------- bond description ----------
	static int totBond;
	static int nBondtype;
	static std::vector<BondType> bondTypelist;
	static std::vector<int> nBond;
	static std::vector<std::vector<Array2>> bondTravlist;

	// ===========================================================
	// ========================= private =========================
	// ===========================================================
private:
	// ---------- molecule string data ----------
	std::string energy_str;
	std::vector<std::string> atom_str;

	// ---------- molcule numerical data ----------
	double Energy;
	Eigen::MatrixXd X;
	Eigen::VectorXd vectorR;

	//Eigen::MatrixXd matrixR;
	//Eigen::MatrixXd matrixR2;
	//std::vector<Eigen::MatrixXd> cos0;

	// ---------- bond data ----------
	std::vector<std::vector<Bond>> bond;

	// ---------- variable usage ---------- 
	static bool ifVectorR;
	static bool ifBond;
	static void BondInfo();
};

// ========== template functions ==========

template<typename Derived>
void Molecule::CalcVectorR(const Eigen::MatrixBase<Derived> & R)
{
	int pos = 0;
	for (int i = 0; i < totAtom - 1; ++i) {
		for (int j = i + 1; j < totAtom; ++j) {
			const_cast<Eigen::MatrixBase<Derived>&>(R)(pos++) = (X.col(i) - X.col(j)).norm();
		}
	}

	return ;
}

// ========================================

class Molecule::BondType
{
public:
	BondType(const int & iE, const int & jE) {
		iElem = (iE < jE) ? iE : jE;
		jElem = (iE < jE) ? jE : iE;
	}
	
	inline bool operator == (const BondType type) { return iElem == type.iElem && jElem == type.jElem; }
	inline bool operator != (const BondType type) { return iElem != type.iElem || jElem != type.jElem; }
	friend inline std::ostream & operator << (std::ostream & os, const BondType & t) { 
		os << Molecule::num2elem[t.iElem] << '-' << Molecule::num2elem[t.jElem];
		return os;
	}

private:
	int iElem;
	int jElem;
};

struct Molecule::Array2
{
	int iAtom;
	int jAtom;

	Array2(const int & i, const int & j) :iAtom(i), jAtom(j) {}
};

class Molecule::Bond
{
public:
	Bond() :len(0.0), iAtom(0), jAtom(0) {}
	Bond(const double & len, const int & i, const int & j) :len(len), iAtom(i), jAtom(j) {}

	friend inline bool operator < (const Bond & a, const Bond & b) {
		return a.len < b.len;
	}

	friend inline bool operator > (const Bond & a, const Bond & b) {
		return a.len > b.len;
	}

	friend inline std::ostream & operator << (std::ostream & os, const Bond & a) {
		os << a.iAtom << ", " << a.jAtom << ", " << a.len;
		return os;
	}

private:
	double len;
	int iAtom;
	int jAtom;
};

#endif // !MOLECULE_H_