///////////////////////////////////////////////////////////////////////////////
// maxprotein.hh
//
// Compute the set of foods that maximizes protein, within a calorie budget,
// with the greedy method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>

// Simple structure for a single protein
struct Protein {
	Protein() {
		description = "";
		sequence = "";
	}
	Protein(std::string desc, std::string seq) {
		description = desc;
		sequence = seq;
	}
	std::string		description;
	std::string 	sequence;
};

// class for BLOSUM penalties.. acts as a matrix holding penalties based 
//     on transitions for one amino acid to another
class BlosumPenaltyArray {
public:
	BlosumPenaltyArray() {
		// nothing here
	}
	~BlosumPenaltyArray() {
		// nothing here
	}
	BlosumPenaltyArray(BlosumPenaltyArray & that) {
		internal_copy(that);
	}
	BlosumPenaltyArray & operator=(BlosumPenaltyArray & that) {
		internal_copy(that);
		return *this;
	}

	int get_penalty(char c1, char c2) {
		return _penaltyMap[c1][c2];
	}

	void set_penalty(char c1, char c2, int penalty) {
		_penaltyMap[c1][c2] = penalty;
	}

	void debug_map() {
		for (auto itr1 = _penaltyMap.begin(); itr1 != _penaltyMap.end(); itr1++) {
			for (auto itr2 = itr1->second.begin(); itr2 != itr1->second.end(); itr2++) {
				std::cout << itr2->second << "  ";
			}
			std::cout << std::endl;
		}
	}

private:
	void internal_copy(BlosumPenaltyArray & that) {
		this->_penaltyMap = that._penaltyMap;
	}

	std::map<char, std::map<char, int>> _penaltyMap;
};


// Alias for a vector of shared pointers to Protein objects.
typedef std::vector<std::shared_ptr<Protein>> ProteinVector;


// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_proteins(ProteinVector & proteins, const std::string& path) 
{
  //std::cout << "Loading proteins from [" << path << "]" << std::endl;
  proteins.clear();
  std::ifstream ifs(path.c_str());
  if (!ifs.is_open() || !ifs.good()) {
    std::cout << "Failed to open [" << path << "]" << std::endl;
    return false;
  }
  int proteinsLoaded = 0;
  bool have_description = false;
  std::shared_ptr<Protein> newProtein = nullptr;
  while (!ifs.eof()) {
    std::string lineBuffer;
    std::getline(ifs, lineBuffer);
    if (ifs.eof()) {
      break;
    }
    if (lineBuffer.size() == 0) {
		continue;
	}
    if (lineBuffer[0] == '>') {
		newProtein = std::shared_ptr<Protein>(new Protein);
		newProtein->description = lineBuffer.substr(1);
        have_description = true;
    } else if (have_description) {
		newProtein->sequence = lineBuffer;
	    proteins.push_back(newProtein);
        proteinsLoaded++;
        have_description = false;
    }
  }

	ifs.close();

  return true;
}

// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool save_proteins(ProteinVector & proteins, const std::string& path) 
{
	std::cout << "Saving proteins from [" << path << "]" << std::endl;
	std::ofstream ofs(path.c_str());
	if (!ofs.is_open() || !ofs.good()) {
		std::cout << "Failed to open [" << path << "]" << std::endl;
		return false;
	}

	for (int i = 0; i < proteins.size(); i++) {
		ofs << proteins[i]->description << std::endl;
		ofs << proteins[i]->sequence.substr(10,10) << std::endl;
	}

	ofs.close();

  return true;
}

// -------------------------------------------------------------------------
// Load the BLOSUM penalties from a standard BLOSUM file (matrix format)
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_blosum_file(BlosumPenaltyArray & bpa, const std::string& path) 
{
  std::ifstream ifs(path.c_str());
  if (!ifs.is_open() || !ifs.good()) {
    std::cout << "Failed to open [" << path << "]" << std::endl;
    return false;
  }

  std::vector<char> aas; // Create vector to hold our Amino Acids

  while (!ifs.eof()) {
    std::string lineBuffer;
    std::getline(ifs, lineBuffer);
    if (ifs.eof()) {
      break;
    }
    if (lineBuffer.size() == 0) {
		continue;
	}

    if (lineBuffer[0] == '$') {
		std::string buf;
		std::stringstream ss(lineBuffer.substr(1)); // Insert the string into a stream
	    while (ss >> buf) {
	        aas.push_back(buf[0]);
		}
		continue;
	}

	int penalty;
	char thisRowChar = lineBuffer[0];
	std::stringstream ss(lineBuffer.substr(1)); // Insert the string into a stream
	int tokenCount = 0;
    while (ss >> penalty) {
        bpa.set_penalty(thisRowChar, aas[tokenCount], penalty);
		tokenCount++;
	}
  }

  return true;
}

// -------------------------------------------------------------------------
int local_alignment(const std::string & string1, 
					const std::string & string2, 
					BlosumPenaltyArray & bpa,
					std::string & matchString1, 
					std::string & matchString2)
{
	int best_score = 0;
	matchString1 = "";
	matchString2 = "";
	return best_score;
}

// -------------------------------------------------------------------------
std::shared_ptr<Protein> local_alignment_best_match(
					ProteinVector & proteins, 
					const std::string & string1,
					BlosumPenaltyArray & bpa,
					std::string & matchString1, 
					std::string & matchString2)
{
	std::shared_ptr<Protein> best_protein = nullptr;
	matchString1 = "";
	matchString2 = "";

	return best_protein;
}


