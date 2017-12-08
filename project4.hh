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
  int best_score = 0, best_i = 0, up = 0, left = 0, diag = 0, max = 0;
  bool done = false;
	matchString1 = "";
	matchString2 = "";
  int len_str1 = string1.length();
  int len_str2 = string2.length();

  int dp[len_str1 + 1][len_str2 + 1]; // Dynamic programming array
  char bp[len_str1 + 1][len_str2 + 1]; // Back pointers array

  // Initialize 2-D arrays
  for(int i = 0; i < len_str1 + 1; i++)
  {
    for(int j = 0; j < len_str2 + 1; j++)
    {
      dp[i][j] = 0;
      bp[i][j] = '?';
    }
  }

  // Set dynamic programming array values including penalties
  for(int i = 1; i < len_str1 + 1; i++)
  {
    for(int j = 1; j < len_str2 + 1; j++)
    {
      up = dp[i - 1][j] + bpa.get_penalty(string1[i - 1], '*');
      left = dp[i][j - 1] + bpa.get_penalty('*', string2[j - 1]);
      diag = dp[i - 1][j - 1] + bpa.get_penalty(string1[i - 1], string2[j - 1]);

      // Set back pointers array values accordingly
      if(left > up)
      {
        if(left > diag)
        {
          bp[i][j] = 'l';
        }
        else
        {
          bp[i][j] = 'd';
        }
      }
      else
      {
        if(up > diag)
        {
          bp[i][j] = 'u';
        }
        else
        {
          bp[i][j] = 'd';
        }
      }

      // Find the largest value and set dynamic array, if none, use zero.
      max = 0;
      if(up > max)
      {
        max = up;
      }
      if(left > max)
      {
        max = left;
      }
      if(diag > max)
      {
        max = diag;
      }
      dp[i][j] = max;
    }
  }

  // Locate best score from bottom row
  best_i = len_str1;
  int best_j = 0;
  for(int j = 1; j < len_str2 + 1; j++)
  {
    if(dp[best_i][j] > best_score)
    {
      best_score = dp[best_i][j];
      best_j = j;
    }
  }

  // Follow back pointers to get alignment
  int i = best_i;
  int j = best_j;

  while(!done)
  {
    if(bp[i][j] == 'u')
    {
      matchString1 += string1[i - 1];
      matchString2 += '*';
      i--;
    }
    else if(bp[i][j] == 'l')
    {
      matchString1 += '*';
      matchString2 += string2[j - 1];
      j--;
    }
    else if(bp[i][j] == 'd')
    {
      matchString1 += string1[i - 1];
      matchString2 += string2[j - 1];
      i--;
      j--;
    }
    else if(bp[i][j] == '?')
    {
      done = true;
    }
  }

  // Return matches in appropriate order
  reverse(matchString1.begin(), matchString1.end());
  reverse(matchString2.begin(), matchString2.end());
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
  int best_score = 0, score = 0, best_i = 0;
  matchString1 = "";
  matchString2 = "";
  std::string runningBest1 = "";
  std::string runningBest2 = "";

  // Traverse proteins vector to locate best match
  for(int i = 0; i < proteins.size(); i++)
  {
    score = local_alignment(proteins[i]->sequence, string1, bpa, matchString1, matchString2);
    if(score > best_score)
    {
      // Store best's details
      best_score = score;
      best_i = i;
      runningBest1 = matchString1;
      runningBest2 = matchString2;
    }
  }
  matchString1 = runningBest1;
  matchString2 = runningBest2;
  std::cout << string1 << ": best score - " << best_score << std::endl;
	return best_protein = proteins[best_i];
}
