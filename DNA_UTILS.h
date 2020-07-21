#pragma once
#ifndef DNA_UTILS_H
#define DNA_UTILS_H

#include <vector>
#include <fstream>
#include <bitset>
#include <iostream>
using namespace std;

class DNA_UTILS
{
private:
	int count_for_1_a[4], count_for_equals[4], count_for_1_b[4];
	
public:
	DNA_UTILS();
	~DNA_UTILS();

	struct differ{
		int gap;
		int freaqency;
		int direct;
		int neighbor;
	};

	struct histdata{
		int hist_a[256];
		int hist_b[256];
		int pos_a[256];
		int pos_b[256];
	};

	struct positiondata{
		int index_a;
		int index_b;
		int a_pos;
		int b_pos;
		int cast;
	};

	histdata genHistogram(vector<bitset<8>>, vector<bitset<8>>);
	string encode(char);
	vector<string> explode(const string&, const char&);			// Equals to the split function in JAVA. Return a vector of strings delimited with given delimiter		
	int findMember(int[], bitset<8>);
	string bitsSelector(string);
	string decode_eight(bitset<8>);
	char decode(string);
	vector<positiondata> addSymbol(histdata, int);
	pair<vector<positiondata>, vector<positiondata>> sort(vector<positiondata>);
	void write2File(vector <histdata>, vector <vector<positiondata>>, vector < pair<vector<positiondata>, vector<positiondata>> >, string);
	float dynamicHistogram(vector < pair<vector<positiondata>, vector<positiondata>> >, string);
	float dynamicHistogram_new(vector < pair<vector<positiondata>, vector<positiondata>> >, string);
	pair<vector<int>,vector<int>> newaddition(vector < pair<vector<positiondata>, vector<positiondata>> >, string);
	vector<float> newaddition2(vector < pair<vector<positiondata>, vector<positiondata>> >, vector<string>,vector<string>);
	int getMemberIndex(vector<differ>, int);
	int getMaximum(int, int);
	int mergesort(positiondata *, int);
	void mergePortion(positiondata *, int, int, positiondata *);
	int getMax(vector<differ>);
	int getSigma(vector<differ>);
	float getNewScore(vector<differ>);
	int getsigmasquared(vector<differ>);
	void printNewScores(vector < vector < float > >, vector<vector<float>> , string, string, string);
	float takemax(float []);
};

#endif