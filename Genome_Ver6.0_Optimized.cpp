// Genome_Ver6.0_Optimized.cpp : Defines the entry point for the console application.
//

//#include <vld.h>
#include <iostream> // injecting libraries
#include <fstream> 
#include <vector>
#include <windows.h>
#include <atlstr.h>
#include <bitset>
#include "DNA_UTILS.h"
#include <time.h>
#include <iomanip>

using namespace std;

DNA_UTILS genome;

vector<int> nucCounter;

pair <vector<bitset<8>>, vector<string>> LoadFile(string); // initiating custom functions
pair <vector<bitset<8>>, vector<string>> LoadFile(string, int);
vector<string> FindFilesFromFolder(string);
wstring s2ws(const string&);

int main()
{
	clock_t t1, t2;
	t1 = clock();
	ofstream ScoreFile;
	vector < string > TestObjects = FindFilesFromFolder("D:\\temp");
	vector < DNA_UTILS::histdata > Set_of_Histograms;
	vector < vector<DNA_UTILS::positiondata> > Set_of_Locations;
	vector < pair< vector< DNA_UTILS::positiondata >, vector< DNA_UTILS::positiondata >> > Set_of_Reffers;
	vector < float > row;
	vector < vector < float > > ScoringMatrix; // declaring scoring matrix for original function 
	vector < vector< vector <float>>> newMatrices;
	vector<vector<float>> newScorings; 
	string line_1 = ",,,", line_2 = ",,,", line_3 = ",,,";

	if (TestObjects.size() >= 2)
	{
		for (int i = 0; i < TestObjects.size() - 1; i++)
		{
			cout << "////////////////////////" << endl;
			pair<vector<bitset<8>>, vector<string>> gene_A = LoadFile(TestObjects[i]);
			vector<string> gene_a_name = genome.explode(gene_A.second[4], ' ');
			if (i == 0)
			{
				line_1 = line_1 + gene_a_name[0] + ",";
				line_2 = line_2 + gene_a_name[1] + ",";
				line_3 = line_3 + gene_A.second[1] + ",";
			}
			for (int j = i + 1; j < TestObjects.size(); j++)
			{

				pair<vector<bitset<8>>, vector<string>> gene_B = LoadFile(TestObjects[j]);
				vector<string> gene_b_name = genome.explode(gene_B.second[4], ' ');
				string filename = gene_A.second[1] + "-" + gene_a_name[0] + " " + gene_a_name[1] + " VS " + gene_B.second[1] + "-" + gene_b_name[0] + " " + gene_b_name[1];
				int shift_ID = 0;
				nucCounter[0] > nucCounter[1] ? shift_ID = 1 : shift_ID = 0;
				int shifter = 0;

				if (i == 0)
				{
					line_1 = line_1 + gene_b_name[0] + ",";
					line_2 = line_2 + gene_b_name[1] + ",";
					line_3 = line_3 + gene_B.second[3] + ",";
				}
				Set_of_Histograms.push_back(genome.genHistogram(gene_A.first, gene_B.first));
				Set_of_Locations.push_back(genome.addSymbol(Set_of_Histograms[shifter], shifter));
				Set_of_Reffers.push_back(genome.sort(Set_of_Locations[shifter]));
				shifter++;
				while(shifter<4){

					// genome shifting for the next iteration
					
					if (shift_ID)
					{
						gene_A.first.clear();
						gene_A.second.clear();
						gene_a_name.clear();
						gene_A = LoadFile(TestObjects[i], shifter);
						gene_a_name = genome.explode(gene_A.second[4], ' ');

					}
					else
					{
						gene_B.first.clear();
						gene_B.second.clear();
						gene_b_name.clear();
						gene_B = LoadFile(TestObjects[j], shifter);
						gene_b_name = genome.explode(gene_B.second[4], ' ');
					}
					
					Set_of_Histograms.push_back(genome.genHistogram(gene_A.first, gene_B.first));
					Set_of_Locations.push_back(genome.addSymbol(Set_of_Histograms[shifter], shifter));
					Set_of_Reffers.push_back(genome.sort(Set_of_Locations[shifter]));
					
					shifter++;
				} 
				nucCounter.erase(nucCounter.begin() + 1);
				genome.write2File(Set_of_Histograms, Set_of_Locations, Set_of_Reffers, filename);
				row.push_back(genome.dynamicHistogram_new(Set_of_Reffers, filename));
				//scorings.push_back(genome.newaddition(Set_of_Reffers, filename));
				//newScorings.push_back(genome.newaddition2(Set_of_Reffers, gene_a_name, gene_b_name ));
				Set_of_Histograms.clear();
				Set_of_Locations.clear();
				Set_of_Reffers.clear();
			} // end of j increment (columns)
			ScoringMatrix.push_back(row);
			gene_A.first.clear();
			gene_A.second.clear();
			gene_a_name.clear();
			nucCounter.clear();
			row.clear();
		} // end of i (rows)
	}
	else
		cout << "There is no sufficient data streams\n" << endl;

	ScoreFile.open("Scoring Matrix.csv");
	ScoreFile << line_1 << endl;
	ScoreFile << line_2 << endl;
	ScoreFile << line_3 << endl;

	vector <string> col_1 = genome.explode(line_1, ',');
	vector <string> col_2 = genome.explode(line_2, ',');
	vector <string> col_3 = genome.explode(line_3, ',');

	for (int i = 0; i < ScoringMatrix.size(); i++)
	{
		ScoreFile << col_1[i] << ",";
		ScoreFile << col_2[i] << ",";
		ScoreFile << col_3[i] << ",";
		for (int j = 0; j <= i; j++)
		{
			ScoreFile << "-,";
		}

		for (int j = 0; j < ScoringMatrix[i].size(); j++)
		{
			ScoreFile << std::fixed << std::setprecision(2) << ScoringMatrix[i][j] << ",";
		}

		ScoreFile << endl;
	}

	ofstream Rfile;
	Rfile.open("RfileNeeded.csv");

	Rfile << "gene1" << ",";
	Rfile << "gene2" << ",";
	Rfile << "distance" << endl;

	for (int i = 0; i < ScoringMatrix.size(); i++)
	{
		
		//ScoreFile << col_2[i] << ",";
		//ScoreFile << col_3[i] << ",";
		for (int j = 0; j < ScoringMatrix[i].size(); j++)
		{
			Rfile << col_1[i] + " "+ col_2[i] + "(" + col_3[i] + ")" << ",";
			Rfile << col_1[i+j+1] + " "+ col_2[i+j+1] + "(" + col_3[i+j+1] + ")" << ",";
			//Rfile << col_1[i+j+1] + " "+ col_2[i+j+1] << ",";

			Rfile << ScoringMatrix[i][j] << endl;
		}

	}

	// adding the code to include sigma and sigmaSquared Scores 
	//genome.printScores(ScoringMatrix,scorings, line_1,line_2,line_3);

	// adding the code to include sigma and sigmaSquared Scores 
	//genome.printNewScores(ScoringMatrix, newScorings, line_1,line_2,line_3);

	ScoreFile.close();
	Rfile.close();
	t2 = clock();
	float diff((float)t2 - (float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
	cout << seconds << endl;
	std::system("pause");
	return 0;
}

pair <vector<bitset<8>>, vector<string>> LoadFile(string filename)
{
	vector<bitset<8>> buffer; // nucleotide buffer
	string stream = "";
	char c;
	int count = 0;
	ifstream fin(filename);

	do
	{
		fin.get(c);
		stream = stream + c;
	} while (c != '\n');

	vector<string> geneIdentifier = genome.explode(stream, '|');

	stream = "";

	if (fin.is_open())
	{
		int i = 0;
		cout << filename + " file Opened successfully!!!" << endl;
		fin.get(c);
		while (!fin.eof())
		{
			if (isalpha(c))			// checking for valid characters
			{
				count++;
				if (i < 4)
				{
					stream = stream + genome.encode(c);
					i++;
				}

				if (i == 4)
				{
					buffer.push_back(bitset<8>(string(stream)));
					stream = "";
					i = 0;
				}
			}
			fin.get(c);
		}
		cout << "char count : " << count << endl;
		nucCounter.push_back(count);
		return make_pair(buffer, geneIdentifier);
	}
	else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
}

pair <vector<bitset<8>>, vector<string>> LoadFile(string filename, int shifter)
{
	vector<bitset<8>> buffer; // nucleotide buffer
	string stream = "";
	char c;
	int count = 0;
	ifstream fin(filename);

	do
	{
		fin.get(c);
		stream = stream + c;
	} while (c != '\n');

	vector<string> geneIdentifier = genome.explode(stream, '|');

	stream = "";

	if (fin.is_open())
	{
		int i = 0;
		cout << filename + " file Opened successfully!!!" << endl;

		while (i <= shifter)
		{
			fin.get(c);
			i++;
		}

		i = 0;

		while (!fin.eof())
		{
			if (isalpha(c))			// checking for valid characters
			{
				count++;
				if (i < 4)
				{
					stream = stream + genome.encode(c);
					i++;
				}

				if (i == 4)
				{
					buffer.push_back(bitset<8>(string(stream)));
					stream = "";
					i = 0;
				}
			}
			fin.get(c);
		}
		cout << "char count : " << count << endl;
		nucCounter.push_back(count);
		return make_pair(buffer, geneIdentifier);
	}
	else //file could not be opened
	{
		cout << "File could not be opened." << endl;
	}
}

vector<string> FindFilesFromFolder(string path)								// This function return the file paths of the given folder
{
	HANDLE            hFile;
	WIN32_FIND_DATA   FindFileData;
	vector<string> fileList;

	string chFolderpath = path + "\\*.fasta";
	hFile = FindFirstFile(s2ws(chFolderpath).c_str(), &FindFileData);

	if (hFile == INVALID_HANDLE_VALUE) {

		cout << "Invalid file handle" << endl;
	}

	do
	{

		CT2CA temp(FindFileData.cFileName);
		string filename(temp);
		string filepath = path + "/" + filename;
	
		fileList.push_back(filepath);

	} while (FindNextFile(hFile, &FindFileData));

	return fileList;
}

wstring s2ws(const string& s)				// Convert string to wstring
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}