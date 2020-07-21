#include "DNA_UTILS.h"
#include "btree.h"
#include "pos_btree.h"
using namespace std;

DNA_UTILS::DNA_UTILS()
{
}


DNA_UTILS::~DNA_UTILS()
{
}


DNA_UTILS::histdata DNA_UTILS::genHistogram(vector<bitset<8>>gene_a, vector<bitset<8>>gene_b)
{
	int chk_pos_a[256], chk_pos_b[256];
	DNA_UTILS::histdata results;

	for (int i = 0; i < 256; i++)
	{
		results.hist_a[i] = 0;
		results.hist_b[i] = 0;
		results.pos_a[i] = 0;
		results.pos_b[i] = 0;
		chk_pos_a[i] = -1;
		chk_pos_b[i] = -1;
	}

	for (int i = 0; i < (int)gene_a.size(); i++)
	{
		results.hist_a[gene_a[i].to_ulong()]++;
		if (chk_pos_a[gene_a[i].to_ulong()] == -1)
		{
			results.pos_a[gene_a[i].to_ulong()] = i * 4;
			chk_pos_a[gene_a[i].to_ulong()] = 0;
		}
	}
	for (int i = 0; i < (int)gene_b.size(); i++)
	{
		results.hist_b[gene_b[i].to_ulong()]++;
		if (chk_pos_b[gene_b[i].to_ulong()] == -1)
		{
			results.pos_b[gene_b[i].to_ulong()] = i * 4;
			chk_pos_b[gene_b[i].to_ulong()] = 0;
		}
	}

	return results;
}

int DNA_UTILS::findMember(int hist[], bitset<8> base)
{
	string base_str = base.to_string<char, char_traits<char>, allocator<char>>();
	for (int i = 0; i < (int)base_str.length(); i += 2)
	{
		string temp_base = base_str;
		string replicants = bitsSelector(base_str.substr(i, 2));
		for (int j = 0; j < (int)replicants.length(); j += 2)
		{
			temp_base.replace(i, 2, replicants.substr(j, 2));
			int temp_index = bitset<8>(string(temp_base)).to_ullong();
			if (hist[temp_index] == 1)
				return temp_index;
		}

	}

	return -1;
}

vector<DNA_UTILS::positiondata> DNA_UTILS::addSymbol(DNA_UTILS::histdata histo, int shift)
{
	DNA_UTILS::positiondata instantData;
	count_for_1_a[shift] = 0, count_for_equals[shift] = 0, count_for_1_b[shift] = 0;
	vector<DNA_UTILS::positiondata> reffer_A;
	pos_btree poTree = pos_btree();
	for (int i = 0; i < 256; i++)
	{
		if (histo.hist_a[i] == 1 && histo.hist_b[i] == 1 /*&& histo.hist_a[i] != 0 && histo.hist_b[i] != 0*/)
		{
			instantData.index_a = i;
			instantData.index_b = i;
			instantData.a_pos = histo.pos_a[i];
			instantData.b_pos = histo.pos_b[i];
			instantData.cast = 0;
			poTree.insert(instantData);
			count_for_equals[shift]++;
		}
		else if (histo.hist_a[i] == 1 && histo.hist_b[i] == 0)
		{
			int in_b = findMember(histo.hist_b, bitset<8>(int(i)));
			if (in_b != -1)
			{
				instantData.index_a = i;
				instantData.index_b = in_b;
				instantData.a_pos = histo.pos_a[i];
				instantData.b_pos = histo.pos_b[in_b];
				instantData.cast = 1;
				poTree.insert(instantData);

			}
		}
		else if (histo.hist_a[i] == 0 && histo.hist_b[i] == 1)
		{
			int in_a = findMember(histo.hist_a, bitset<8>(int(i)));
			if (in_a != -1)
			{
				instantData.index_a = in_a;
				instantData.index_b = i;
				instantData.a_pos = histo.pos_a[in_a];
				instantData.b_pos = histo.pos_b[i];
				instantData.cast = 1;
				poTree.insert(instantData);
			}
		}

		if (histo.hist_a[i] == 1)
			count_for_1_a[shift]++;

		if (histo.hist_b[i] == 1)
			count_for_1_b[shift]++;
	}

	reffer_A = poTree.printTree();
	return reffer_A;
}

pair<vector<DNA_UTILS::positiondata>, vector<DNA_UTILS::positiondata>> DNA_UTILS::sort(vector<DNA_UTILS::positiondata> reffer_A)
{
	vector<DNA_UTILS::positiondata> reffer_B = reffer_A;

	mergesort(reffer_B.data(), reffer_B.size());
	return make_pair(reffer_A, reffer_B);

}

int DNA_UTILS::mergesort(DNA_UTILS::positiondata *input, int size)
{
	DNA_UTILS::positiondata *scratch = (DNA_UTILS::positiondata *)malloc(size * sizeof(DNA_UTILS::positiondata));

		if (scratch != NULL)
		{
			mergePortion(input, 0, size, scratch);
			free(scratch);
			return 1;
		}
		else
		{
			return 0;
		}

	
}

void DNA_UTILS::mergePortion(DNA_UTILS::positiondata *input, int left, int right, DNA_UTILS::positiondata *scratch)
{

	if (right == left + 1)
	{
		return;
	}
	else
	{
		int i = 0;
		int length = right - left;
		int midpoint_distance = length / 2;

		int l = left, r = left + midpoint_distance;


		mergePortion(input, left, left + midpoint_distance, scratch);
		mergePortion(input, left + midpoint_distance, right, scratch);


		for (i = 0; i < length; i++)
		{
			if (l < left + midpoint_distance && (r == right || getMaximum(input[l].b_pos, input[r].b_pos) == input[l].b_pos))
			{
				scratch[i] = input[l];
				l++;
			}
			else
			{
				scratch[i] = input[r];
				r++;
			}
		}

		for (i = left; i < right; i++)
		{
			input[i] = scratch[i - left];
		}
	}
}

void DNA_UTILS::write2File(vector <DNA_UTILS::histdata>histograms, vector <vector<DNA_UTILS::positiondata>> locations,
	vector < pair<vector<DNA_UTILS::positiondata>, vector<DNA_UTILS::positiondata>> > reffers, string filename)
{
	ofstream histFile, resFile, sortFile_A, sortFile_B;
	histFile.open(".\\Histograms\\" + filename + ".csv");
	resFile.open(".\\Results_unsorted\\" + filename + ".csv");
	sortFile_A.open(".\\Results_sorted_preffered_A\\" + filename + ".csv");
	sortFile_B.open(".\\Results_sorted_preffered_B\\" + filename + ".csv");

	histFile << ",Sequance A,,,,Sequance B" << endl;
	histFile << "Shift,Index,Capacity,Bar Chart,Position,Capacity,Bar Chart,Position" << endl;

	resFile << "Shift,Index A,Symbol A,Position A,Index B,Symbol B,Position B,Number of ones in A,Number of ones in B,Number of equals" << endl;

	sortFile_A << "Shift,Index A,Symbol A,Position A,Index B,Symbol B,Position B,Number of ones in A, Number of ones in B,Number of equals" << endl;
	sortFile_B << "Shift,Index A,Symbol A,Position A,Index B,Symbol B,Position B,Number of ones in A, Number of ones in B,Number of equals" << endl;

	for (int it = 0; it < 4; it++)
	{
		// Begin write of histogram file
		for (int i = 0; i < 256; i++)
		{
			if (i == 0)
				histFile << it;
			histFile << "," << i << "," << histograms[it].hist_a[i] << ",";
			if (histograms[it].hist_a[i] > 0)
			{
				for (int j = 1; j <= histograms[it].hist_a[i]; j++)
					histFile << "*";
			}
			if (histograms[it].pos_a[i] > 0)
				histFile << "," << histograms[it].pos_a[i];
			else
				histFile << ",-";
			histFile << "," << histograms[it].hist_b[i] << ",";
			if (histograms[it].hist_b[i] > 0)
			{
				for (int j = 1; j <= histograms[it].hist_b[i]; j++)
					histFile << "*";
			}
			if (histograms[it].pos_b[i] > 0)
				histFile << "," << histograms[it].pos_b[i];
			else
				histFile << ",-";
			histFile << endl;
		}

		// Begin write of unsorted results file

		for (int i = 0; i < locations[it].size(); i++)
		{
			if (i == 0)
				resFile << it;

			resFile << "," << locations[it][i].index_a << "," << decode_eight(bitset<8>(int(locations[it][i].index_a))) << "," << locations[it][i].a_pos << ","
				<< locations[it][i].index_b << "," << decode_eight(bitset<8>(int(locations[it][i].index_b))) << "," << locations[it][i].b_pos;

			if (i == 0)
				resFile << "," << count_for_1_a[it] << "," << count_for_1_b[it] << "," << count_for_equals[it];
			resFile << endl;

		}

		// Begin write of sorted results according to sequence A file

		for (int i = reffers[it].first.size() - 1; i >= 0; i--)
		{
			if (i == reffers[it].first.size() - 1)
				sortFile_A << it;
			sortFile_A << "," << reffers[it].first[i].index_a << "," << decode_eight(bitset<8>(int(reffers[it].first[i].index_a))) << "," << reffers[it].first[i].a_pos << ","
				<< reffers[it].first[i].index_b << "," << decode_eight(bitset<8>(int(reffers[it].first[i].index_b))) << "," << reffers[it].first[i].b_pos;

			if (i == reffers[it].first.size() - 1)
				sortFile_A << "," << count_for_1_a[it] << "," << count_for_1_b[it] << "," << count_for_equals[it];
			sortFile_A << endl;
		}

		// Begin write of sorted results according to sequence B file

		for (int i = reffers[it].second.size() - 1; i >= 0; i--)
		{
			if (i == reffers[it].second.size() - 1)
				sortFile_B << it;
			sortFile_B << "," << reffers[it].second[i].index_a << "," << decode_eight(bitset<8>(int(reffers[it].second[i].index_a))) << "," << reffers[it].second[i].a_pos << ","
				<< reffers[it].second[i].index_b << "," << decode_eight(bitset<8>(int(reffers[it].second[i].index_b))) << "," << reffers[it].second[i].b_pos;

			if (i == reffers[it].second.size() - 1)
				sortFile_B << "," << count_for_1_a[it] << "," << count_for_1_b[it] << "," << count_for_equals[it];
			sortFile_B << endl;
		}

	}

	histFile.close();
	resFile.close();
	sortFile_A.close();
	sortFile_B.close();
}

float DNA_UTILS::dynamicHistogram(vector < pair<vector<DNA_UTILS::positiondata>, vector<DNA_UTILS::positiondata>> > reffers, string filename)
{
	ofstream diffHisto_A, diffHisto_B;
	diffHisto_A.open(".\\Difference Histogram\\" + filename + ".csv");
	diffHisto_B.open(".\\Difference Histogram\\Results_preffered_B\\" + filename + ".csv");
	vector<string> spaces = explode(filename, ' ');
	vector<string> dashes_A = explode(spaces[0], '-');
	vector<string> dashes_B = explode(spaces[3], '-');
	diffHisto_A << "GI,Genus,Specie,GI,Genus,Specie" << endl;
	diffHisto_A << dashes_A[0] << "," << dashes_A[1] << "," << spaces[1] << "," << dashes_B[0] << "," << dashes_B[1] << "," << spaces[4] << endl;
	diffHisto_B << "GI,Genus,Specie,GI,Genus,Specie" << endl;
	diffHisto_B << dashes_A[0] << "," << dashes_A[1] << "," << spaces[1] << "," << dashes_B[0] << "," << dashes_B[1] << "," << spaces[4] << endl;

	diffHisto_A << "Shift,Difference,Frequency,Direct,Neighbor" << endl;
	diffHisto_B << "Shift,Difference,Frequency,Direct,Neighbor" << endl;

	vector<DNA_UTILS::differ> preff_A;
	vector<DNA_UTILS::differ> preff_B;
	DNA_UTILS::differ temp;
	int index;
	int scorchart[4][3];
	int BinVal = 0, maxFreq = 0;
	int sigma, sigmaSquared; 
	for (int it = 0; it < 4; it++)
	{
		btree tree_a = btree();

		for (int i = 0; i < reffers[it].first.size(); i++)
		{
			temp.gap = reffers[it].first[i].a_pos - reffers[it].first[i].b_pos;
			temp.freaqency = 1;
			if (reffers[it].first[i].cast == 0)
			{
				temp.direct = 1;
				temp.neighbor = 0;
			}
			else
			{
				temp.direct = 0;
				temp.neighbor = 1;
			}
			tree_a.insert(temp);
		}

		preff_A = tree_a.printTree();

		btree tree_b = btree();

		
		for (int i = 0; i < reffers[it].second.size(); i++)
		{
			temp.gap = reffers[it].second[i].a_pos - reffers[it].second[i].b_pos;
			temp.freaqency = 1;
			if (reffers[it].second[i].cast == 0)
			{
				temp.direct = 1;
				temp.neighbor = 0;
			}
			else
			{
				temp.direct = 0;
				temp.neighbor = 1;
			}
			tree_b.insert(temp);
		}

		preff_B = tree_b.printTree();
		
		int index;

		if (preff_A.size() > 0)
		{
			index = getMax(preff_A);
			//sigma = getSigma(preff_A);
			//sigmaSquared = getsigmasquared(preff_A);

			if (index == -1)
			{
				scorchart[it][0] = 0;
				scorchart[it][1] = 0;
				scorchart[it][2] = 0;
			}
			else
			{
				if (maxFreq <= preff_A[index].freaqency)
				{
					maxFreq = preff_A[index].freaqency;
					BinVal = abs(preff_A[index].gap);
				}
				scorchart[it][0] = preff_A[index].freaqency;
				scorchart[it][1] = preff_A[index].direct;
				scorchart[it][2] = preff_A[index].neighbor;
			}
		}
		else
		{
			scorchart[it][0] = 0;
			scorchart[it][1] = 0;
			scorchart[it][2] = 0;
		}

		for (int i = 0; i < preff_A.size(); i++)
		{
			if (i == 0)
				diffHisto_A << it;

			diffHisto_A << "," << preff_A[i].gap << "," << preff_A[i].freaqency << "," << preff_A[i].direct << "," << preff_A[i].neighbor << endl;
		}

		for (int i = 0; i < preff_B.size(); i++)
		{
			if (i == 0)
				diffHisto_B << it;

			diffHisto_B << "," << preff_B[i].gap << "," << preff_B[i].freaqency << "," << preff_B[i].direct << "," << preff_B[i].neighbor << endl;
		}

		preff_A.clear();
		preff_B.clear();
	}

	diffHisto_A.close();
	diffHisto_B.close();
	
	int sum_of_freq = scorchart[0][0] + scorchart[1][0] + scorchart[2][0] + scorchart[3][0];
	int sum_of_direct = scorchart[0][1] + scorchart[1][1] + scorchart[2][1] + scorchart[3][1];
	int sum_of_neighbor = scorchart[0][2] + scorchart[1][2] + scorchart[2][2] + scorchart[3][2];

	float result = ((4 * sum_of_direct + 3 * sum_of_neighbor) - sum_of_neighbor - BinVal) / 4.0;
	
	return result;

}

float DNA_UTILS::dynamicHistogram_new(vector < pair<vector<DNA_UTILS::positiondata>, vector<DNA_UTILS::positiondata>> > reffers, string filename)
{
	ofstream diffHisto_A, diffHisto_B;
	diffHisto_A.open(".\\Difference Histogram\\" + filename + ".csv");
	diffHisto_B.open(".\\Difference Histogram\\Results_preffered_B\\" + filename + ".csv");
	vector<string> spaces = explode(filename, ' ');
	vector<string> dashes_A = explode(spaces[0], '-');
	vector<string> dashes_B = explode(spaces[3], '-');
	diffHisto_A << "GI,Genus,Specie,GI,Genus,Specie" << endl;
	diffHisto_A << dashes_A[0] << "," << dashes_A[1] << "," << spaces[1] << "," << dashes_B[0] << "," << dashes_B[1] << "," << spaces[4] << endl;
	diffHisto_B << "GI,Genus,Specie,GI,Genus,Specie" << endl;
	diffHisto_B << dashes_A[0] << "," << dashes_A[1] << "," << spaces[1] << "," << dashes_B[0] << "," << dashes_B[1] << "," << spaces[4] << endl;

	diffHisto_A << "Shift,Difference,Frequency,Direct,Neighbor" << endl;
	diffHisto_B << "Shift,Difference,Frequency,Direct,Neighbor" << endl;

	vector<DNA_UTILS::differ> preff_A;
	vector<DNA_UTILS::differ> preff_B;
	DNA_UTILS::differ temp;
	int index;
	float score[4];
	int BinVal = 0, maxFreq = 0;
	int sigma, sigmaSquared; 
	for (int it = 0; it < 4; it++)
	{
		btree tree_a = btree();

		for (int i = 0; i < reffers[it].first.size(); i++)
		{
			temp.gap = reffers[it].first[i].a_pos - reffers[it].first[i].b_pos;
			temp.freaqency = 1;
			if (reffers[it].first[i].cast == 0)
			{
				temp.direct = 1;
				temp.neighbor = 0;
			}
			else
			{
				temp.direct = 0;
				temp.neighbor = 1;
			}
			tree_a.insert(temp);
		}

		preff_A = tree_a.printTree();

		btree tree_b = btree();

		
		for (int i = 0; i < reffers[it].second.size(); i++)
		{
			temp.gap = reffers[it].second[i].a_pos - reffers[it].second[i].b_pos;
			temp.freaqency = 1;
			if (reffers[it].second[i].cast == 0)
			{
				temp.direct = 1;
				temp.neighbor = 0;
			}
			else
			{
				temp.direct = 0;
				temp.neighbor = 1;
			}
			tree_b.insert(temp);
		}

		preff_B = tree_b.printTree();
			
		score[it] =  getNewScore(preff_A);
			
		for (int i = 0; i < preff_A.size(); i++)
		{
			if (i == 0)
				diffHisto_A << it;

			diffHisto_A << "," << preff_A[i].gap << "," << preff_A[i].freaqency << "," << preff_A[i].direct << "," << preff_A[i].neighbor << endl;
		}

		for (int i = 0; i < preff_B.size(); i++)
		{
			if (i == 0)
				diffHisto_B << it;

			diffHisto_B << "," << preff_B[i].gap << "," << preff_B[i].freaqency << "," << preff_B[i].direct << "," << preff_B[i].neighbor << endl;
		}

		preff_A.clear();
		preff_B.clear();
	}

	diffHisto_A.close();
	diffHisto_B.close();

	float result = takemax(score);
	
	return result;

}



/*
pair<vector<int>,vector<int>> DNA_UTILS::newaddition(vector < pair<vector<DNA_UTILS::positiondata>, vector<DNA_UTILS::positiondata>> > reffers, string filename)
{
	

	vector<DNA_UTILS::differ> preff_A;
	DNA_UTILS::differ temp;
	vector<int> sigmas, squaredSigmas;
	for (int it = 0; it < 4; it++)
	{
		btree tree_a = btree();

		for (int i = 0; i < reffers[it].first.size(); i++)
		{
			temp.gap = reffers[it].first[i].a_pos - reffers[it].first[i].b_pos;
			temp.freaqency = 1;
			if (reffers[it].first[i].cast == 0)
			{
				temp.direct = 1;
				temp.neighbor = 0;
			}
			else
			{
				temp.direct = 0;
				temp.neighbor = 1;
			}
			tree_a.insert(temp);
		}
		preff_A = tree_a.printTree();

		if (preff_A.size() > 0)
		{
			sigmas.push_back(getSigma(preff_A));
			squaredSigmas.push_back(getsigmasquared(preff_A));
			
		}
		else
		{
			cout<< "not enough files" << endl;
		}

		preff_A.clear();
	}

	return make_pair(sigmas, squaredSigmas);

}

vector<float>DNA_UTILS::newaddition2(vector < pair<vector<DNA_UTILS::positiondata>, vector<DNA_UTILS::positiondata>> > reffers, vector<string> filename_a, vector<string> filename_b)
{
	
	vector<DNA_UTILS::differ> preff_A;
	DNA_UTILS::differ temp;
	vector<float> score;

	for (int it = 0; it < 4; it++)
	{
		btree tree_a = btree();

		for (int i = 0; i < reffers[it].first.size(); i++)
		{
			temp.gap = reffers[it].first[i].a_pos - reffers[it].first[i].b_pos;
			temp.freaqency = 1;
			if (reffers[it].first[i].cast == 0)
			{
				temp.direct = 1;
				temp.neighbor = 0;
			}
			else
			{
				temp.direct = 0;
				temp.neighbor = 1;
			}
			tree_a.insert(temp);
		}
		preff_A = tree_a.printTree();

		if (preff_A.size() > 0)
		{
			//sigmas.push_back(getSigma(preff_A));
			score.push_back(getNewScore(preff_A));
			
		}
		else
		{
			cout<< "not enough files" << endl;
		}

		preff_A.clear();
	}
	// check whethr its intergenus, samegenus or same spacie
	//cout<< filename_a[0] << " " <<filename_a[1] << " " <<filename_b[0] << " " <<filename_b[1] <<endl;
	//getchar();
	if (filename_a[0]==filename_b[0] && filename_a[1]==filename_b[1]){
		//samespacie
		score.push_back(1);
	}else if(filename_a[0]==filename_b[0]){
		//same genus
		score.push_back(2);
	}else{
		//intergenus 
		score.push_back(3);
	}

	return score;

}

void DNA_UTILS::printNewScores(vector<vector<float>> ScoringMatrix, vector<vector<float>> scorings, string line1, string line2, string line3)
{
	ofstream ScoreFile,ScoreFileSigma2;
	ScoreFile.open("Scoring Matrix new.csv");
	int size = ScoringMatrix.size();
	for(int it = 0;it <4 ; it++){
		ScoreFile << line1 << endl;
		ScoreFile << line2 << endl;
		ScoreFile << line3 << endl;
		
		vector <string> col_1 = explode(line1, ',');
		vector <string> col_2 = explode(line2, ',');
		vector <string> col_3 = explode(line3, ',');
		
		int k = 0;	
		
		for (int i = 0; i < size; i++)
		{
			int j;
			ScoreFile << col_1[i] << ",";
			ScoreFile << col_2[i] << ",";
			ScoreFile << col_3[i] << ",";

			//cout<<ScoringMatrix.size() << endl;

			for (j = 0; j <= i; j++)
			{
				ScoreFile << "-,";
			}

			for (; j <size+1; j++)
			{

				ScoreFile << scorings[k][it] << ",";
				k++;
			}
			ScoreFile << endl;
		}
	
		ScoreFile << endl;
		ScoreFile << endl;
	
	}

	ScoreFile << line1 << endl;
	ScoreFile << line2 << endl;
	ScoreFile << line3 << endl;

	vector <string> col_1 = explode(line1, ',');
	vector <string> col_2 = explode(line2, ',');
	vector <string> col_3 = explode(line3, ',');

		int k = 0;	
		//int size = ScoringMatrix.size();
		for (int i = 0; i < size; i++)
		{
			int j;
			ScoreFile << col_1[i] << ",";
			ScoreFile << col_2[i] << ",";
			ScoreFile << col_3[i] << ",";

			//cout<<ScoringMatrix.size() << endl;

			for (j = 0; j <= i; j++)
			{
				ScoreFile << "-,";
			}

			for (; j <size+1; j++)
			{

				ScoreFile << takemax(scorings[k][0],scorings[k][1],scorings[k][2],scorings[k][3]) << ",";
				k++;
			}
			ScoreFile << endl;
		}

	ScoreFile.close();

	ofstream Specie_cat;
	Specie_cat.open("catagorizing.csv");
	int values = scorings.size();
		//ScoreFile << line1 << endl;
		//ScoreFile << line2 << endl;
		//ScoreFile << line3 << endl;
	
		//vector <string> col_1 = explode(line1, ',');
		//vector <string> col_2 = explode(line2, ',');
		//vector <string> col_3 = explode(line3, ',');
			
		//int size = ScoringMatrix.size();
	Specie_cat << "same spacie values" << endl;
	for (int i = 0; i < values; i++)
		{
			if(scorings[i][4]==1){
				Specie_cat << takemax(scorings[i][0],scorings[i][1],scorings[i][2],scorings[i][3])<<",";
			}

		}
		Specie_cat << endl;
		Specie_cat << "same genus values" << endl;
		for (int i = 0; i < values; i++)
		{
			if(scorings[i][4]==2){
			Specie_cat << takemax(scorings[i][0],scorings[i][1],scorings[i][2],scorings[i][3])<<",";
			}

		}
		
		Specie_cat << endl;
		Specie_cat << "inter genus values" << endl;
		for (int i = 0; i < values; i++)
		{
			if(scorings[i][4]==3){
			Specie_cat << takemax(scorings[i][0],scorings[i][1],scorings[i][2],scorings[i][3])<<",";
			}

		}

		Specie_cat << endl;
		Specie_cat << endl;
		Specie_cat.close();
}


	/*ScoreFileSigma2.open("Scoring Matrix new2.csv");
	for(int it = 0;it <4 ; it++){
		ScoreFileSigma2 << line1 << endl;
		ScoreFileSigma2 << line2 << endl;
		ScoreFileSigma2 << line3 << endl;
	
		vector <string> col_1 = explode(line1, ',');
		vector <string> col_2 = explode(line2, ',');
		vector <string> col_3 = explode(line3, ',');
		
		int k = 0;	
		int size = ScoringMatrix.size();
		for (int i = 0; i < size; i++)
		{
			int j;
			ScoreFileSigma2 << col_1[i] << ",";
			ScoreFileSigma2 << col_2[i] << ",";
			ScoreFileSigma2 << col_3[i] << ",";

			//cout<<ScoringMatrix.size() << endl;

			for (j = 0; j <= i; j++)
			{
				ScoreFileSigma2 << "-,";
			}

			for (; j <size+1; j++)
			{

				ScoreFileSigma2 << scorings[k].second[it] << ",";
				k++;
			}
			ScoreFileSigma2 << endl;
		}
	
		ScoreFileSigma2 << endl;
		ScoreFileSigma2 << endl;
	
	}
	ScoreFileSigma2.close();
	*/
	


float DNA_UTILS::takemax(float a[4]){
	
	float max = a[0];
	for(int i = 1;i< 4; i++){
		if(max<a[i]){
			max = a[i];
		}

	}
	return max;
}

/*

int DNA_UTILS::getMemberIndex(vector<DNA_UTILS::differ> base, int data)
{
	for (int i = 0; i < base.size(); i++)
	{
		if (base[i].gap == data)
			return i;
	}

	return -1;
} 
*/


int DNA_UTILS::getMax(vector<DNA_UTILS::differ> arg)
{
	int max = arg[0].freaqency;
	int i, maxIndex;
	int cond = 0;
	for (i = 0; i < arg.size(); i++)
	{
		if (max <= arg[i].freaqency)
		{
			if (max == arg[i].freaqency)
				cond = 1;
			max = arg[i].freaqency;
			maxIndex = i;
		}
	}

	if (cond == 1 && max < 5)
		maxIndex = -1;

	return maxIndex;

}

// use this
float DNA_UTILS::getNewScore(vector<DNA_UTILS::differ> arg)
{
	float Wg = 1 , Wi = 1, WM = 5; 
	float total = 0 ;
	int freqTotal = 0;
	for (int i = 0; i < arg.size(); i++)
	{
		total = total + arg[i].freaqency*(WM*(6*arg[i].direct + 1.5*arg[i].neighbor)-Wi*arg[i].neighbor-Wg*abs(arg[i].gap));
		freqTotal = freqTotal + arg[i].freaqency;
		
	}
	float score = 1.0*total/freqTotal;
	return score;
}

/*
int DNA_UTILS::getSigma(vector<DNA_UTILS::differ> arg)
{
	
	int sigma = 0 ;
	for (int i = 0; i < arg.size(); i++)
	{
		
		sigma = sigma + arg[i].freaqency * abs(arg[i].gap); 
	}
	return sigma;

}

int DNA_UTILS::getsigmasquared(vector<DNA_UTILS::differ> arg)
{
	
	int sigmaSquared = 0 ;
	for (int i = 0; i < arg.size(); i++)
	{	
		sigmaSquared = sigmaSquared + arg[i].freaqency * (abs(arg[i].gap) * abs(arg[i].gap)); 
	}
	return sigmaSquared;

}
*/

string DNA_UTILS::bitsSelector(string sample)
{
	if (sample == "11")
		return "100100";
	else if (sample == "10")
		return "110100";
	else if (sample == "01")
		return "111000";
	else
		return "111001";
}

string DNA_UTILS::encode(char sample)
{
	if (sample == 'A')
		return "11";
	else if (sample == 'C')
		return "01";
	else if (sample == 'G')
		return "10";
	else
		return "00";
}

char DNA_UTILS::decode(string twobits_str)
{
	if (twobits_str == "11")
		return 'A';
	else if (twobits_str == "01")
		return 'C';
	else if (twobits_str == "10")
		return 'G';
	else
		return 'T';
}

string DNA_UTILS::decode_eight(bitset<8> eightbits)
{
	string eightbits_str = eightbits.to_string<char, char_traits<char>, allocator<char> >();
	string buffer = "";
	for (int i = 0; i < eightbits_str.length(); i += 2)
	{
		buffer += decode(eightbits_str.substr(i, 2));
	}
	return buffer;
}

vector<string> DNA_UTILS::explode(const string& str, const char& ch)
{
	string next;
	vector<string> result;
	// For each character in the string
	for (string::const_iterator it = str.begin(); it != str.end(); it++)
	{
		// If we've hit the terminal character
		if (*it == ch)
		{
			// If we have some characters accumulated
			if (!next.empty())
			{
				// Add them to the result vector
				result.push_back(next);
				next.clear();
			}
		}
		else
		{
			// Accumulate the next character into the sequence
			next += *it;
		}
	}

	if (!next.empty())
		result.push_back(next);

	return result;
}

int DNA_UTILS::getMaximum(int x, int y)
{
	if (x > y)
	{
		return x;
	}
	else
	{
		return y;
	}
}