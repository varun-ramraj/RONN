/*=========================================================
**
** This program is designed for running RONNv3.2 for detecting
** disordered residues in proteins.
** Date:	March 2012
** Designer:	Zheng Rong Yang (Exeter University)
** Programmer:	Varun Ramraj (University of Oxford)
** Contact:	varun@strubi.ox.ac.uk
**
** Modified by Shyam Saladi (saladi@caltech.edu, California Institute of Technology), 01-09-15
**
==========================================================*/


#include "callBBF.h"

using namespace std;

//@@@@@ the path can be changed
string myPath="/usr/local/RONNv3_2";

int WINSIZE = 19;

int main(int argc,char **argv)
{
	//number of models
	int nM=10;
	
	stringstream *convert;
	FILE *fp;

	//$VR$: VERY IMPORTANT. CONTROLS COST FUNCTION
	double disorder_weight=0.53;
	
	string fModel, fPDF;
	int r, m, nR;
	
	double mean;
	vector< vector<double> > X;

	string query;

	//$VR$: Debug output
	//cout << "The cost coefficient or disorder weight is: " << disorder_weight << endl;

	string fastaheader;	
	string line;
    // print the FASTA header
	getline(std::cin, line);

	fastaheader = line;

	//cout << line;
	//cout << "\n";
	
	// get the sequence
    while (cin >> line)
    {
        query = query + line;
    }
	nR=query.size();

	if (nR < WINSIZE)
	{
	    string qscore = ">";
	    cout << fastaheader << endl;

	    for (int qq = 0; qq < nR; qq++)
	    {
		qscore += "?";
	    }

	    cout << qscore << endl;

	    cout << ">$" << endl;

	    return 0;

	}

	for(m=0; m<nM; m++)
	{
		//sprintf(fModel,"%s/c%d/model.rec",myPath,m);
		convert = new stringstream();
		(*convert) << m;
		string mString = convert->str();
	
		fModel = myPath + "/c" + mString + "/model.rec";
		
	   	//TODO: get rid of FILE* stuff
		if(!(fp=fopen(fModel.c_str(),"r"))) 
		{ 
			cout << "Can't find " << fModel << endl;
			return -1;
		}
		
		fclose(fp);

		//sprintf(fPDF,"%s/c%d/pdfs.rec",myPath,m);
		fPDF = myPath + "/c" + mString + "/pdfs.rec";
		
		//TODO: get rid of FILE* stuff
		if(!(fp=fopen(fPDF.c_str(),"r"))) 
		{ 
			cout << "Can't find " << fPDF << endl; 
			return -1;
		}
		fclose(fp);

        vector<double> scores;
		int retval = callBBF_driver(query, fModel, fPDF, disorder_weight, scores);

		//if (retval == 0)  
		//  printf("................Completing testing ............! Model: %d\n", m);

        for(r=0;r<nR;r++)
		{
			X.push_back(vector<double>()); 
			X[r].push_back(scores[r]);
		}
	}


	//print out headers
	
	string seqstring = ">";
	string headstring = ">";

	for (r = 0; r < nR; r++)
	{
	    mean = 0;

	    for (m = 0; m < nM; m++)
	    {
		mean += X[r][m];

	    }

	    mean /= (double)nM;

	    //4 conditions
	    if (mean < 0.4) {  headstring += " "; }
	    else if (mean >= 0.4 && mean < 0.5) {  headstring += "-"; }
	    else if (mean >= 0.5 && mean < 0.6) {  headstring += "="; }
	    else if (mean >= 0.6) {  headstring += "#";  }

	    seqstring += query[r];
	}	
	
	headstring += "\n>$";

	cout << fastaheader << endl;
	cout << seqstring << endl;	
	cout << headstring << endl;	
	
	
	
	


	//print out per residue score	
	for(r=0;r<nR;r++)
	{
		mean=0;
		
		for(m=0; m<nM; m++)
			mean+=X[r][m];
		
		mean/=(double)nM;
		
		printf("%c\t%f\n",query[r],mean);
		
		
	}
	
	delete convert;
	
	return 0;
}
