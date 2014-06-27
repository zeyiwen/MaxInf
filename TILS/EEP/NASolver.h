#ifndef VORNLC_H
#define VORNLC_H
#include "Utility.h"

struct Facility
{
	float x;
	float y;
};

struct Client
{
	float x;
	float y;
	float nn;
};

struct Candidate
{
	float x;
	float y;
	int count;
};

class NASolver
{
private:
	//number of objects
	int nMa;
	int nFa;
	int nCa;
	//array for storing points from files
	Client dataM[MAX_SIZE_DATAM];
	Facility dataF[MAX_SIZE_DATAF];
	Candidate dataC[MAX_SIZE_DATAC];
	//timer
	double dDuration;
public:
	NASolver(string fileM, string fileF, string fileC, int ma, int fa, int ca)
	{
		nMa = ma;
		nFa = fa;
		nCa = ca;
		GetData(fileM, fileF, fileC);
	}
	void Solve(); //use to solve the query
	void GetData(string fileM, string fileF, string fileC); //read data from files
	void PrintDuration();
	void PrintAnswerSet();
	void PrintComputation();
};
void NASolver::PrintComputation()
{
	unsigned long long nTotalComputation;
	nTotalComputation = nMa * nFa;
	nTotalComputation += nMa * nCa;
	cout << "Total Computation: " << nTotalComputation << endl;
}
void NASolver::PrintAnswerSet()
{
	for(int k = 0; k < K; k++)
	{
		cout<< k + 1 << ". \t" << dataC[k].count<<'\t' <<"[ "<< dataC[k].x << " , " << dataC[k].y << " ]" << '\n';
	}
	//cout << "C " << nCa << " F " << nFa << " M " << nMa << endl;
}
void NASolver::PrintDuration()
{
	cout << "NA CPU time: " << dDuration << " ms" << endl;
}
void NASolver::GetData(string fileM, string fileF, string fileC)
{
	const char * fn;
	fstream filestr;
	fn = fileM.c_str();
	filestr.open(fn, fstream::in);
	int id, count = 0;
	char dot;
	while(filestr>>dataM[count].x>>dot>>dataM[count].y>>dot>>id)
		count++;
	filestr.close();
	filestr.clear();

	count = 0;
	fn = fileF.c_str();
	filestr.open(fn, fstream::in);
	while(filestr>>dataF[count].x>>dot>>dataF[count].y>>dot>>id)
		count++;
	filestr.close();
	filestr.clear();

	count = 0;
	fn = fileC.c_str();
	filestr.open(fn, fstream::in);
	while(filestr>>dataC[count].x>>dot>>dataC[count].y>>dot>>id)
	{
		dataC[count].count = 0;
		count++;
	}
	filestr.close();
	filestr.clear();
}
void NASolver::Solve()
{
	start_timer();
	int nMpoint = nMa;
	int nFpoint = nFa;
	int nCpoint = nCa;

	double distTemp;
	double curDist;

	double distTemp2;
	
	for(int i=0; i<nMpoint; i++)
	{
		//compute nearest facility distance for each m
		curDist = (dataM[i].x - dataF[0].x)*(dataM[i].x - dataF[0].x) + 
			(dataM[i].y - dataF[0].y)*(dataM[i].y - dataF[0].y);
		for(int j=1; j<nFpoint; j++)
		{
			distTemp =(dataM[i].x - dataF[j].x)*(dataM[i].x - dataF[j].x) + 
				(dataM[i].y - dataF[j].y)*(dataM[i].y - dataF[j].y);
			if(distTemp < curDist)
				curDist = distTemp;
		}
		//compute influence value for each candidate point
		for(int j=0; j<nCpoint; j++)
		{
			distTemp2 =((dataM[i].x - dataC[j].x)*(dataM[i].x - dataC[j].x) + 
				(dataM[i].y - dataC[j].y)*(dataM[i].y - dataC[j].y));
			if (distTemp2 <= curDist)
				dataC[j].count++;
		}
	}

	float fTemp[2];
	int nTemp;
	int nIndex;
	//sort candidate location by influence value
	for(int n = 0;n < nCpoint; n++)
	{
		nTemp = dataC[n].count;
		nIndex = n;
		for(int l = n + 1;l < nCpoint; l++)
		{
			if(nTemp < dataC[l].count)
			{
				nTemp = dataC[l].count;
				nIndex = l;
			}
		}

		fTemp[0] = dataC[n].x;
		fTemp[1] = dataC[n].y;
		nTemp = dataC[n].count;
		dataC[n].x = dataC[nIndex].x;
		dataC[n].y = dataC[nIndex].y;
		dataC[n].count = dataC[nIndex].count;
		dataC[nIndex].x = fTemp[0];
		dataC[nIndex].y = fTemp[1];
		dataC[nIndex].count = nTemp;
	}
	dDuration = stop_timer(0);
}
#endif VORNLC_H