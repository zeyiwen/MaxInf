#ifndef EEP_H
#define EEP_H
#include "Utility.h"
#include "EE.h"
using namespace EEP_SPACE;

class EEPSolver
{
public:
	EEPSolver(string fileM, string fileF, string fileC, int ma, int fa, int ca);
	~EEPSolver(){};
	void Solve(); //use to solve the query
	void GetData(string fileM, string fileF, string fileC, int ma, int fa, int ca); //read data from files and build Rtrees
	void PrintDuration();
	void PrintDistComputation();
	void PrintAnswerSet();
	void PrintIO();
private:
	//three Rtrees
	aRTree<int, float, 2> treeM;
	RTree_EEP<int, float, 2> treeF;
	RTree_EEP<int, float, 2> treeC;
	//number of objects
	int nMa;
	int nFa;
	int nCa;
	//CPU time
	double dDuration;
	//distance computation counter
	long long nTotalCounter;
	long long nExactCounter;
	//IO counter;
	long long nAccessM;
	long long nAccessF;
	long long nAccessC;
	//answer set
	list<EntryC*> answerSet;
};
EEPSolver::EEPSolver(string fileM, string fileF, string fileC, int ma, int fa, int ca)
{
	nMa = ma;
	nFa = fa;
	nCa = ca;
	GetData(fileM, fileF, fileC, ma, fa, ca);
}
void EEPSolver::PrintIO()
{

}
void EEPSolver::PrintAnswerSet()
{
	EntryC *pEnC;
	int nOutput = 1;
	for(int i = 0; i < K; i++)
	{
		pEnC = answerSet.front();
		cout << nOutput << ". minInf: " << pEnC->minInf << " | " << "maxInf: " << pEnC->maxInf << '\t' 
			<<"[ "<< pEnC->pC->m_rect.m_min[0] << " , " << pEnC->pC->m_rect.m_min[1] << " ]"<< '\n';
		nOutput++;
		answerSet.pop_front();
	}
	//cout << "C " << nCa << " F " << nFa << " M " << nMa << endl;
}
void EEPSolver::PrintDuration()
{
	cout << "EEP CPU time :" << dDuration << " ms" << endl;
}
void EEPSolver::PrintDistComputation()
{
	cout << "Total Computation: " << counter << endl;
	cout << "Pruning Computation: " << counter - ExactCounter << endl;
	cout << "Exact Computation: " << ExactCounter <<  endl;
}
void EEPSolver::Solve()
{
	start_timer();
	EEP(treeM.m_root, treeF.m_root, treeC.m_root, K, nMa, nFa, nCa);
	dDuration = stop_timer(0);
	answerSet = result_c2;
	result_c2.clear();
}
void EEPSolver::GetData(string fileM, string fileF, string fileC, int ma, int fa, int ca)
{
	int nSizeofNode = 0;
	nSizeofNode = sizeof(Node);
	int nSizeofaNode = 0;
	nSizeofaNode = sizeof(aNode);
	const char * fn;

	static Point dataM[MAX_SIZE_DATAM];
	static Point dataF[MAX_SIZE_DATAF];
	static Point dataC[MAX_SIZE_DATAC];

	//read file
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
		count++;
	filestr.close();
	filestr.clear();

	float arrary_Max[2] = {0,0};
	int dataID = 0;
	arrary_Max[0] = 0, arrary_Max[1] = 0;
	dataID = 0;

	for(int i = 0;i < fa;i++, dataID++)
	{
		arrary_Max[0] = dataF[i].x;
		arrary_Max[1] = dataF[i].y;
		treeF.Insert(arrary_Max,arrary_Max,dataID);
	}

	for(int i = 0;i < ma;i++, dataID++)
	{		
		arrary_Max[0] = dataM[i].x;
		arrary_Max[1] = dataM[i].y;
		treeM.Insert(arrary_Max,arrary_Max,dataID);
	}

	arrary_Max[0] = 0, arrary_Max[1] = 0;
	dataID = 0;
	for(int i = 0;i < ca;i++, dataID++)
	{
		arrary_Max[0] = dataC[i].x;
		arrary_Max[1] = dataC[i].y;
		treeC.Insert(arrary_Max,arrary_Max,dataID);
	}
}


#endif EEP_H