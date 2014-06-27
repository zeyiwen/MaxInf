#ifndef SOLVER_H
#define SOLVER_H
//common header
#include "RTree.h"
#include <math.h>
#include <time.h>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "constant.h"
using namespace std;

//Constant
#define PARAN		7 //3
#define RANGE		1024
#define DIMENSION	2


//struct
struct Point
{
	double x;
	double y;
};

///----------------------------------------
/// Abstract Class Solver
/// general utilities and common interface
///----------------------------------------

class Solver
{
public:
	Solver(){};
	virtual ~Solver(){};
	virtual void Solve() = 0;
	virtual void Print() = 0;

protected:
	double GetDistBtwnS(const Point pa, const Point pb)
	{return (pa.x - pb.x)*(pa.x - pb.x) + (pa.y - pb.y)*(pa.y - pb.y);};

	double GetDistBtwn(const Point pa, const Point pb)
	{return sqrt(GetDistBtwnS(pa, pb));};

	double GetDistBtwnS(const double ax, const double ay, const double bx, const double by)
	{return (ax-bx)*(ax-bx) + (ay-by)*(ay-by);};

	double GetDistBtwn(const double ax, const double ay, const double bx, const double by)
	{return sqrt(GetDistBtwnS(ax, ay, bx, by));};

	double GetArea(const RTree<unsigned int, double, DIMENSION, double>::Rect rect)
	{return (rect.m_max[0] - rect.m_min[0]) * (rect.m_max[1] - rect.m_min[1]);};

	double GetAreaRoot(const RTree<unsigned int, double, DIMENSION, double>::Rect rect)
	{return sqrt((rect.m_max[0] - rect.m_min[0])) * sqrt((rect.m_max[1] - rect.m_min[1]));};
};

///----------------------------------------
/// Abstract Class MemorySolver
///----------------------------------------
class MemorySolver: public Solver
{
public:
	MemorySolver(string fileC, string fileF, string fileM)
	{
		ReadFile(fileC, dataC);
		ReadFile(fileF, dataF);
		ReadFile(fileM, dataM);
		influence.resize(FILE_C_SIZE);
		for(int c = 0; c < FILE_C_SIZE; c++) influence.push_back(0);
	}
	virtual ~MemorySolver(){};
	void Print();

protected:
	vector<int> influence;
	vector<Point> dataC;
	vector<Point> dataF;
	vector<Point> dataM;
	long during;
	long start;
	long end;

	void ReadFile(string, vector<Point>&);
	void StartTiming(){start = clock();};
	void EndTiming(){end = clock(); during = end - start;};
};

void MemorySolver::ReadFile(string filename, vector<Point>& data)
{
	const char* fn;
	fstream filestr;
	fn = filename.c_str();
	char dot;
	int id;
	Point eachPoint;
	filestr.open(fn, fstream::in);
	while(filestr>>eachPoint.x>>dot>>eachPoint.y>>dot>>id)
	{
		data.push_back(eachPoint);
	}
	filestr.close();
}

void MemorySolver::Print()
{
	multimap<int, int> index;
	for(int c = 0; c < (int)influence.size(); c++)
		index.insert(pair<int, int>(influence[c], c));
	multimap<int, int>::iterator itR = index.end();
	if(influence.size() < K)
	{
		for(int i = 0; i < (int)influence.size(); i++)
		{
			itR--;
			cout<<i+1<<'\t'<<itR->first;
			cout<<'\t'<<"[ "<<dataC[itR->second].x<<" , "<<dataC[itR->second].y<<" ]"<<endl;
		}
		cout<<"All rests are zeros"<<endl;
	}
	for(int i = 0; i < K; i++)
	{
		itR--;
		cout<<i+1<<'\t'<<itR->first;
		cout<<'\t'<<"[ "<<dataC[itR->second].x<<" , "<<dataC[itR->second].y<<" ]"<<endl;
	}
	cout<<"C"<<FILE_C_SIZE<<" F"<<FILE_F_SIZE<<" M"<<FILE_M_SIZE<<" in "<<during<<" ms."<<endl;
}

///----------------------------------------
/// Abstract Class ExternalSolver
///----------------------------------------
class ExternalSolver: public Solver
{
public:
	ExternalSolver(string fileM, string fileF, string fileC, int ma, int fa, int ca)
	{
		ReadRTree(fileC, treeC);
		ReadRTree(fileF, treeF);
		ReadRTree(fileM, treeM);
		countPA = 0;
		countIO = 0;
		paraC = ca;
		paraF = fa;
		paraM = ma;
	}
	void Print();
protected:
	RTree<unsigned int, double, DIMENSION, double> treeC;
	RTree<unsigned int, double, DIMENSION, double> treeF;
	RTree<unsigned int, double, DIMENSION, double> treeM;
	int countPA;
	int countIO;
	int paraC;
	int paraF;
	int paraM;
	multimap<int, Point> answer;

	void ReadRTree(const string, RTree<unsigned int, double, DIMENSION, double>&);

	double GetMinMaxDistBtwn(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	bool GetMinMaxMidBtwn(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect,
		Point*);
	bool IsContained(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	bool IsIntersect(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	
	//debug finished
	int CheckOrientation(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	void GetTypicalRatio(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const int,
		double *);
	bool GetIntersect(const Point, const double, const Point, const double, Point &);
	bool CheckPointAndLine(const Point, const double, const Point);

public:
	double GetMinDistBtwn(const Point, const RTree<unsigned int, double, DIMENSION, double>::Rect);
	double GetMinMaxDistBtwn(const Point, const RTree<unsigned int, double, DIMENSION, double>::Rect);

protected:
	//debug finished
	void SplitRect(vector<RTree<unsigned int,double,DIMENSION,double>::Rect >*, 
				RTree<unsigned int,double,DIMENSION,double>::Rect *, 
				vector<RTree<unsigned int, double, DIMENSION, double>::Rect> &);

	//debug finished
	void SplitHelper(const RTree<unsigned int, double, DIMENSION, double>::Rect *,
		const RTree<unsigned int, double, DIMENSION, double>::Rect *,
		vector<RTree<unsigned int, double, DIMENSION, double>::Rect> &);
};

void ExternalSolver::ReadRTree(const string file, 
							   RTree<unsigned int,double,2,double>& tree)
{
	fstream filestr;
	char dot;
	double eachPoint[DIMENSION];
	int id;
	const char * fn = file.c_str();
	filestr.open(fn, fstream::in);
	while(filestr>>eachPoint[0]>>dot>>eachPoint[1]>>dot>>id)
		tree.Insert(eachPoint, eachPoint, id);
	filestr.close();
}

void ExternalSolver::Print()
{
	multimap<int, Point>::iterator itA = answer.end();
	for(int i = 0; i < K; i++)
	{
		itA--;
		cout<<i+1<<'\t'<<itA->first<<'\t';
		cout<<"[ "<<itA->second.x<<" , "<<itA->second.y<<" ]."<<endl;
	}
}


// returns -1 if rb is contained in ra, 0 if ra and rb intersect, and minMax(ra, rb) dist otherwise
double ExternalSolver::GetMinMaxDistBtwn(const RTree<unsigned int, double, DIMENSION, double>::Rect ra,
					  const RTree<unsigned int, double, DIMENSION, double>::Rect rb)
{
	if(IsContained(rb, ra))
		return -1.0;
	if(IsIntersect(rb, ra))
		return 0.0;
	bool xintersect = !((rb.m_min[0] > ra.m_max[0])||(rb.m_max[0] < ra.m_min[0]));
	bool yintersect = !((rb.m_min[1] > ra.m_max[1])||(rb.m_max[1] < ra.m_min[1]));
	//anti-clockwise
	if(rb.m_max[0] < ra.m_min[0] && rb.m_max[1] < ra.m_min[1])
		return GetDistBtwn(ra.m_min[0], ra.m_min[1], rb.m_min[0], rb.m_min[1]);
	else if(xintersect && rb.m_max[1] < ra.m_min[1])
		return (ra.m_min[1] - rb.m_max[1]);
	else if(rb.m_min[0] > ra.m_max[0] && rb.m_max[1] < ra.m_min[1])
		return GetDistBtwn(ra.m_max[0], ra.m_min[1], rb.m_max[0], rb.m_min[1]);
	else if(yintersect && rb.m_min[0] > ra.m_max[0])
		return (rb.m_min[0] - ra.m_max[0]);
	else if(rb.m_min[0] > ra.m_max[0] && rb.m_min[1] > ra.m_max[1])
		return GetDistBtwn(ra.m_max[0], ra.m_max[1], rb.m_max[0], rb.m_max[1]);
	else if(xintersect && rb.m_min[1] > ra.m_max[1])
		return (rb.m_min[1] - ra.m_max[1]);
	else if(rb.m_max[0] < ra.m_min[0] && rb.m_min[1] > ra.m_max[1])
		return GetDistBtwn(ra.m_min[0], ra.m_max[1], rb.m_min[0], rb.m_max[1]);
	else if(yintersect && rb.m_max[0] < ra.m_min[0])
		return (ra.m_min[0] - rb.m_max[0]);
	else return -2.0;//error
}

bool ExternalSolver::IsContained(const RTree<unsigned int, double, DIMENSION, double>::Rect ra,
								 const RTree<unsigned int, double, DIMENSION, double>::Rect rb)
{
	if(ra.m_min[0] > rb.m_min[0]
	&& ra.m_max[0] < rb.m_max[0]
	&& ra.m_min[1] > rb.m_min[1]
	&& ra.m_max[1] < rb.m_max[1])
		return true;
	else return false;
}

bool ExternalSolver::IsIntersect(const RTree<unsigned int, double, DIMENSION, double>::Rect ra,
								 const RTree<unsigned int, double, DIMENSION, double>::Rect rb)
{
	RTree<unsigned int, double, DIMENSION,double>::Rect intersectRect;
	intersectRect.m_min[0] = max(ra.m_min[0], rb.m_min[0]);
	intersectRect.m_min[1] = max(ra.m_min[1], rb.m_min[1]);
	intersectRect.m_max[0] = min(ra.m_max[0], rb.m_max[0]);
	intersectRect.m_max[1] = min(ra.m_max[1], rb.m_max[1]);

	if(intersectRect.m_min[0] > intersectRect.m_max[0]
	|| intersectRect.m_min[1] > intersectRect.m_max[1])
		return false;
	else if(IsContained(ra, rb) || IsContained(rb, ra))
		return false;
	else
		return true;
}

/// ------------------------------------------------------------------
/// function CheckOrientation
/// 0 -- intersects, contained, or contains;
/// 1 -- SW; 2 -- S; 3 -- SE; 4 -- E; 5 -- NE; 6 -- N; 7 -- NW; 8 -- W
/// ------------------------------------------------------------------
int ExternalSolver::CheckOrientation(const RTree<unsigned int, double, DIMENSION, double>::Rect rTarget,
									 const RTree<unsigned int, double, DIMENSION, double>::Rect rReference)
{
	if(IsIntersect(rTarget, rReference) || IsContained(rTarget, rReference) || IsContained(rReference, rTarget))
	{
		return 0;
	}
	bool xintersect = !((rTarget.m_min[0] > rReference.m_max[0])||(rTarget.m_max[0] < rReference.m_min[0]));
	bool yintersect = !((rTarget.m_min[1] > rReference.m_max[1])||(rTarget.m_max[1] < rReference.m_min[1]));
	if(rTarget.m_max[0] < rReference.m_min[0] && rTarget.m_max[1] < rReference.m_min[1] )
		return 1;
	if(xintersect && rTarget.m_max[1] < rReference.m_min[1])
		return 2;
	if(rTarget.m_min[0] > rReference.m_max[0] && rTarget.m_max[1] < rReference.m_min[1])
		return 3;
	if(yintersect && rTarget.m_min[0] > rReference.m_max[0])
		return 4;
	if(rTarget.m_min[0] > rReference.m_max[0] && rTarget.m_min[1] > rReference.m_max[1])
		return 5;
	if(xintersect && rTarget.m_min[1] > rReference.m_max[1])
		return 6;
	if(rTarget.m_max[0] < rReference.m_min[0] && rTarget.m_min[1] > rReference.m_max[1])
		return 7;
	if(yintersect && rTarget.m_max[0] < rReference.m_min[0])
		return 8;
	else return 0;
}

///-----------------------------------------------------
///if intersects, contained, or contains, returns false
///otherwise return true and p as the mid Point
///-----------------------------------------------------
bool ExternalSolver::GetMinMaxMidBtwn(const RTree<unsigned int, double, DIMENSION, double>::Rect rTarget,
									   const RTree<unsigned int, double, DIMENSION, double>::Rect rReference,
									   Point* p)
{
	if(IsContained(rTarget, rReference) || IsContained(rReference, rTarget) || IsIntersect(rTarget, rReference))
		return false;
	bool xintersect = !((rTarget.m_min[0] > rReference.m_max[0])||(rTarget.m_max[0] < rReference.m_min[0]));
	bool yintersect = !((rTarget.m_min[1] > rReference.m_max[1])||(rTarget.m_max[1] < rReference.m_min[1]));
	//anti-clockwise
	if(rTarget.m_max[0] < rReference.m_min[0] && rTarget.m_max[1] < rReference.m_min[1]) // 1 -- SW
	{
		p[0].x = (rTarget.m_min[0] + rReference.m_min[0]) / 2;
		p[0].y = (rTarget.m_min[1] + rReference.m_min[1]) / 2;
	}
	else if(xintersect && rTarget.m_max[1] < rReference.m_min[1]) // 2 -- S
	{
		p[0].x = (rTarget.m_min[0] + rReference.m_min[0]) / 2;                   //0 - min
		p[0].y = (rTarget.m_min[1] + rReference.m_min[1]) / 2;
		p[1].x = (rTarget.m_max[0] + rReference.m_max[0]) / 2;					//1 - max
		p[1].y = p[0].y;
	}
	else if(rTarget.m_min[0] > rReference.m_max[0] && rTarget.m_max[1] < rReference.m_min[1]) // 3 -- SE 
	{
		p[0].x = (rTarget.m_max[0] + rReference.m_max[0]) / 2;
		p[0].y = (rTarget.m_min[1] + rReference.m_min[1]) / 2;
	}
	else if(yintersect && rTarget.m_min[0] > rReference.m_max[0]) // 4 -- E
	{
		p[0].x = (rTarget.m_max[0] + rReference.m_max[0]) / 2;					
		p[0].y = (rTarget.m_max[1] + rReference.m_max[1]) / 2;					// 0 - max
		p[1].x = p[0].x;
		p[1].y = (rTarget.m_min[1] + rReference.m_min[1]) / 2;					// 1 - min
	}
	else if(rTarget.m_min[0] > rReference.m_max[0] && rTarget.m_min[1] > rReference.m_max[1]) // 5 -- NE
	{
		p[0].x = (rTarget.m_max[0] + rReference.m_max[0]) / 2;
		p[0].y = (rTarget.m_max[1] + rReference.m_max[1]) / 2;
	}
	else if(xintersect && rTarget.m_min[1] > rReference.m_max[1]) // 6 -- N
	{
		p[0].x = (rTarget.m_min[0] + rReference.m_min[0]) / 2;					// 0 - min
		p[0].y = (rTarget.m_max[1] + rReference.m_max[1]) / 2;
		p[1].x = (rTarget.m_max[0] + rReference.m_max[0]) / 2;					// 1 - max
		p[1].y = p[0].y;
	}
	else if(rTarget.m_max[0] < rReference.m_min[0] && rTarget.m_min[1] > rReference.m_max[1]) // 7 -- NW
	{
		p[0].x = (rTarget.m_min[0] + rReference.m_min[0]) / 2;
		p[0].y = (rTarget.m_max[1] + rReference.m_max[1]) / 2;
	}
	else if(yintersect && rTarget.m_max[0] < rReference.m_min[0]) // 8 -- W
	{
		p[0].x = (rTarget.m_min[0] + rReference.m_min[0]) / 2;
		p[0].y = (rTarget.m_max[1] + rReference.m_max[1]) / 2;					// 0 - max
		p[1].x = p[0].x;
		p[1].y = (rTarget.m_min[1] + rReference.m_min[1]) / 2;					// 1 - min
	}
	return true;
}

// notice orientation is rb on ra, namely, rb relatively to reference ra.
void ExternalSolver::GetTypicalRatio(const RTree<unsigned int,double,2,double>::Rect rTarget, 
									 const RTree<unsigned int,double,2,double>::Rect rReference, 
									 const int orientation, 
									 double *ratio)
{
	switch(orientation)
	{
	case 1:
		ratio[0] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		ratio[1] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		break;
	case 2:// S ratio[1], ratio[2]
		ratio[2] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		ratio[3] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_max[1] - rReference.m_min[1]); 
		break;
	case 3:
		ratio[4] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		ratio[5] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		break;
	case 4:// E ratio[6], ratio[7]
		//remind special cases
		//assert(rReference.m_min[1] != rTarget.m_max[1]);
		if(rReference.m_min[1] == rTarget.m_max[1])
			ratio[6] = 0; // to represent the +infinity
		else
			ratio[6] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		//assert(rReference.m_max[1] != rTarget.m_min[1]);
		if(rReference.m_max[1] == rTarget.m_min[1])
			ratio[7] = 0; // to represent the +infinity
		else
			ratio[7] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		break;
	case 5:
		ratio[8] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		ratio[9] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		break;
	case 6:// N ratio[10], ratip[11]
		ratio[10] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		ratio[11] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		break;
	case 7:
		ratio[12] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		ratio[13] = (rReference.m_max[0] - rTarget.m_min[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		break;
	case 8:// W ratio[14], ratio[15]
		//remind special cases
		//assert(rReference.m_min[1] != rTarget.m_max[1]);
		if(rReference.m_min[1] == rTarget.m_max[1])
			ratio[14] = 0; // to represent the  +infinity
		else
			ratio[14] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_max[1] - rReference.m_min[1]);
		//assert(rReference.m_max[1] != rTarget.m_min[1]);
		if(rReference.m_max[1] == rTarget.m_min[1])
			ratio[15] = 0;// to respresent the  +infinity
		else
			ratio[15] = (rReference.m_min[0] - rTarget.m_max[0]) / (rTarget.m_min[1] - rReference.m_max[1]);
		break;
	default:
		break;
	}
}

//void ExternalSolver::GetIntersect(const Point pa,
bool ExternalSolver::GetIntersect(const Point pa, 
								  const double ra, 
								  const Point pb, 
								  const double rb,
								  Point & intersect)
{
	//assert(rb!=ra)
	if(rb==ra)
		return false;;
	intersect.x = (pa.y - pb.y - pa.x*ra + pb.x*rb) / (rb - ra);
	intersect.y = ra*(intersect.x - pa.x)+pa.y;
	return true;
}

bool ExternalSolver::CheckPointAndLine(const Point mid, 
									   const double ratio,
									   const Point p)
{
	return p.y >= (ratio*(p.x - mid.x) + mid.y);
}

void ExternalSolver::SplitRect(vector<RTree<unsigned int,double,DIMENSION,double>::Rect>* sourceR, 
							   RTree<unsigned int,double,DIMENSION,double>::Rect* spliterR, 
							   vector<RTree<unsigned int, double, DIMENSION, double>::Rect>& targetR)
{
	for(vector<RTree<unsigned int, double, DIMENSION, double>::Rect>::iterator itS = sourceR->begin();
		itS != sourceR->end();
		itS++)
	{
		SplitHelper(&(*itS), spliterR, targetR);
	}
}

void ExternalSolver::SplitHelper(const RTree<unsigned int, double, DIMENSION, double>::Rect * source ,
								 const RTree<unsigned int, double, DIMENSION, double>::Rect * spliter,
								 vector<RTree<unsigned int, double, DIMENSION, double>::Rect>& answer)
{
	RTree<unsigned int, double, DIMENSION, double>::Rect overlap;
	overlap.m_min[0] = max(source->m_min[0], spliter->m_min[0]);
	overlap.m_min[1] = max(source->m_min[1], spliter->m_min[1]);
	overlap.m_max[0] = min(source->m_max[0], spliter->m_max[0]);
	overlap.m_max[1] = min(source->m_max[1], spliter->m_max[1]);
	if(overlap.m_min[0] >= overlap.m_max[0] ||
		overlap.m_min[1] >= overlap.m_max[1])
	{
		answer.push_back(*source);
		return;
	}
	//if overlap == source
	if(overlap.m_min[0] == source->m_min[0] &&
		overlap.m_min[1] == source->m_min[1] &&
		overlap.m_max[0] == source->m_max[0] &&
		overlap.m_max[1] == source->m_max[1])
	{
		return;
	}
	// if overlap == spliter
	if(overlap.m_min[0] == spliter->m_min[0] &&
		overlap.m_min[1] == spliter->m_min[1] &&
		overlap.m_max[0] == spliter->m_max[0] &&
		overlap.m_max[1] == spliter->m_max[1])
	{
		RTree<unsigned int, double, DIMENSION, double>::Rect tempRect1, tempRect2, tempRect3, tempRect4;
		// 1 left
		tempRect1.m_min[0] = source->m_min[0];
		tempRect1.m_min[1] = source->m_min[1];
		tempRect1.m_max[0] = spliter->m_min[0];
		tempRect1.m_max[1] = source->m_max[1];
		// make sure it's a rectangle
		if(tempRect1.m_max[0] - tempRect1.m_min[0] > EPSILON &&
			tempRect1.m_max[1] - tempRect1.m_min[1] > EPSILON )
		answer.push_back(tempRect1);
		// 2 lower
		tempRect2.m_min[0] = spliter->m_min[0];
		tempRect2.m_min[1] = source->m_min[1];
		tempRect2.m_max[0] = spliter->m_max[0];
		tempRect2.m_max[1] = spliter->m_min[1];
		if(tempRect2.m_max[0] - tempRect2.m_min[0] > EPSILON &&
			tempRect2.m_max[1] - tempRect2.m_min[1] > EPSILON )
		answer.push_back(tempRect2);
		// 3 upper
		tempRect3.m_min[0] = spliter->m_min[0];
		tempRect3.m_min[1] = spliter->m_max[1];
		tempRect3.m_max[0] = spliter->m_max[0];
		tempRect3.m_max[1] = source->m_max[1];
		if(tempRect3.m_max[0] - tempRect3.m_min[0] > EPSILON &&
			tempRect3.m_max[1] - tempRect3.m_min[1] > EPSILON )
		answer.push_back(tempRect3);
		// 4 right
		tempRect4.m_min[0] = spliter->m_max[0];
		tempRect4.m_min[1] = source->m_min[1];
		tempRect4.m_max[0] = source->m_max[0];
		tempRect4.m_max[1] = source->m_max[1];
		if(tempRect4.m_max[0] - tempRect4.m_min[0] > EPSILON &&
			tempRect4.m_max[1] - tempRect4.m_min[1] > EPSILON )
		answer.push_back(tempRect4);
		return;
	}
	// how many corners of spliter is contained in source
	int countContainedCorner = 0;
	vector<int> containedRecord;
	if( source->m_max[0] - overlap.m_min[0] > EPSILON &&
		overlap.m_min[0] - source->m_min[0] > EPSILON &&
		source->m_max[1] - overlap.m_min[1] > EPSILON &&
		overlap.m_min[1] - source->m_min[1] > EPSILON)	
	{
		countContainedCorner++;
		containedRecord.push_back(0);
	}
	if( source->m_max[0] - overlap.m_max[0] > EPSILON &&
		overlap.m_max[0] - source->m_min[0] > EPSILON &&
		source->m_max[1] - overlap.m_max[1] > EPSILON &&
		overlap.m_max[1] - source->m_min[1] > EPSILON)
	{
		countContainedCorner++;
		containedRecord.push_back(2);
	}
	if( source->m_max[0] - overlap.m_min[0] > EPSILON &&
		overlap.m_min[0] - source->m_min[0] > EPSILON &&
		source->m_max[1] - overlap.m_max[1] > EPSILON &&
		overlap.m_max[1] - source->m_min[1] > EPSILON)
	{
		countContainedCorner++;
		containedRecord.push_back(3);
	}
	if( source->m_max[0] - overlap.m_max[0] > EPSILON &&
		overlap.m_max[0] - source->m_min[0] > EPSILON &&
		source->m_max[1] - overlap.m_min[1] > EPSILON &&
		overlap.m_min[1] - source->m_min[1] > EPSILON)
	{
		countContainedCorner++;
		containedRecord.push_back(1);
	}
	if(countContainedCorner == 1)
	{
		int containedCorner = containedRecord[0];
		RTree<unsigned int, double, DIMENSION, double>::Rect tempRectA, tempRectB;
		switch(containedCorner)
		{
			case 0:
					tempRectA.m_min[0] = source->m_min[0];
					tempRectA.m_min[1] = source->m_min[1];
					tempRectA.m_max[0] = overlap.m_min[0];
					tempRectA.m_max[1] = source->m_max[1];
					tempRectB.m_min[0] = overlap.m_min[0];
					tempRectB.m_min[1] = source->m_min[1];
					tempRectB.m_max[0] = source->m_max[0];
					tempRectB.m_max[1] = overlap.m_min[1];
					break;
				case 1:
					tempRectA.m_min[0] = source->m_min[0];
					tempRectA.m_min[1] = source->m_min[1];
					tempRectA.m_max[0] = overlap.m_max[0];
					tempRectA.m_max[1] = overlap.m_min[1];
					tempRectB.m_min[0] = overlap.m_max[0];
					tempRectB.m_min[1] = source->m_min[1];
					tempRectB.m_max[0] = source->m_max[0];
					tempRectB.m_max[1] = source->m_max[1];
					break;
				case 2:
					tempRectA.m_min[0] = source->m_min[0];
					tempRectA.m_min[1] = overlap.m_max[1];
					tempRectA.m_max[0] = overlap.m_max[0];
					tempRectA.m_max[1] = source->m_max[1];
					tempRectB.m_min[0] = overlap.m_max[0];
					tempRectB.m_min[1] = source->m_min[1];
					tempRectB.m_max[0] = source->m_max[0];
					tempRectB.m_max[1] = source->m_max[1];
					break;
				case 3:
					tempRectA.m_min[0] = source->m_min[0];
					tempRectA.m_min[1] = source->m_min[1];
					tempRectA.m_max[0] = overlap.m_min[0];
					tempRectA.m_max[1] = source->m_max[1];
					tempRectB.m_min[0] = overlap.m_min[0];
					tempRectB.m_min[1] = overlap.m_max[1];
					tempRectB.m_max[0] = source->m_max[0];
					tempRectB.m_max[1] = source->m_max[1];
					break;
				default:
					break;
				}
		answer.push_back(tempRectA);
		answer.push_back(tempRectB);
		return;
	}
	if(countContainedCorner == 2)
	{
		RTree<unsigned int, double, DIMENSION, double>::Rect tempRectA, tempRectB, tempRectC;
		int containedA = containedRecord[0];
		int containedB = containedRecord[1];
		if((containedA == 0 && containedB == 1) ||
			(containedA == 1 && containedB == 0))
		{
			tempRectA.m_min[0] = source->m_min[0];
			tempRectA.m_min[1] = source->m_min[1];
			tempRectA.m_max[0] = overlap.m_min[0];
			tempRectA.m_max[1] = source->m_max[1];
			tempRectB.m_min[0] = overlap.m_min[0];
			tempRectB.m_min[1] = source->m_min[1];
			tempRectB.m_max[0] = overlap.m_max[0];
			tempRectB.m_max[1] = overlap.m_min[1];
			tempRectC.m_min[0] = overlap.m_max[0];
			tempRectC.m_min[1] = source->m_min[1];
			tempRectC.m_max[0] = source->m_max[0];
			tempRectC.m_max[1] = source->m_max[1];
		}
		else if ((containedA == 1 && containedB == 2) ||
				(containedA == 2 && containedB == 1))
		{
			tempRectA.m_min[0] = source->m_min[0];
			tempRectA.m_min[1] = overlap.m_max[1];
			tempRectA.m_max[0] = source->m_max[0];
			tempRectA.m_max[1] = source->m_max[1];
			tempRectB.m_min[0] = overlap.m_max[0];
			tempRectB.m_min[1] = overlap.m_min[1];
			tempRectB.m_max[0] = source->m_max[0];
			tempRectB.m_max[1] = overlap.m_max[1];
			tempRectC.m_min[0] = source->m_min[0];
			tempRectC.m_min[1] = source->m_min[1];
			tempRectC.m_max[0] = source->m_max[0];
			tempRectC.m_max[1] = overlap.m_min[1];
		}
		else if((containedA == 2 && containedB == 3) || 
			(containedA == 3 && containedB == 2))
		{
			tempRectA.m_min[0] = source->m_min[0];
			tempRectA.m_min[1] = source->m_min[1];
			tempRectA.m_max[0] = overlap.m_min[0];
			tempRectA.m_max[1] = source->m_max[1];
			tempRectB.m_min[0] = overlap.m_min[0];
			tempRectB.m_min[1] = overlap.m_max[1];
			tempRectB.m_max[0] = overlap.m_max[0];
			tempRectB.m_max[1] = source->m_max[1];
			tempRectC.m_min[0] = overlap.m_max[0];
			tempRectC.m_min[1] = source->m_min[1];
			tempRectC.m_max[0] = source->m_max[0];
			tempRectC.m_max[1] = source->m_max[1];
		}
		else if((containedA == 0 && containedB == 3) || 
			(containedA == 3 && containedB == 0))
		{
			tempRectA.m_min[0] = source->m_min[0];
			tempRectA.m_min[1] = overlap.m_max[1];
			tempRectA.m_max[0] = source->m_max[0];
			tempRectA.m_max[1] = source->m_max[1];
			tempRectB.m_min[0] = source->m_min[0];
			tempRectB.m_min[1] = overlap.m_min[1];
			tempRectB.m_max[0] = overlap.m_min[0];
			tempRectB.m_max[1] = overlap.m_max[1];
			tempRectC.m_min[0] = source->m_min[0];
			tempRectC.m_min[1] = source->m_min[1];
			tempRectC.m_max[0] = source->m_max[0];
			tempRectC.m_max[1] = overlap.m_min[1];
		}
		answer.push_back(tempRectA);
		answer.push_back(tempRectB);
		answer.push_back(tempRectC);
		return;
	}
	if(countContainedCorner == 0)
	{
		RTree<unsigned int, double, DIMENSION, double>::Rect tempRect;
		if( overlap.m_min[1] - source->m_min[1] > EPSILON)
		{
			tempRect.m_min[0] = source->m_min[0];
			tempRect.m_min[1] = source->m_min[1];
			tempRect.m_max[0] = source->m_max[0];
			tempRect.m_max[1] = overlap.m_min[1];
			answer.push_back(tempRect);
		}
		if( source->m_max[1] - overlap.m_max[1] > EPSILON )
		{
			tempRect.m_min[0] = source->m_min[0];
			tempRect.m_min[1] = overlap.m_max[1];
			tempRect.m_max[0] = source->m_max[0];
			tempRect.m_max[1] = source->m_max[1];
			answer.push_back(tempRect);
		}
		if( overlap.m_min[0] - source->m_min[0] > EPSILON)
		{
			tempRect.m_min[0] = source->m_min[0];
			tempRect.m_min[1] = source->m_min[1];
			tempRect.m_max[0] = overlap.m_min[0];
			tempRect.m_max[1] = source->m_max[1];
			answer.push_back(tempRect);
		}
		if( source->m_max[0] - overlap.m_max[0] > EPSILON)
		{
			tempRect.m_min[0] = overlap.m_max[0];
			tempRect.m_min[1] = source->m_min[1];
			tempRect.m_max[0] = source->m_max[0];
			tempRect.m_max[1] = source->m_max[1];
			answer.push_back(tempRect);
		}
		return;
	}
}

double ExternalSolver::GetMinDistBtwn(const Point p, 
									  const RTree<unsigned int, double, DIMENSION, double>::Rect rect)
{
	if(p.x <= rect.m_max[0]
	&& p.x >= rect.m_min[0]
	&& p.y <= rect.m_max[1]
	&& p.y >= rect.m_min[1])
	{
		return 0.0;
	}
	else if(p.x <= rect.m_max[0] && p.x >= rect.m_min[0])
	{
		if(p.y < rect.m_min[1])
			return (rect.m_min[1] - p.y);
		else// if(p.y > rect.m_max[1])
			return (p.y - rect.m_max[1]);
	}
	else if(p.y <= rect.m_max[1] && p.y >= rect.m_min[1])
	{
		if(p.x < rect.m_min[0])
			return (rect.m_min[0] - p.x);
		else// if(p.x > rect.m_max[0])
			return (p.x - rect.m_max[0]);
	}
	else
	{
		if(p.x < rect.m_min[0] && p.y < rect.m_min[1])
			return GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]);
		if(p.x > rect.m_max[0] && p.y < rect.m_min[1])
			return GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1]);
		if(p.x > rect.m_max[0] && p.y > rect.m_max[1])
			return GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1]);
		if(p.x < rect.m_min[0] && p.y > rect.m_max[1])
			return GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1]);
		else
			return -1.0;
	}
}

double ExternalSolver::GetMinMaxDistBtwn(const Point p, 
									  const RTree<unsigned int,double,DIMENSION,double>::Rect rect)
{
	if(p.x <= rect.m_max[0]
	&& p.x >= rect.m_min[0]
	&& p.y <= rect.m_max[1]
	&& p.y >= rect.m_min[1])
	{
		return 0.0;
	}
	else if(p.x <= rect.m_min[0] && p.y <= rect.m_min[1])
		return min(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1]));
	else if(p.x > rect.m_min[0] && p.x < rect.m_max[0] && p.y < rect.m_min[1])
		return min(max(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1])),
		min(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1])));
	else if(p.x >= rect.m_max[0] && p.y <= rect.m_min[1])
		return min(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1]));
	else if(p.x > rect.m_max[0] && p.y > rect.m_min[1] && p.y < rect.m_max[1])
		return min(max(GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1])),
		min(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1])));
	else if(p.x >= rect.m_max[0] && p.y >= rect.m_max[1])
		return min(GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1]));
	else if(p.x < rect.m_max[0] && p.x > rect.m_min[0] && p.y > rect.m_max[1])
		return min(max(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1])),
		min(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1])));
	else if(p.x <= rect.m_min[0] && p.y >= rect.m_max[1])
		return min(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1]));
	else if(p.x < rect.m_min[0] && p.y > rect.m_min[1] && p.y < rect.m_max[1])
		return min(max(GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_min[1]), GetDistBtwn(p.x, p.y, rect.m_min[0], rect.m_max[1])),
		min(GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_max[1]), GetDistBtwn(p.x, p.y, rect.m_max[0], rect.m_min[1])));
	else
		//error
		return -1.0;
}

struct PlanAnswer
{
	int influence;
	set<int> setC;
};

class PlanAnswerCompare
{
public:
	bool operator()(const PlanAnswer pa, const PlanAnswer pb)
	{
		return pa.influence < pb.influence;
	}
};

class PlanAnswerCompareB
{
public:
	bool operator()(const PlanAnswer pa, const PlanAnswer pb)
	{
		return pa.influence > pb.influence;
	}
};

class PlanSolver:public Solver
{
public:
	PlanSolver(string fileC, string fileF, string fileM)
	{
		ReadFile(fileC, dataC);
		ReadFile(fileF, dataF);
		ReadFile(fileM, dataM);
		answer.reserve(K);
	};
	virtual ~PlanSolver(){};
	void Print();
	void BriefPrint();

protected:
	vector<PlanAnswer> answer;
	vector<Point> dataC;
	vector<Point> dataF;
	vector<Point> dataM;
	long during;
	long start;
	long end;

	void ReadFile(string, vector<Point>&);
	void StartTiming(){start = clock();};
	void EndTiming(){end = clock(); during = end - start;};
};

void PlanSolver::ReadFile(string filename, vector<Point>& data)
{
	const char* fn;
	fstream filestr;
	fn = filename.c_str();
	char dot;
	int id;
	Point eachPoint;
	filestr.open(fn, fstream::in);
	while(filestr>>eachPoint.x>>dot>>eachPoint.y>>dot>>id)
	{
		data.push_back(eachPoint);
	}
	filestr.close();
}

void PlanSolver::Print()
{
	PlanAnswerCompareB comp;
	sort(answer.begin(), answer.end(), comp);
	for(int i = 0; i < K; i++)
	{
		cout<<i+1<<". "<<answer[i].influence<<"  ";
		for(set<int>::iterator it = answer[i].setC.begin(); it != answer[i].setC.end(); it++)
		{
			cout<<"[ "<<dataC[*it].x<<" , "<<dataC[*it].y<<" ]  ";
		}
		cout<<endl;
	}
	cout<<"C"<<FILE_C_SIZE<<" F"<<FILE_F_SIZE<<" M"<<FILE_M_SIZE<<endl;
	cout<<"N="<<PARAN<<" K="<<K<<endl;
	cout<<"NLCPlan costs "<<during<<" ms."<<endl;
}

void PlanSolver::BriefPrint()
{
	PlanAnswerCompareB comp;
	sort(answer.begin(), answer.end(), comp);
	int answerSize = answer.size();
	/*
	for(int i = 0; i < answerSize; i++)
	{
		cout<<i+1<<". "<<answer[i].influence<<"  ";
		for(set<int>::iterator it = answer[i].setC.begin(); it != answer[i].setC.end(); it++)
		{
			cout<<*it<<'\t';
		}
		cout<<endl;
	}
	*/
	cout<<"C"<<FILE_C_SIZE<<" F"<<FILE_F_SIZE<<" M"<<FILE_M_SIZE<<endl;
	cout<<"N="<<PARAN<<" K="<<K<<endl;
	cout<<"NLCPlan costs "<<during<<" ms."<<endl;
}

#endif SOLVER_H