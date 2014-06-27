#include "Solver.h"
///------------------------------------------------------
// 1. check out the unstable feature
// 2. check out the min max strategy
// 3. check out the pruning strategy
///------------------------------------------------------


///------------------------------------------
/// Bounding Influece Pruning Solver
///------------------------------------------

#define ALPHA 2

	//static RTree<unsigned int, double, DIMENSION, double>::Rect relevantM [MA];
	static RTree<unsigned int, double, DIMENSION, double>::Rect eachRelevantM [FILE_M_SIZE];
	static unsigned int tagRelevantM[FILE_M_SIZE];
	static RTree<unsigned int, double, DIMENSION, double>::Rect eachRelevantF [FILE_F_SIZE];
	static unsigned int tagRelevantF[FILE_F_SIZE];
struct Slot
{
	RTree<unsigned int, double, DIMENSION, double>::Branch *branch;
	vector<RTree<unsigned int, double, DIMENSION, double>::Branch *> innerRelevant;
	vector<RTree<unsigned int, double, DIMENSION, double>::Branch *> outerRelevant;
	int infUpperBound;
	RTree<unsigned int, double, DIMENSION, double>::Rect rectM;
	vector<RTree<unsigned int, double, DIMENSION, double>::Rect> finalRectM;
	double key;
};

class SlotCompare
{
public:
	bool operator()(const Slot sa, const Slot sb)
	{
		return sa.key < sb.key;
	}
};

class BranchPCompare
{
public:
	BranchPCompare(RTree<unsigned int, double, DIMENSION, double>::Rect *r)
	{
		targetRect = r;
	};
	bool operator()(const RTree<unsigned int, double, DIMENSION, double>::Branch * ba,
		const RTree<unsigned int, double, DIMENSION, double>::Branch * bb)
	{
		return GetMinMaxDistBtwn(ba->m_rect, *targetRect) < GetMinMaxDistBtwn(bb->m_rect, *targetRect);
	};
protected:
	RTree<unsigned int, double, DIMENSION, double>::Rect *targetRect;

	double GetMinMaxDistBtwn(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	bool IsContained(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	bool IsIntersect(const RTree<unsigned int, double, DIMENSION, double>::Rect,
		const RTree<unsigned int, double, DIMENSION, double>::Rect);
	double GetDistBtwnS(const double ax, const double ay, const double bx, const double by)
	{return (ax-bx)*(ax-bx) + (ay-by)*(ay-by);};
	double GetDistBtwn(const double ax, const double ay, const double bx, const double by)
	{return sqrt(GetDistBtwnS(ax, ay, bx, by));};
};

// returns -1 if rb is contained in ra, 0 if ra and rb intersect, and minMax(ra, rb) dist otherwise
double BranchPCompare::GetMinMaxDistBtwn(const RTree<unsigned int, double, DIMENSION, double>::Rect ra,
					  const RTree<unsigned int, double, DIMENSION, double>::Rect rb)
{
	if(IsContained(rb, ra))
		return -1.0;
	if(IsIntersect(ra, rb))
		return 0.0;
	bool xintersect = !((rb.m_min[0] > ra.m_max[0])||(rb.m_max[0] < ra.m_min[0]));
	bool yintersect = !((rb.m_min[1] > ra.m_max[1])||(rb.m_max[1] < ra.m_min[1]));
	//anti-clockwise
	if(rb.m_max[0] < ra.m_min[0] && rb.m_max[1] < ra.m_max[1])
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

bool BranchPCompare::IsContained(const RTree<unsigned int, double, DIMENSION, double>::Rect ra,
								 const RTree<unsigned int, double, DIMENSION, double>::Rect rb)
{
	if(ra.m_min[0] > rb.m_min[0] &&
		ra.m_min[1] > rb.m_min[1] &&
		ra.m_max[0] < rb.m_min[0] &&
		ra.m_max[1] < rb.m_max[1])
		return true;
	else return false;
}

bool BranchPCompare::IsIntersect(const RTree<unsigned int, double, DIMENSION, double>::Rect ra,
								 const RTree<unsigned int, double, DIMENSION, double>::Rect rb)
{
	bool xintersect = !((rb.m_min[0] > ra.m_max[0])||(rb.m_max[0] < ra.m_min[0]));
	bool yintersect = !((rb.m_min[1] > ra.m_max[1])||(rb.m_max[1] < ra.m_min[1]));
	if(xintersect && yintersect)
		return true;
	else return false;
}

class BranchSizeCompare
{
public:
	bool operator()(const RTree<unsigned int, double, DIMENSION, double>::Branch* ba,
		const RTree<unsigned int, double, DIMENSION, double>::Branch* bb)
	{
		double areaA = (ba->m_rect.m_max[0] - ba->m_rect.m_min[0])*(ba->m_rect.m_max[1] - ba->m_rect.m_min[1]);
		double areaB = (bb->m_rect.m_max[0] - bb->m_rect.m_min[0])*(bb->m_rect.m_max[1] - bb->m_rect.m_min[1]);
		return areaA < areaB;
	}
};

/// --------------------------------------------------------------------------------------------------------
/// Class BIPSolver
/// The implementation of Bounding Influence Pruning algorithm
/// Inherits ExternalSolver class
/// Constructor parameters : fileC, fileF, fileM
/// Member functions:
///		public:			void Solve()				Solves the given problem
///		protected:		void InitializeHeap()		Constructs the running heap for C
///						void DoTop()				Does work for the top Slot on the heap
///						void ComputeLeaf()			Computes the exact influence value for c in Slot	
///						void Update()				Computes values in Slot for internal nodes (MBRc)
///						void InitializeFSet()		Generates the initial expanded F sets
///						void ExpandInnerFor()		Expands an intersected (inner) relevant nodes
///														for a given internal node (MBRc)
///						void ExpandInitialFor()		Expands a record on currentFSet for a given
///														initial node (MBRc) on the initial heap
///						void FindIntersect()		Finds intersection points for current internal
///														node (MBRc) based on a given outer relevant node
///						void IsPruned()				Attempts to prune a given outer relevant node
///														based on intersection point found
///						void EstimateUpperBound()	Estimate the upper bounding influence value for an
///														internal node (MBRc) based on intersection points
///														(and maybe domainance points or even mutual pruning)
/// --------------------------------------------------------------------------------------------------------

class BIPSolver: public ExternalSolver
{
public:
	//BIPSolver arguemant table modified by ZY
	BIPSolver(string fileM, string fileF, string fileC, int ma, int fa, int ca):ExternalSolver(fileM, fileF, fileC, ma, fa, ca)
	{
		rootC = treeC.SooneGetRoot();
		rootF = treeF.SooneGetRoot();
		rootM = treeM.SooneGetRoot();
		treeC.Count();
		treeF.Count();
		treeM.Count();
		heap.reserve(paraC);
		infThreshold = 0;
		countC = 0;
		countPrune = 0;
		during  = 0;
		duringLeaf = 0;
		endBuilding = clock();
		countNodeC = treeC.SooneCount(countLeafC);
		countNodeM = treeM.SooneCount(countLeafM);
		countNodeF = treeF.SooneCount(countLeafF);
		countAccessC = 0;
		countAccessF = 0;
		countAccessM = 0;
		countComputeC = 0;
		countAccessDataF = 0;
		countAccessDataM = 0;
		counterIOC = 0;
		counterIOF = 0;
		counterIOM = 0;
		counterComputeC = 0;
		counterComputeF = 0;
		counterComputeM = 0;
		heapSize = 0;
		heapMaxSize = 0;
		heapCounter = 0;
		heapSizeSum = 0;
		slotSizeSum = 0;
		heapCounter = 0;
		slotCounter = 0;
		exactComputationCounter = 0;
		pruningComputationCounter = 0;
	};
	//added by ZY
	void PrintDuration(){cout<<"BIP: "<< GetDuring()<<" ms"<<endl;}
	void PrintComputation();
	//end added
	void Solve();
	int GetEndDuring(){return endBuilding;};
	int GetDuring(){return during;};
	int GetLeafDuring(){return duringLeaf;};
	int CountNodeC(){return countNodeC;};
	int CountNodeM(){return countNodeM;};
	int CountNodeF(){return countNodeF;};
	int CountLeafC(){return countLeafC;};
	int CountLeafF(){return countLeafF;};
	int CountLeafM(){return countLeafM;};
	int CountAccessC(){return countAccessC;};
	//int CountComputeC(){return countComputeC;};
	int CountAccessF(){return countAccessF;};
	int CountAccessM(){return countAccessM;};
	int CountAccessDataF(){return countAccessDataF;};
	int CountAccessDataM(){return countAccessDataM;};
	long long CountIOC(){return counterIOC;};
	long long CountIOF(){return counterIOF;};
	long long CountIOM(){return counterIOM;};
	long long CountComputeC(){return counterComputeC;};
	long long CountComputeF(){return counterComputeF;};
	long long CountComputeM(){return counterComputeM;};
	long long GetHeapAverageSize(){return (int)(heapSizeSum / heapCounter)*(slotSizeSum / slotCounter);};
	long long GetHeapMaxSize(){return heapMaxSize *(slotSizeSum/ slotCounter);};
	long long GetExactComputationCounter(){return exactComputationCounter;};
	long long GetPruningComputationCounter(){return pruningComputationCounter;};
	//int CountNodeCA(){return treeC->a}
protected:
	RTree<unsigned int, double, DIMENSION, double>::Node * rootC;
	RTree<unsigned int, double, DIMENSION, double>::Node * rootF;
	RTree<unsigned int, double, DIMENSION, double>::Node * rootM;
	SlotCompare comp;
	BranchSizeCompare compBS;	
	
	vector<Slot> heap;
	set<RTree<unsigned int, double, DIMENSION, double>::Branch *> currentFSet;
	
	int infThreshold;
	long start;
	long end;
	long during;
	long startLeaf;
	long endLeaf;
	long duringLeaf;
	long endBuilding;
	//int fa [CA];
	//int ma [CA]; 
	int countC;
	int countPrune;
	int countNodeC;
	int countNodeM;
	int countNodeF;
	int countLeafC;
	int countLeafM;
	int countLeafF;
	int countAccessC;
	int countAccessM;
	int countAccessF;
	int countComputeC;
	int countAccessDataF;
	int countAccessDataM;
	long long counterIOC;
	long long counterIOF;
	long long counterIOM;
	long long counterComputeC;
	long long counterComputeF;
	long long counterComputeM;
	long long heapSizeSum;
	long long slotSizeSum;
	int heapSize;
	int heapCounter;
	int slotCounter;
	int heapMaxSize;

	long long exactComputationCounter;
	long long pruningComputationCounter;

	void InitializeHeap();
	void InitializeFSet();
	
	void DoTop();


	bool PruneLeaf(Slot *);
	void EnhancedComputeLeaf(Slot *);
	void BIPSolver::ExpandRectForLeaf(Slot *, RTree<unsigned int, double, DIMENSION, double>::Rect,
		RTree<unsigned int, double, DIMENSION, double>::Branch *);
	bool PruneMBRF(RTree<unsigned int, double, DIMENSION, double>::Branch *, Point *, double);

	void Update(Slot *);
	void ExpandFor(Slot *, RTree<unsigned int, double, DIMENSION, double>::Branch*);
	void ExpandForInner(Slot *, RTree<unsigned int, double, DIMENSION, double>::Branch*);
	//void ExpandForOuter(Slot *, RTree<unsigned int, double, DIMENSION, double>::Branch*);
	bool FindIntersect(
		vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator, 
		Slot *,
		Point *,
		Point *,
		Point *,
		double *,
		bool *);
	bool IsPruned(
		vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator, 
		Slot *,
		Point*,
		bool*,
		double*);
	int EstimateUpperBound(
		Slot *, 
		Point*, 
		Point *, 
		bool*);


	int countF(RTree<unsigned int, double, DIMENSION, double>::Branch *);
};
//added by ZY
void BIPSolver::PrintComputation()
{
	cout << "Total Computation: " << GetExactComputationCounter() +  GetPruningComputationCounter() << endl;
	cout << "Pruning Computation: " << GetPruningComputationCounter() <<endl;
	cout << "Exact Computation: " << GetExactComputationCounter() <<endl;
}
//end added

void BIPSolver::Solve()
{	
	start = clock();
	answer.clear();
	during = 0;
	duringLeaf = 0;
	InitializeHeap();
	while(heap.size() != 0 && heap.front().infUpperBound > infThreshold)
	{
		DoTop();
	}
	end= clock();
	during = end - start;
}

void BIPSolver::InitializeHeap()
{
	bool workDone = false;		
	//count the access
	countAccessC++;
	for(int i = 0; i < rootC->m_count; i++)
	{
		Slot eachSlot;
		eachSlot.branch = &(rootC->m_branch[i]);
		eachSlot.key = GetArea(eachSlot.branch->m_rect);
		heap.push_back(eachSlot);

	}
	make_heap(heap.begin(), heap.end(), comp);
	while(!workDone)
	{
		Slot curSlot = heap.front();
		pop_heap(heap.begin(), heap.end(), comp);
		heap.pop_back();

		if(GetAreaRoot(curSlot.branch->m_rect) > (RANGE / ALPHA))
		{			
			//count the access
			countAccessC++;
			for(int i = 0; i < curSlot.branch->m_child->m_count; i++)
			{
				Slot eachSlot;
				eachSlot.branch = &(curSlot.branch->m_child->m_branch[i]);
				eachSlot.key = GetArea(eachSlot.branch->m_rect);
				heap.push_back(eachSlot);
				push_heap(heap.begin(), heap.end(), comp);
			}
		}
		else
		{
			workDone = true;
			heap.push_back(curSlot);
			push_heap(heap.begin(), heap.end(), comp);
		}
	}
	InitializeFSet();
	int testF = 0;
	set<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itCF;

	// initialize MBRs of c
	for(vector<Slot>::iterator itV = heap.begin();
		itV != heap.end();
		itV++)
	{
		Update(&(*itV));
	}
	pruningComputationCounter += heap.size();
	make_heap(heap.begin(), heap.end(), comp);
	heapSizeSum = heap.size();
	heapMaxSize = heap.size();
	heapCounter++;
}

int BIPSolver::countF(RTree<unsigned int, double, DIMENSION, double>::Branch * branch)
{
	if(branch->m_child->IsLeaf())
		return branch->m_child->m_count;
	else
	{
		int count = 0;
		for(int i = 0; i < branch->m_child->m_count; i++)
		{
			count+= countF(&(branch->m_child->m_branch[i]));
		}
		return count;
	}
}

void BIPSolver::InitializeFSet()
{
	vector<RTree<unsigned int, double, DIMENSION, double>::Branch *> heapF;
	//count the access
	countAccessF++;
	for(int i = 0; i < rootF->m_count; i++)
	{
		heapF.push_back(&(rootF->m_branch[i]));
	}
	make_heap(heapF.begin(), heapF.end(), compBS);
	bool workDone = false;
	while(!workDone)
	{
		RTree<unsigned int, double, DIMENSION, double>::Branch * topBranch = heapF.front();
		pop_heap(heapF.begin(), heapF.end(), compBS);
		heapF.pop_back();

		if(GetAreaRoot(topBranch->m_rect) > (RANGE) / ALPHA)
		{				
			//count the access
			countAccessF++;
			for(int i =0; i < topBranch->m_child->m_count; i++)
			{
				heapF.push_back(&(topBranch->m_child->m_branch[i]));
				push_heap(heapF.begin(), heapF.end(), compBS);
			}
		}
		else
		{
			workDone = true;
			heapF.push_back(topBranch);
			push_heap(heapF.begin(), heapF.end(), compBS);
		}
	}
	for(vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator it = heapF.begin();
		it != heapF.end();
		it++)
	{
		currentFSet.insert(*it);
	}
}

void BIPSolver::DoTop()
{
	counterIOC++;
	int tempSize = heap.size();
	heapSizeSum += tempSize;
	heapMaxSize = tempSize > heapMaxSize ? tempSize : heapMaxSize;
	heapCounter++;
	slotCounter++;
	Slot curC = heap.front();
	slotSizeSum += sizeof(curC.branch) + sizeof(curC.key) + 
		sizeof(curC.innerRelevant.size() * sizeof(curC.branch)) +
		sizeof(curC.outerRelevant.size() * sizeof(curC.branch));
	pop_heap(heap.begin(), heap.end(), comp);
	heap.pop_back();
	// check to see if it is a leaf
	if(curC.branch->m_rect.m_min[0] == curC.branch->m_rect.m_max[0]
	&& curC.branch->m_rect.m_min[1] == curC.branch->m_rect.m_max[1])
	{
		EnhancedComputeLeaf(&curC);
		//count the access
		//countComputeC++;
		//countAccessC++;
	}
	else
	{	
		//count the access
		countAccessC++;
		for(int i = 0; i < curC.branch->m_child->m_count; i++)
		{
			Slot eachSlot;
			eachSlot.outerRelevant = curC.outerRelevant;
			eachSlot.innerRelevant = curC.innerRelevant;
			eachSlot.branch = &(curC.branch->m_child->m_branch[i]);
			Update(&eachSlot);
			if(eachSlot.branch->m_rect.m_min[0] == eachSlot.branch->m_rect.m_max[0] &&
				eachSlot.branch->m_rect.m_min[1] == eachSlot.branch->m_rect.m_max[1])
			{
				EnhancedComputeLeaf(&eachSlot);
				continue;
			}
			heap.push_back(eachSlot);
			pruningComputationCounter++;
			push_heap(heap.begin(), heap.end(), comp);
		}
	}
}

bool BIPSolver::PruneLeaf(Slot * slotC)
{
	bool foundInFour = false;
	bool foundOne = false;
	bool foundTwo = false;
	bool foundThree = false;
	bool foundFour = false;
	Point midOne, midTwo, midThree, midFour;
	double ratioOne, ratioTwo, ratioThree, ratioFour;
	int slotCSize = slotC->outerRelevant.size();
	counterIOF += slotCSize;//the first access to all outer relevant F is counted
	for(int i = 0; i < slotCSize; i++)
	{
		if(!foundOne && slotC->outerRelevant[i]->m_rect.m_min[0] >= slotC->branch->m_rect.m_min[0]
		&& slotC->outerRelevant[i]->m_rect.m_min[1] >= slotC->branch->m_rect.m_min[1])
		{
			foundOne = true;
			midOne.x = (slotC->outerRelevant[i]->m_rect.m_max[0] + slotC->branch->m_rect.m_min[0]) / 2;
			midOne.y = (slotC->outerRelevant[i]->m_rect.m_max[1] + slotC->branch->m_rect.m_min[1]) / 2;
			ratioOne = (midOne.x - slotC->branch->m_rect.m_min[0])/(slotC->branch->m_rect.m_min[1] - midOne.y);
		}
		else if(!foundTwo && slotC->outerRelevant[i]->m_rect.m_max[0] <= slotC->branch->m_rect.m_min[0]
		&& slotC->outerRelevant[i]->m_rect.m_min[1] >= slotC->branch->m_rect.m_min[1])
		{
			foundTwo = true;
			midTwo.x = (slotC->outerRelevant[i]->m_rect.m_min[0] + slotC->branch->m_rect.m_min[0]) / 2;
			midTwo.y = (slotC->outerRelevant[i]->m_rect.m_max[1] + slotC->branch->m_rect.m_min[1]) / 2;
			ratioTwo = (midTwo.x - slotC->branch->m_rect.m_min[0])/(slotC->branch->m_rect.m_min[1] - midTwo.y);
		}
		else if(!foundThree && slotC->outerRelevant[i]->m_rect.m_max[0] <= slotC->branch->m_rect.m_min[0]
		&& slotC->outerRelevant[i]->m_rect.m_max[1] <= slotC->branch->m_rect.m_min[1])
		{
			foundThree = true;
			midThree.x = (slotC->outerRelevant[i]->m_rect.m_min[0] + slotC->branch->m_rect.m_min[0]) / 2;
			midThree.y = (slotC->outerRelevant[i]->m_rect.m_min[1] + slotC->branch->m_rect.m_min[1]) / 2;
			ratioThree = (midThree.x - slotC->branch->m_rect.m_min[0])/(slotC->branch->m_rect.m_min[1] - midThree.y);
		}
		else if(!foundFour && slotC->outerRelevant[i]->m_rect.m_min[0] >= slotC->branch->m_rect.m_min[0]
		&& slotC->outerRelevant[i]->m_rect.m_max[1] <= slotC->branch->m_rect.m_min[1])
		{
			foundFour = true;
			midFour.x = (slotC->outerRelevant[i]->m_rect.m_max[0] + slotC->branch->m_rect.m_min[0]) / 2;
			midFour.y = (slotC->outerRelevant[i]->m_rect.m_min[1] + slotC->branch->m_rect.m_min[1]) / 2;
			ratioFour = (midFour.x - slotC->branch->m_rect.m_min[0])/(slotC->branch->m_rect.m_min[1] - midFour.y);
		}
		foundInFour = foundOne && foundTwo && foundThree && foundFour;
		if(foundInFour)
			break;
	}
	if(!foundInFour)
		return false;
	else
	{
		Point cur;
		cur.x = slotC->branch->m_rect.m_min[0];
		cur.y = slotC->branch->m_rect.m_min[1];
		Point vOne, vTwo, vThree, vFour;
		double tempA = midTwo.y - midOne.y + ratioOne*midOne.x - ratioTwo*midTwo.x;
		double tempB = ratioOne - ratioTwo;
		vOne.x = tempA / tempB;
		vOne.y = ratioOne*(vOne.x - midOne.x) + midOne.y;
		tempA = midThree.y - midTwo.y + ratioTwo*midTwo.x - ratioThree*midThree.x;
		tempB = ratioTwo - ratioThree;
		vTwo.x = tempA / tempB;
		vTwo.y = ratioTwo*(vTwo.x - midTwo.x) + midTwo.y;
		tempA = midFour.y - midThree.y + ratioThree*midThree.x - ratioFour*midFour.x;
		tempB = ratioThree - ratioFour;
		vThree.x = tempA / tempB;
		vThree.y = ratioThree*(vThree.x - midThree.x) + midThree.y;
		tempA = midOne.y - midFour.y + ratioFour*midFour.x - ratioOne*midOne.x;
		tempB = ratioFour - ratioOne;
		vFour.x = tempA / tempB;
		vFour.y = ratioFour*(vFour.x - midFour.x) + midFour.y;
		pruningComputationCounter += 4;
		double rOne = GetDistBtwn(cur, vOne) / 2;
		double rTwo = GetDistBtwn(cur, vTwo) / 2;
		double rThree = GetDistBtwn(cur, vThree) / 2;
		double rFour = GetDistBtwn(cur, vFour) / 2;
		RTree<unsigned int, double, DIMENSION, double>::Rect appRect, refRect;
		appRect.m_min[0] = (vTwo.x - rTwo) < (vThree.x - rThree) ? (vTwo.x - rTwo) : (vThree.x - rThree);
		appRect.m_min[1] = (vThree.y - rThree) < (vFour.y - rFour) ? (vThree.y - rThree) : (vFour.y - rFour);
		appRect.m_max[0] = (vOne.x + rOne) > (vFour.x + rFour) ? (vOne.x + rOne) : (vFour.x + rFour);
		appRect.m_max[1] = (vOne.y + rOne) > (vTwo.y + rTwo) ? (vOne.y + rOne) : (vTwo.y + rTwo);
		appRect.m_min[0] = appRect.m_min[0] < (-1)*RANGE ? (-1)*RANGE : appRect.m_min[0];
		appRect.m_min[1] = appRect.m_min[1] < (-1)*RANGE ? (-1)*RANGE : appRect.m_min[1];
		appRect.m_max[0] = appRect.m_max[0] > RANGE ? RANGE : appRect.m_max[0];
		appRect.m_max[1] = appRect.m_max[1] > RANGE ? RANGE : appRect.m_max[1];
		refRect.m_min[0] = (2*appRect.m_min[0] - cur.x) < (-1)*RANGE ? (-1)*RANGE : (2*appRect.m_min[0] - cur.x);
		refRect.m_min[1] = (2*appRect.m_min[1] - cur.y) < (-1)*RANGE ? (-1)*RANGE : (2*appRect.m_min[1] - cur.y);
		refRect.m_max[0] = (2*appRect.m_max[0] - cur.x) > RANGE ? RANGE : (2*appRect.m_max[0] - cur.x);
		refRect.m_max[1] = (2*appRect.m_max[1] - cur.x) > RANGE ? RANGE : (2*appRect.m_max[1] - cur.y);
		//slotC->finalRectM.clear();
		//slotC->finalRectM.push_back(appRect);
		vector<RTree<unsigned int, double, DIMENSION, double>::Branch*> temp = slotC->outerRelevant;
		slotC->outerRelevant.clear();
		for(int i = 0; i < slotCSize; i++)
		{
			ExpandRectForLeaf(slotC, refRect, temp[i]);
		}
		return true;
	}
}

void BIPSolver::ExpandRectForLeaf(Slot * slotC, RTree<unsigned int, double, DIMENSION, double>::Rect refRect,
								  RTree<unsigned int, double, DIMENSION, double>::Branch * branch)
{
	if(branch->m_rect.m_min[0] >= refRect.m_min[0]
	&& branch->m_rect.m_min[1] >= refRect.m_min[1]
	&& branch->m_rect.m_max[0] <= refRect.m_max[0]
	&& branch->m_rect.m_max[1] <= refRect.m_max[1])
	{
		slotC->outerRelevant.push_back(branch);
	}
	else 
	{
		Point overlapMin, overlapMax;
		overlapMin.x = max(branch->m_rect.m_min[0], refRect.m_min[0]);
		overlapMin.y = max(branch->m_rect.m_min[1], refRect.m_min[1]);
		overlapMax.x = min(branch->m_rect.m_max[0], refRect.m_max[0]);
		overlapMax.y = min(branch->m_rect.m_max[1], refRect.m_max[1]);
		if(overlapMin.x <= overlapMax.x
			&& overlapMin.y <= overlapMax.y)
		{
			if(overlapMin.x == overlapMax.x && overlapMin.y == overlapMax.y)
			{
				slotC->outerRelevant.push_back(branch);
			}
			else
			{
				//counterIOF++;
				int size = branch->m_child->m_count;					
				//count the access
				countAccessF++;
				if(branch->m_child->IsLeaf())// leaf then only one access to all records
					counterIOF++;
				else //internal than multiple accesses
					counterIOF += size;
				for(int i = 0; i < size; i++)
				{
					ExpandRectForLeaf(slotC, refRect, &branch->m_child->m_branch[i]);
				}
			}

		}
	}
}

void BIPSolver::EnhancedComputeLeaf(Slot * slotC)
{
	
	//assert(slotC->innerRelevant.size() == 0);
	PruneLeaf(slotC);
	startLeaf = clock();
	int amountRelevantM;
	vector<RTree<unsigned int, double, DIMENSION, double>::Rect> relevantM;
	int rectSize = slotC->finalRectM.size();
	int countMRecord = 0;
	for(int i = 0; i < rectSize; i++)
	{
		amountRelevantM = treeM.SooneSearchOverlap(slotC->finalRectM[i].m_min, slotC->finalRectM[i].m_max, eachRelevantM, tagRelevantM, countMRecord);
		counterIOM += countMRecord;
		for(int j = 0; j < amountRelevantM; j++)
		{
			relevantM.push_back(eachRelevantM[j]);
		}
	}
	int fSize = slotC->outerRelevant.size();
	Point pointM;
	int mSize = relevantM.size();
	if(mSize < infThreshold)
	{
		Point currentC;
		currentC.x = slotC->branch->m_rect.m_min[0];
		currentC.y = slotC->branch->m_rect.m_min[1];
		answer.insert(pair<int, Point>(0, currentC));
		return;
	}
	int influenceValue = 0;
	countAccessDataM += mSize;
	for(int m = 0; m < mSize; m++)
	{
		counterComputeM++;
		influenceValue++;
		pointM.x = relevantM[m].m_min[0];
		pointM.y = relevantM[m].m_min[1];
		exactComputationCounter++;
		double dist = GetDistBtwn(slotC->branch->m_rect.m_min[0], slotC->branch->m_rect.m_min[1],
			pointM.x, pointM.y);
		for(int i = 0; i < fSize; i++)
		//for(int f = 0; f < FA; f++)
		{
			if(PruneMBRF(slotC->outerRelevant[i], &pointM, dist))
			{
				influenceValue--;
				break;
			}
		}
		if( (mSize - m + influenceValue) < infThreshold)
			break;
	}
	endLeaf = clock();
	Point currentC;
	currentC.x = slotC->branch->m_rect.m_min[0];
	currentC.y = slotC->branch->m_rect.m_min[1];
	counterComputeC++;
	answer.insert(pair<int, Point>(influenceValue, currentC));
	if(answer.size() < K)
		infThreshold = 0;
	else
	{

		multimap<int, Point>::iterator itA = answer.end();
		for(int i = 0; i < K; i++)
			itA--;
		infThreshold = itA->first;
	}
	duringLeaf += endLeaf - startLeaf;
}

bool BIPSolver::PruneMBRF(RTree<unsigned int, double, DIMENSION, double>::Branch * branch,
						  Point *pointM, double dist)
{
	if(branch->m_rect.m_min[0] == branch->m_rect.m_max[0]
	&& branch->m_rect.m_min[1] == branch->m_rect.m_max[1])
	{
		countAccessDataF++;
		counterComputeF++;
		exactComputationCounter++;
		return (GetDistBtwn(pointM->x, pointM->y, branch->m_rect.m_min[0], branch->m_rect.m_min[1]) < dist);
	}
	pruningComputationCounter++;
	if(GetMinDistBtwn(*pointM, branch->m_rect) > dist)
		return false;
	else
	{	
		//count the access
		countAccessF++;
		if(branch->m_child->IsLeaf())// if it is a leaf, than one access
			counterIOF++;
		else // if it is a internal, than multiple accesses
			counterIOF += branch->m_child->m_count;
		for(int i = 0; i < branch->m_child->m_count; i++)
		{

			if(PruneMBRF(&(branch->m_child->m_branch[i]), pointM, dist))
				return true;
		}
		return false;
	}
}

void BIPSolver::Update(Slot * it)
{
	if(it->outerRelevant.empty())
	{
		set<RTree<unsigned int, double, DIMENSION, double>::Branch *> iterationSet = currentFSet;
		for(set<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itF = iterationSet.begin();
			itF != iterationSet.end();
			itF++)
		{
			if(IsIntersect((*itF)->m_rect, it->branch->m_rect)
				|| IsContained(it->branch->m_rect, (*itF)->m_rect))
				ExpandFor(it, *itF);
			else if(IsContained((*itF)->m_rect, it->branch->m_rect))
			{
				it->innerRelevant.push_back(*itF);
			}
			else
			{
				it->outerRelevant.push_back(*itF);
			}
		}
	}

	else
	{
		vector<RTree<unsigned int, double, DIMENSION, double>::Branch *> iterationVector = it->innerRelevant;
		it->innerRelevant.clear();
		for(vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itIR = iterationVector.begin();
			itIR != iterationVector.end();
			itIR++)
		{
			ExpandForInner(it, *itIR);
		}
	}

	// for outers, first use them to generate intersect, then check which of them can be pruned
	BranchPCompare * compBP = new BranchPCompare(&(it->branch->m_rect));
	pruningComputationCounter += it->outerRelevant.size();
	sort(it->outerRelevant.begin(), it->outerRelevant.end(), *compBP);
	// 12 intersects are stored in anti-clockwise direction from the S(3), E(3), N(3), W(3)
	Point foundIntersect[12];
	// 8 regions are stored in anti-clockwise direction, from SW, S, SE, E, NE, N, NW, W
	bool foundInRegion[8] = {false};
	// 4 nearest points in SW, SE, NE, and NW
	Point nnPoint[4];
	// 16 typical ratios
	double foundRatio[16];
	// 12 middle points
	Point midPoint[12];
	//vector<RTree<unsigned int, double, DIMENSION, double>::Branch *> notValid;
	// find intersects first
	vector<bool> validTag (it->outerRelevant.size(), true);
	int count = 0;
	for(vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itOuter = it->outerRelevant.begin();
		itOuter != it->outerRelevant.end();
		itOuter++)
	{
		bool foundAll = true;
		for(int i = 0; i < 8; i++) foundAll = foundAll&&foundInRegion[i];
		if(!foundAll)
		{
			counterIOF++;// find the nearest relevant F
			bool isValid = FindIntersect(itOuter, it, foundIntersect, nnPoint, midPoint, foundRatio, foundInRegion);
			if(!isValid)
				validTag[count] = false;
		}
		else
			validTag[count] = false;
		count++;
	}

	vector<RTree<unsigned int, double, DIMENSION, double>::Branch *> iterationVector = it->outerRelevant;
	it->outerRelevant.clear();
	count = 0;
	for(vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itOR = iterationVector.begin();
		itOR != iterationVector.end();
		itOR++)
	{
		if(!validTag[count])
		{
			counterIOF++; // prune irrelevant F
			bool pruned = IsPruned(itOR, it, foundIntersect,foundInRegion, foundRatio);
			//bool pruned = false;
			if(pruned)
			{
				//if((*itOR)->m_rect.m_min[0] == (*itOR)->m_rect.m_max[0])
					//countPrune++;
				//else countPrune += countF(*itOR);
			}
			else
				it->outerRelevant.push_back(*itOR);
		}
		else
			it->outerRelevant.push_back(*itOR);
		count++;
	}
	// use found intersects to estimate the upper bound influence for it
	// notice this method also set the rectM in it
	int infUpperBound = EstimateUpperBound(it, foundIntersect, nnPoint, foundInRegion);
	it->infUpperBound = infUpperBound;
	it->key = 1 * infUpperBound;
}


void BIPSolver::ExpandFor(Slot *it, RTree<unsigned int, double, DIMENSION, double>::Branch *branch)
{
	// this should not happen, since a point cannot contain a rect or intersect with a rect
	if(branch->m_rect.m_min[0]==branch->m_rect.m_max[0]&&branch->m_rect.m_min[1]==branch->m_rect.m_max[1])
	{
		it->innerRelevant.push_back(branch);
		return;
	}

	//expand branch
	//currentFSet.erase(currentFSet.find(branch));		
	//count the access
	countAccessF++;
	if(branch->m_child->IsLeaf())// if it is a leaf node, than count only one access
		counterIOF++;
	else	// if it is an internal node, than count multiple accesses
		counterIOF += branch->m_child->m_count;
	for(int i = 0; i < branch->m_child->m_count; i++)
	{

		RTree<unsigned int, double, DIMENSION, double>::Branch* childBranch = &branch->m_child->m_branch[i];
		//currentFSet.insert(childBranch);
		if(IsIntersect(childBranch->m_rect, it->branch->m_rect))
			ExpandFor(it, childBranch);
		else if(IsContained(childBranch->m_rect, it->branch->m_rect))
			it->innerRelevant.push_back(childBranch);
		else if(IsContained(it->branch->m_rect, childBranch->m_rect))
			ExpandFor(it, childBranch);
		else
			it->outerRelevant.push_back(childBranch);
	}
}


void BIPSolver::ExpandForInner(Slot *it, RTree<unsigned int, double, DIMENSION, double>::Branch *branch)
{ 
	if(IsContained(branch->m_rect, it->branch->m_rect))
	{
		it->innerRelevant.push_back(branch);
	}
	else if(IsContained(it->branch->m_rect, branch->m_rect)
		|| IsIntersect(it->branch->m_rect, branch->m_rect))
	{
		if(branch->m_rect.m_min[0] == branch->m_rect.m_max[0]
		&& branch->m_rect.m_min[1] == branch->m_rect.m_max[1])
		{
			it->innerRelevant.push_back(branch);
			return;
		}
		//count the access
		countAccessF++;
		if(branch->m_child->IsLeaf()) // if it is a leaf, than only count one access for its children
			counterIOF++;
		else // if it is an internal node, than count multiple accesses for all its children
			counterIOF += branch->m_child->m_count;
		for(int i = 0; i < branch->m_child->m_count; i++)
		{
			//currentFSet.insert(&(branch->m_child->m_branch[i]));
			ExpandForInner(it, &(branch->m_child->m_branch[i]));
		}
	}
	else
	{
		it->outerRelevant.push_back(branch);
	}
}

bool BIPSolver::FindIntersect(vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itF,
							  Slot* itC,
							  Point* intersect,
							  Point* nn,
							  Point* mid,
							  double* ratio,
							  bool* foundRegion)
{
	//first check orientation
	// itF is target, itC is reference
	int orientation = CheckOrientation((*itF)->m_rect, itC->branch->m_rect);
	// should not intersect or contained, contains, since itF is outer relevant nodes of itC
	assert(orientation != 0);
	if(!foundRegion[orientation-1])
	{
		// this is a valid f MBR, providing new information
		// itF is target, itC is reference
		GetTypicalRatio((*itF)->m_rect, itC->branch->m_rect, orientation, ratio);

		if(orientation % 2 == 0) // S, E, N, W
		{
			Point midPoint [2];
			pruningComputationCounter++;
			bool legal = GetMinMaxMidBtwn((*itF)->m_rect, itC->branch->m_rect, midPoint);
			assert(legal);
			if(orientation == 2)
			{
				mid[1] = midPoint[0];
				mid[2] = midPoint[1];
			}
			else if(orientation == 4)
			{
				mid[4] = midPoint[0];
				mid[5] = midPoint[1];
			}
			else if(orientation == 6)
			{
				mid[7] = midPoint[0];
				mid[8] = midPoint[1];
			}
			else if(orientation == 8)
			{
				mid[10] = midPoint[0];
				mid[11] = midPoint[1];
			}
		}
		else // SW, SE, NE, NW
		{
			Point midPoint;
			bool legal = GetMinMaxMidBtwn((*itF)->m_rect, itC->branch->m_rect , &midPoint);
			assert(legal);
			nn[orientation / 2] = midPoint;
			if(orientation == 1)
				mid[0] = midPoint;
			else if(orientation == 3)
				mid[3] = midPoint;
			else if(orientation == 5)
				mid[6] = midPoint;
			else if(orientation == 7)
				mid[9] = midPoint;
		}

		//calculate intersect
		switch(orientation)
		{
		case 1://SW
			if(foundRegion[1])//SW and S
				GetIntersect(mid[0],ratio[0], mid[2], ratio[3], intersect[0]);
			if(foundRegion[7])//SW and W
			{
				if(ratio[14] == 0)// remind ratio[14] might be +infinity
				{
					intersect[11].x = mid[10].x;
					intersect[11].y = ratio[0]*(intersect[11].x - mid[0].x) + mid[0].y;
				}
				else
					GetIntersect(mid[0], ratio[1], mid[10], ratio[14], intersect[11]);
			}
			break;
		case 2://S
			GetIntersect(mid[1], ratio[2], mid[2], ratio[3], intersect[1]);
			if(foundRegion[0])//SW and S
				GetIntersect(mid[0], ratio[0], mid[2], ratio[3], intersect[0]);
			if(foundRegion[2])//S and SE
				GetIntersect(mid[1], ratio[2], mid[3], ratio[5], intersect[2]);
			break;
		case 3://SE
			if(foundRegion[1])//S and SE
				GetIntersect(mid[1], ratio[2], mid[3], ratio[5], intersect[2]);
			if(foundRegion[3])//SE and E
			{
				if(ratio[6] == 0) // remind ratio[6] might be +infinity
				{
					intersect[3].x = mid[4].x;
					intersect[3].y = ratio[4]*(intersect[3].x - mid[3].x) + mid[3].y;
				}
				else
					GetIntersect(mid[3], ratio[4], mid[4], ratio[6], intersect[3]);
			}
			break;
		case 4://E
			if(ratio[7] == 0)
			{
				intersect[4].x = mid[4].x;
				intersect[4].y = ratio[7]*(intersect[4].x - mid[5].x)+mid[5].y;
			}
			if(ratio[6] == 0)
			{
				intersect[4].x = mid[5].x;
				intersect[4].y = ratio[6]*(intersect[4].x - mid[4].x)+mid[4].y;
			}
			else
				GetIntersect(mid[4], ratio[6], mid[5], ratio[7], intersect[4]);
			if(foundRegion[2])//SE and E
			{
				if(ratio[6] == 0)
				{
					intersect[3].x = mid[4].x;
					intersect[3].y = ratio[4]*(intersect[3].x - mid[3].x) + mid[3].y;
				}
				else
					GetIntersect(mid[3], ratio[4], mid[4], ratio[6], intersect[3]);
			}
			if(foundRegion[4])//E and NE
			{
				if(ratio[7] == 0)
				{
					intersect[5].x = mid[5].x;
					intersect[5].y = ratio[9]*(intersect[5].x - mid[6].x) + mid[6].y;
				}
				else
					GetIntersect(mid[5], ratio[7], mid[6], ratio[9], intersect[5]);
			}
			break;
		case 5://NE
			if(foundRegion[3])//NE and E
			{
				if(ratio[7] == 0)
				{
					intersect[5].x = mid[5].x;
					intersect[5].y = ratio[9]*(intersect[5].x - mid[6].x) + mid[6].y;				
				}
				else
					GetIntersect(mid[5], ratio[7], mid[6], ratio[9], intersect[5]);
			}
			if(foundRegion[5])//NE and N
			{
				GetIntersect(mid[6], ratio[8], mid[7], ratio[10], intersect[6]);
			}
			break;
		case 6:
			GetIntersect(mid[7], ratio[10], mid[8], ratio[11], intersect[7]);
			if(foundRegion[4])//NE and N
				GetIntersect(mid[6], ratio[8], mid[7], ratio[10], intersect[6]);
			if(foundRegion[6])//N and NW
				GetIntersect(mid[8], ratio[11], mid[9], ratio[13], intersect[8]);
			break;
		case 7:
			if(foundRegion[5])//N and NW
				GetIntersect(mid[9], ratio[13], mid[8], ratio[11], intersect[8]);
			if(foundRegion[7])//NW and W
			{
				if(ratio[15] == 0)
				{
					intersect[9].x = mid[11].x;
					intersect[9].y = ratio[12]*(intersect[9].x - mid[9].x) + mid[9].y;
				}
				else
					GetIntersect(mid[9], ratio[12], mid[11], ratio[15], intersect[9]);
			}
			break;
		case 8:
			if(ratio[15] == 0)
			{
				intersect[10].x = mid[11].x;
				intersect[10].y = ratio[14]*(intersect[10].x - mid[10].x) + mid[10].y;
			}
			if(ratio[14] == 0)
			{
				intersect[10].x = mid[10].x;
				intersect[10].y = ratio[15]*(intersect[10].x - mid[11].x) + mid[11].y;
			}
			GetIntersect(mid[10], ratio[14], mid[11], ratio[15], intersect[10]);
			if(foundRegion[6])//W and NW
			{	
				if(ratio[15] == 0)
				{
					intersect[9].x = mid[11].x;
					intersect[9].y = ratio[12]*(intersect[9].x - mid[9].x) + mid[9].y;
				}
				else
					GetIntersect(mid[11], ratio[15], mid[9], ratio[12], intersect[9]);
			}
			if(foundRegion[0])//W and SW
			{
				if(ratio[14] == 0)
				{
					intersect[11].x = mid[10].x;
					intersect[11].y = ratio[1]*(intersect[11].x - mid[0].x) + mid[0].y;
				}
				else
					GetIntersect(mid[10], ratio[14], mid[0], ratio[1], intersect[11]);
			}
			break;
		default:
			break;
		}
		foundRegion[orientation-1] = true;
		return true;
	}
	else
		return false;
}

bool BIPSolver::IsPruned(vector<RTree<unsigned int, double, DIMENSION, double>::Branch *>::iterator itF,
						 Slot * itC,
						 Point* intersect,
						 bool* foundRegion,
						 double* foundRatio)
{
	bool isPruned = false;
	// itF is target, itC is reference
	
	int orientation = CheckOrientation(itC->branch->m_rect, (*itF)->m_rect);		
	double ratio[16];
	GetTypicalRatio(itC->branch->m_rect, (*itF)->m_rect, orientation, ratio);
	Point midPoint[2];
	bool legal = GetMinMaxMidBtwn(itC->branch->m_rect, (*itF)->m_rect,  midPoint);
	assert(legal);
	bool lineA = false;
	bool lineB = false;
	double xmin = intersect[9].x < intersect[11].x ? intersect[9].x : intersect[11].x;
	double xmax = intersect[3].x > intersect[5].x ? intersect[3].x : intersect[5].x;
	double ymin = intersect[0].y < intersect[2].y ? intersect[0].y : intersect[2].y;
	double ymax = intersect[6].y > intersect[8].y ? intersect[6].y : intersect[8].y;
	xmin = xmin < (-1)*RANGE ? (-1)*RANGE : xmin;
	ymin = ymin < (-1)*RANGE ? (-1)*RANGE : ymin;
	xmax = xmax > RANGE ? RANGE : xmax;
	ymax = ymax > RANGE ? RANGE : ymax;

	switch(orientation)
	{
	case 5:// F in in SW of C
		// check S and W
		if(foundRegion[1] && foundRegion[7] && foundRegion[2] && foundRegion[6]
		&& midPoint[0].x <= xmin && midPoint[0].y <= ymin)
			isPruned = true;
		else if(foundRegion[1] && foundRegion[7])
		{
			lineA = CheckPointAndLine(midPoint[0], ratio[0], intersect[11]);
			lineB = CheckPointAndLine(midPoint[0], ratio[1], intersect[0]);
			// ratio[1] checks S
			if(foundRegion[2])
				lineA = lineA && CheckPointAndLine(midPoint[0], ratio[1], intersect[2]);
			else
				lineA = lineA && ratio[1] <= foundRatio[2];// negative, steeper
			// ratio[0] checks W
			if(foundRegion[6])
				lineB = lineB && CheckPointAndLine(midPoint[0], ratio[0], intersect[9]);
			else
			{
				if(foundRatio[15] == 0)
					lineB = lineB && (ratio[0] <= 0);
				else
					lineB = lineB && ratio[0] >= foundRatio[15];// negative, more subdued
			}
			if(lineA && lineB)
			isPruned = true;
		}
		break;
	
	case 6:// F is in the S of C
		// check S, W and E
		if( foundRegion[0] && foundRegion[2] && foundRegion[7] && foundRegion[3])
		{
			lineA = CheckPointAndLine(midPoint[1], ratio[3], intersect[0])
				&& CheckPointAndLine(midPoint[1], ratio[3], intersect[2])
				&& CheckPointAndLine(midPoint[1], ratio[3], intersect[11]);
			lineB = CheckPointAndLine(midPoint[0], ratio[2], intersect[0])
				&& CheckPointAndLine(midPoint[0], ratio[2], intersect[2])
				&& CheckPointAndLine(midPoint[0], ratio[2], intersect[3]);
			if(foundRegion[6])
				lineA = lineA && CheckPointAndLine(midPoint[1], ratio[3], intersect[9]);
			else 
			{
				if(foundRatio[15] == 0)
					lineA = lineA && (ratio[3] <= 0);
				else
					lineA = lineA && (ratio[3] >= foundRatio[15]);
			}
			if(foundRegion[4])
				lineB = lineB && CheckPointAndLine(midPoint[0], ratio[2], intersect[5]);
			else
			{
				if(foundRatio[7] == 0)
					lineB = lineB && (ratio[2] >= 0);
				else
					lineB = lineB && (ratio[2] <= foundRatio[7]);
			}
			if(lineA && lineB)
				isPruned = true;
		}
		break;
	
	case 7:// F is in the SE of C
		//check S and E
		if(foundRegion[1] && foundRegion[3] && foundRegion[0] && foundRegion[4]
		&& midPoint[0].x >= xmax && midPoint[0].y <= ymin)
			isPruned = true;
		else if(foundRegion[1] && foundRegion[3])
		{
			// ratio[4] checks S
			lineA = CheckPointAndLine(midPoint[0], ratio[5], intersect[2]);			
			// ratio[5] check E
			lineB = CheckPointAndLine(midPoint[0], ratio[4], intersect[3]);

			if(foundRegion[0])
				lineA = lineA && CheckPointAndLine(midPoint[0], ratio[5], intersect[0]);
			else
				lineA = lineA && (ratio[5] >= foundRatio[3]);
			if(foundRegion[4])
				lineB = lineB && CheckPointAndLine(midPoint[0], ratio[4], intersect[5]);
			else
			{
				if(foundRatio[7] == 0)
					lineB = lineB && (ratio[4] >= 0);
				else
					lineB = lineB && (ratio[4] <= foundRatio[7]);
			}
			if(lineA && lineB)
				isPruned = true;
		}
		break;
		
	case 8://F is in the E of C
		//checks E, N and S
		if(foundRegion[2] && foundRegion[4] && foundRegion[1] && foundRegion[5])
		{
			if(ratio[6] == 0)
			{
				lineA = midPoint[0].x >= xmax;
			}			
			else
			{
				// midPoint[0] and ratio[7] checks E and S
				lineA = CheckPointAndLine(midPoint[0], ratio[6], intersect[3])
					&& CheckPointAndLine(midPoint[0], ratio[6], intersect[5])
					&& CheckPointAndLine(midPoint[0], ratio[6], intersect[2]);
			}
			if(ratio[7] == 0)
			{
				lineB = midPoint[1].x >= xmax;
			}
			else
			{
				// midPoint[1] and ratio[6] checks E and N
				lineB = CheckPointAndLine(midPoint[1], ratio[7], intersect[3])
					|| CheckPointAndLine(midPoint[1], ratio[7], intersect[5])
					|| CheckPointAndLine(midPoint[1], ratio[7], intersect[6]);
			}
			if(foundRegion[0])
				lineA = lineA && CheckPointAndLine(midPoint[0], ratio[6], intersect[0]);
			else
				lineA = lineA && (ratio[6] >= foundRatio[3]);
			if(foundRegion[6])
				lineB = lineB || CheckPointAndLine(midPoint[1], ratio[7], intersect[8]);
			else
				lineB = lineB || (ratio[7] > foundRatio[11]);
			if(lineA && (!lineB))
				isPruned = true;
		}
		break;
		
	case 1://F is in the NE of C
		//checks E and N
		if(foundRegion[3] && foundRegion[5] && foundRegion[2] && foundRegion[6]
		&& midPoint[0].x >= xmax && midPoint[0].y >= ymax)
			isPruned = true;
		else if(foundRegion[3] && foundRegion[5])
		{
			// ratio[8] checks E
			lineA = CheckPointAndLine(midPoint[0], ratio[8], intersect[5]);
			// ratio[9] checks N 
			lineB = CheckPointAndLine(midPoint[0], ratio[9], intersect[6]);
			if(foundRegion[2])
				lineA = lineA || CheckPointAndLine(midPoint[0], ratio[8], intersect[3]);
			else
			{
				if(foundRatio[6] == 0)
					lineA = lineA || (ratio[8] > 0);
				else
					lineA = lineA || (ratio[8] < foundRatio[6]); // ratio[8] > foundRatio[6]
			}
			if(foundRegion[6])
				lineB = lineB || CheckPointAndLine(midPoint[0], ratio[9], intersect[8]);
			else
				lineB = lineB || (ratio[9] > foundRatio[11]); //ratio[9] < foundRatio[11]
			if((!lineA) && (!lineB))
				isPruned = true;
		}
		break;
		
	case 2://F is in the N of C
		//checks N, E and W
		if(foundRegion[4] && foundRegion[6] && foundRegion[3] && foundRegion[7])
		{
			// midPoint[0] and ratio[10] checks N and W
			lineA = CheckPointAndLine(midPoint[0], ratio[11], intersect[6])
				|| CheckPointAndLine(midPoint[0], ratio[11], intersect[8])
				|| CheckPointAndLine(midPoint[0], ratio[11], intersect[5]);
			// midPoint[1] and ratio[11] checks N and E
			lineB = CheckPointAndLine(midPoint[1], ratio[10], intersect[6])
				|| CheckPointAndLine(midPoint[1], ratio[10], intersect[8])
				|| CheckPointAndLine(midPoint[1], ratio[10], intersect[9]);
			if(foundRegion[0])
				lineA = lineA || CheckPointAndLine(midPoint[1], ratio[10], intersect[11]);
			else
				lineA = lineA || ratio[10] > foundRatio[14];//ratio[10] < foundRatio[14]
			if(foundRegion[2])
				lineB = lineB || CheckPointAndLine(midPoint[0], ratio[11], intersect[3]);
			else
			{
				if(foundRatio[6] == 0)
					lineB = lineB || (ratio[11] > 0);
				else
					lineB = lineB || ratio[11] < foundRatio[6];//ratio[11] > foundRatio[6]
			}
			if((!lineA) && (!lineB))
				isPruned = true;
		}
		break;
		
	case 3:// F is on the NW of C
		//checks N and W
		if(foundRegion[5] && foundRegion[7] && foundRegion[4] && foundRegion[0]
		&& midPoint[0].x <= xmin && midPoint[0].y >= ymax)
			isPruned = true;
		else if(foundRegion[5] && foundRegion[7])
		{
			// ratio[12] checks N
			lineA = CheckPointAndLine(midPoint[0], ratio[12], intersect[8]);
			// ratio[13] checks W
			lineB = CheckPointAndLine(midPoint[0], ratio[13], intersect[9]);

			if(foundRegion[4])
				lineA = lineA || CheckPointAndLine(midPoint[0], ratio[12], intersect[6]);
			else
				lineA = lineA || ratio[12] < foundRatio[10]; //ratio[12] > foundRatio[10]
			if(foundRegion[0])
				lineB = lineB || CheckPointAndLine(midPoint[0], ratio[13], intersect[11]);
			else
			{
				if(foundRatio[14] == 0)
					lineB = lineB || ratio[13] < 0;
				else
					lineB = lineB || ratio[13] > foundRatio[14];//ratio[13] < foundRatio[14];
			}
			if((!lineA) && (!lineB))
				isPruned = true;
		}
		break;
	case 4:// F is on the W of C
		// checks W, S and N
		if(foundRegion[6] && foundRegion[0] && foundRegion[5] && foundRegion[1])
		{
			if(ratio[14] == 0)
			{
				lineA = midPoint[0].x <= xmin;
			}
			else
			{
				lineA = CheckPointAndLine(midPoint[0], ratio[14], intersect[9])
					&& CheckPointAndLine(midPoint[0], ratio[14], intersect[11])
					&& CheckPointAndLine(midPoint[0], ratio[14], intersect[0]);
			}
			if(ratio[15] == 0)
			{
				lineB = midPoint[1].x <= xmin;
			}
			else
			{	lineB = CheckPointAndLine(midPoint[1], ratio[15], intersect[9])
					|| CheckPointAndLine(midPoint[1], ratio[15], intersect[11])
					|| CheckPointAndLine(midPoint[1], ratio[15], intersect[8]);
			}	
			if(foundRegion[2])
				lineA = lineA && CheckPointAndLine(midPoint[0], ratio[14], intersect[2]);
			else
				lineA = lineA && ratio[14] < foundRatio[2];//ratio[14] > foundRatio[10]
			if(foundRegion[4]) 
				lineB = lineB || CheckPointAndLine(midPoint[1], ratio[15], intersect[6]);
			else
				lineB = lineB || ratio[15] < foundRatio[10];
			if(lineA && (!lineB))
				isPruned = true;
		}
		break;
	default:
		break;
	}
	
	return isPruned;
}

int BIPSolver::EstimateUpperBound(Slot* itC, Point* intersect, Point * nn, bool* foundRegion)
{
	double xmin, xmax, ymin, ymax;
	if(foundRegion[6] && foundRegion[7] && foundRegion[0])
		xmin = intersect[9].x > intersect[11].x ? intersect[11].x : intersect[9].x;
	else
		xmin = (-1) * RANGE;
	if(foundRegion[0] && foundRegion[1] && foundRegion[2])
		ymin = intersect[0].y > intersect[2].y ? intersect[2].y : intersect[0].y;
	else
		ymin = (-1) * RANGE;
	if(foundRegion[2] && foundRegion[3] && foundRegion[4])
		xmax = intersect[3].x > intersect[5].x ? intersect[3].x : intersect[5].x;
	else
		xmax = RANGE;
	if(foundRegion[4] && foundRegion[5] && foundRegion[6])
		ymax = intersect[6].y > intersect[8].y ? intersect[6].y : intersect[8].y;
	else
		ymax = RANGE;
	xmax = xmax > RANGE ? RANGE : xmax;
	ymax = ymax > RANGE ? RANGE : ymax;
	xmin = xmin < RANGE *(-1) ? RANGE*(-1) : xmin;
	ymin = ymin < RANGE * (-1) ? RANGE * (-1) : ymin;
	// store the rect in the rectM
	itC->rectM.m_min[0] = xmin;
	itC->rectM.m_min[1] = ymin;
	itC->rectM.m_max[0] = xmax;
	itC->rectM.m_max[1] = ymax;

	// estimate upper bound, notice we first perform pruning
	RTree<unsigned int, double, DIMENSION, double>::Rect rect;
	rect.m_min[0] = xmin;
	rect.m_min[1] = ymin;
	rect.m_max[0] = xmax;
	rect.m_max[1] = ymax;

	// a rectangle overlapped by 4 other rects on the corner, can split into 5 sub rects
	vector<RTree<unsigned int, double, DIMENSION, double>::Rect> refinedRect;
	refinedRect.push_back(rect);
	RTree<unsigned int, double, DIMENSION, double>::Rect domainanceRect[4];
	if(foundRegion[0])
	{
		domainanceRect[0].m_min[0] = (-1) * RANGE;
		domainanceRect[0].m_min[1] = (-1) * RANGE;
		domainanceRect[0].m_max[0] = nn[0].x;
		domainanceRect[0].m_max[1] = nn[0].y;
		vector<RTree<unsigned int, double, DIMENSION, double>::Rect> answer;
		SplitRect(&refinedRect, &domainanceRect[0], answer);
		refinedRect = answer;
	}
	if(foundRegion[2])
	{
		domainanceRect[1].m_min[0] = nn[1].x;
		domainanceRect[1].m_min[1] = (-1) * RANGE;
		domainanceRect[1].m_max[0] = RANGE;
		domainanceRect[1].m_max[1] = nn[1].y;
		vector<RTree<unsigned int, double, DIMENSION, double>::Rect> answer;
		SplitRect(&refinedRect, &domainanceRect[1], answer);
		refinedRect = answer;
	}
	if(foundRegion[4])
	{
		domainanceRect[2].m_min[0] = nn[2].x;
		domainanceRect[2].m_min[1] = nn[2].y;	
		domainanceRect[2].m_max[0] = RANGE;
		domainanceRect[2].m_max[1] = RANGE;
		vector<RTree<unsigned int, double, DIMENSION, double>::Rect> answer;
		SplitRect(&refinedRect, &domainanceRect[2], answer);
		refinedRect = answer;
	}
	if(foundRegion[6])
	{
		domainanceRect[3].m_min[0] = (-1) * RANGE;
		domainanceRect[3].m_min[1] = nn[3].y;
		domainanceRect[3].m_max[0] = nn[3].x;
		domainanceRect[3].m_max[1] = RANGE;
		vector<RTree<unsigned int, double, DIMENSION, double>::Rect> answer;
		SplitRect(&refinedRect, &domainanceRect[3], answer);
		refinedRect = answer;
	}

	int inf = 0;
	for(vector<RTree<unsigned int, double, DIMENSION, double>::Rect>::iterator itR = refinedRect.begin();
		itR != refinedRect.end();
		itR++)
	{
		int access = 0;
		inf += treeM.SooneAggregateCount(*itR, rootM, access);
		itC->finalRectM.push_back((*itR));
		//count the access
		countAccessM += access;
		counterIOM += access;
	}
	itC->infUpperBound = inf;
	return inf;
}