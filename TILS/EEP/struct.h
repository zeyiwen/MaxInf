#ifndef STRUCT_H
#define STRUCT_H
/* This file define all the data structure used by EEP */

#include "Utility.h"
#include "constant.h"
#include <list>
using namespace std;

//define namespace
namespace EEP_SPACE
{


class aNode; //node for aRtree
class Node;

class Rect
{
public:
	Rect(){}
	~Rect(){}
	Rect(const Rect & r)
	{
		for(int i = 0; i < 2; i++)
		{
			m_min[i] = r.m_min[i];
			m_max[i] = r.m_max[i];
		} 
	}
	float m_min[2];                      ///< Min dimensions of bounding box 
	float m_max[2];                      ///< Max dimensions of bounding box 
};

class Branch
{
public:
	Branch(){}
	~Branch(){}
	Branch(const Branch &A)
	{
		m_rect = A.m_rect;
		m_child = A.m_child;
		m_data = A.m_data;
	}
	Rect m_rect;                                  ///< Bounds
	union
	{
		Node* m_child;                              ///< Child node
		int m_data;                            ///< Data Id or Ptr
	};
};

class aBranch //branch for aRtree
{
public:
	aBranch(){}
	~aBranch(){}
	aBranch(const aBranch &A)
	{
		if(A.amount_of_nodes < 0)
			amount_of_nodes = 0;
		else
			amount_of_nodes = A.amount_of_nodes;
		m_rect = A.m_rect;
		m_child = A.m_child;
		m_data = A.m_data;
	}
	Rect m_rect;                                  ///< Bounds
	int amount_of_nodes;
	union
	{
		aNode* m_child;                              ///< Child node
		int m_data;                            ///< Data Id or Ptr
	};
	
};

//define entries used by EEP. These entries are keep in the three lists
class EntryM;

class EntryC
{
public:
	EntryC(){}
	~EntryC(){}
	EntryC(const EntryC &A)
	{
		minInf = A.minInf;
		maxInf = A.maxInf;
		L_m = A.L_m;
		pC = A.pC;
	}
	Branch *pC;
	int minInf;
	int maxInf;
	list<EntryM*> L_m; 
	friend bool operator > (EntryC &A, EntryC& B)
	{
		if(A.maxInf > B.maxInf)
			return true;
		else if(A.maxInf == B.maxInf)
		{
			if(A.minInf > B.minInf)
				return true;
			else
				return false;
		}
		else
			return false;
	}

	friend bool operator >= (EntryC &A, EntryC& B)
	{
		if(A.maxInf > B.maxInf)
			return true;
		else if(A.maxInf == B.maxInf)
		{
			if(A.minInf >= B.minInf)
				return true;
			else
				return false;
		}
		else
			return false;
	}

	friend bool operator < (EntryC &A, EntryC &B)
	{
		if(A.maxInf < B.maxInf)
			return true;
		else if(A.maxInf == B.maxInf)
		{
			if(A.minInf < B.minInf)
				return true;
			else 
				return false;
		}
		else
			return false;
	}
	friend bool operator <= (EntryC& A, EntryC& B)
	{

		if(A.maxInf < B.maxInf)
			return true;
		else if(A.maxInf == B.maxInf)
		{
			if(A.minInf <= B.minInf)
				return true;
			else
				return false;
		}
		else
			return false;
	}
};
class EntryF
{
public:
	EntryF(){}
	~EntryF(){}
	EntryF(const EntryF &A)
	{
		L_m = A.L_m;
		pF = A.pF;
		square = A.square;
	}
	list<EntryM*> L_m;
	float square;
	Branch *pF;
};
class EntryM
{
public:
	EntryM(){}
	~EntryM(){}
	EntryM(const EntryM &A)
	{
		if(A.amount_of_nodes < 0)
			amount_of_nodes = 0;
		else
			amount_of_nodes = A.amount_of_nodes;
		min_dist = A.min_dist;
		minExist_dist = A.minExist_dist;
		L_c = A.L_c;
		L_f = A.L_f;
		square = A.square;
		pM = A.pM;
	}
	float min_dist;
	float minExist_dist;
	float square;
	int amount_of_nodes;
	list<EntryF*> L_f;
	list<EntryC*> L_c;
	aBranch *pM;
};

class ListBranch
{
public:
	ListBranch(){}
	ListBranch(const ListBranch& N)
	{
		m_next = N.m_next;
		m_branch = N.m_branch;
	}
	ListBranch* m_next;                             ///< Next in list
	Branch* m_branch;                                 ///< Node
};


/// Node for each Rtree branch level in EEP
class Node
{
public:
	Node(){}
	Node(Node const & N)
	{
		m_level = N.m_level;
		m_count = N.m_count;
		for(int i = 0; i < FANOUT; i++)
		{
			m_branch[i] = N.m_branch[i];
		}
	}
	bool IsInternalNode()                         { return (m_level > 0); } // Not a leaf, but a internal node
	bool IsLeaf()                                 { return (m_level == 0); } // A leaf, contains data

	int m_count;                                  ///< Count
	int m_level;                                  ///< Leaf is zero, others positive
	Branch m_branch[FANOUT];                    ///< Branch
};

/// Node for each aRtree branch level
class aNode
{
public:
	aNode(){}
	aNode(aNode const & N)
	{
		m_level = N.m_level;
		m_count = N.m_count;
		total_count = N.total_count;
		for(int i = 0; i < FANOUT; i++)
		{
			m_branch[i] = N.m_branch[i];
		}
	}
	bool IsInternalNode()                         { return (m_level > 0); } // Not a leaf, but a internal node
	bool IsLeaf()                                 { return (m_level == 0); } // A leaf, contains data

	int m_count;                                  ///< Count
	int m_level;                                  ///< Leaf is zero, others positive
	int total_count;
	aBranch m_branch[FANOUT];                    ///< Branch
};

class ListNode
{
public:
	ListNode(){}
	ListNode(const ListNode& N)
	{
		m_next = N.m_next;
		m_node = N.m_node;
	}
	ListNode* m_next;                             ///< Next in list
	Node* m_node;                                 ///< Node
};

struct ListPoint
{
	Point mPoint;
	ListPoint *next;
};

}//end EEPSPACE

#endif STRUCT_H