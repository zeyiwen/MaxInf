/**************************************************************/
#ifndef SORT_H
#define SORT_H
#include "Utility.h"

class Sort
{
public:
	bool InsertionSort(list<EntryC*> &Q_c);
};

bool Sort::InsertionSort(list<EntryC*> &Q_c)
{
	list<EntryC*>::iterator first, it_current, pos; 
	int nSize = Q_c.size();
	first = Q_c.begin();
	for(int i = 0; i < nSize; i++)
	{
		it_current = first;
		pos = first;
		first++;
		for(int j = i; j > 0; j--)
		{
			pos--;
			if(j == 1)
			{
				if(*(*it_current) <= *(*pos))
				{
					if(j == i)
					{
						break;
					}
					pos++;
					Q_c.insert(pos,*it_current);
					Q_c.erase(it_current);
					break;
				}
				else
				{
					Q_c.push_front(*it_current);
					Q_c.erase(it_current);
					break;
				}
			}
			if(*(*it_current) <= *(*pos))
			{
				if(i == j)
				{
					break;
				}
				else
				{
					pos++;
					Q_c.insert(pos,*it_current);
					Q_c.erase(it_current);
					break;
				}
			}
		}
	}
	return true;
}

#endif SORT_H
/**********************************************************/