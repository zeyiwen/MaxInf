#ifndef INFLUENCE_H
#define INFLUENCE_H

#include <queue>
#include <assert.h>
#include "Utility.h"

//comparator
bool C_less(EntryC*&elem1, EntryC*&elem2)
{
	return (elem1->maxInf <= elem2->maxInf);
}

//expand list C to at least k entries
bool Expand(list<EntryC*>& Q_c2, int& k)
{
	EntryC *pTempEnC, *pEnC;
	pEnC = new EntryC;
	list<EntryC*> Qtemp2;
	int nSize = 0, counter = 0, nEntry = 0;
	nSize = Q_c2.size();
	int number_of_child = 0;
	for(int i = 0; nEntry < k && i < nSize; i++)
	{	
		pTempEnC = Q_c2.back();
		Q_c2.pop_back();
		if(pTempEnC->pC->m_child == NULL)
		{
			Qtemp2.push_back(pTempEnC);
			continue;
		}
		number_of_child = pTempEnC->pC->m_child->m_count;
		nEntry += number_of_child;

		for(int i = 0; i < number_of_child; i++)
		{
			if(pTempEnC->pC->m_child->m_level == 0)
				pTempEnC->pC->m_child->m_branch[i].m_child = NULL;
			pEnC->pC = &(pTempEnC->pC->m_child->m_branch[i]);
			
			pEnC->L_m = pTempEnC->L_m;
			pEnC->maxInf = pTempEnC->maxInf;
			pEnC->minInf = pTempEnC->minInf;
			Qtemp2.push_back(pEnC);
			pEnC = new EntryC;
		}
	}
	while(nEntry < k && !Qtemp2.empty())
	{
		pTempEnC = Qtemp2.front();
		Qtemp2.pop_front();
		if(pTempEnC->pC->m_child == NULL)
		{
			Q_c2.push_back(pTempEnC);
			continue;
		}
		number_of_child = pTempEnC->pC->m_child->m_count;
		nEntry += number_of_child;

		for(int i = 0; i < number_of_child; i++)
		{
			if(pTempEnC->pC->m_child->m_level == 0)
				pTempEnC->pC->m_child->m_branch[i].m_child = NULL;
			pEnC->maxInf = pTempEnC->maxInf;
			pEnC->minInf = pTempEnC->minInf;
			pEnC->pC = &(pTempEnC->pC->m_child->m_branch[i]);
			Q_c2.push_back(pEnC);
			pEnC = new EntryC;
		}
	}
	while(!Qtemp2.empty())
	{
		pTempEnC = Qtemp2.front();
		Qtemp2.pop_front();
		Q_c2.push_back(pTempEnC);
	}
	return 1;
}

//update c in the prepration stage
bool UpdateC(list<EntryC*> &Q_c2, list<EntryF*> &Q_f, list<EntryM*> &Q_m, long nMa)
{
	int number_of_m = 0, number_of_f = 0, number_of_c = 0;
	float mindist_F_M = 0, minExistDist_F_M = 0, mindist_C_M = 0, maxdist_C_M = 0, temp_mindist = 0, temp_minExistDist = 0;
	int id = 0;
	DistCompute help;
	EntryC *pEnC, *pTempEnC;

	number_of_m = Q_m.size();

	list<EntryC*> temp_Qc2;

	number_of_c = Q_c2.size();
	number_of_f = Q_f.size();

	list<EntryM*>::iterator it_m2;
	list<EntryF*>::iterator it_f2;
	it_m2 = Q_m.begin();

	for(int i = 0; i < number_of_m; i++)
	{
		//initialize L_m, minInf and maxInf
		for(int k = 0; k < number_of_c; k++)
		{
			pTempEnC = Q_c2.back();
			pTempEnC->maxInf = nMa;
			pTempEnC->minInf = 0;
			pTempEnC->L_m.push_back(*it_m2);

			temp_Qc2.push_back(pTempEnC);
			(*it_m2)->L_c.push_back(pTempEnC);
			Q_c2.pop_back();
		}
		while(!temp_Qc2.empty())
		{
			pEnC = temp_Qc2.back();
			temp_Qc2.pop_back();
			Q_c2.push_back(pEnC);
		}
		//initialize L_f and f's L_m
		
		it_f2 = Q_f.begin();
		(*it_f2)->L_m.push_back(*it_m2);//push m into f's L_m

		mindist_F_M = help.minDist((*it_m2)->pM->m_rect, (*it_f2)->pF->m_rect);
		minExistDist_F_M = help.minExistDist((*it_m2)->pM->m_rect, (*it_f2)->pF->m_rect);
		(*it_m2)->L_f.push_back(*it_f2);
		it_f2++;
		for(int j = 1; j < number_of_f; j++)
		{
			(*it_m2)->L_f.push_back(*it_f2);
			(*it_f2)->L_m.push_back(*it_m2);//push m into f's L_m

			temp_mindist = help.minDist((*it_m2)->pM->m_rect, (*it_f2)->pF->m_rect);
			temp_minExistDist = help.minExistDist((*it_m2)->pM->m_rect, (*it_f2)->pF->m_rect);
			if(temp_mindist < mindist_F_M)
			{
				mindist_F_M = temp_mindist;
			}
			if(temp_minExistDist < minExistDist_F_M)
			{
				minExistDist_F_M = temp_minExistDist;
			}
			it_f2++;
		}
		(*it_m2)->min_dist = mindist_F_M;
		(*it_m2)->minExist_dist = minExistDist_F_M;
		it_m2++;
	}
	return true;
}

//update c caused by expanding imp(m)
int UpdateC_M(EntryC*& pEnC, EntryM*& pEnM)
{
	int nIsUpdate = 0;
	float mindist_C_M = 0, maxdist_C_M = 0;
	DistCompute help;
	int nSizeofLc;
	list<EntryC*>::iterator it_c2;
	//initialize L_m, minInf and maxInf
	mindist_C_M = help.minDist(pEnC->pC->m_rect, pEnM->pM->m_rect);
	maxdist_C_M = help.maxDist(pEnC->pC->m_rect, pEnM->pM->m_rect);
	//m is only affected by C
	if(pEnM->min_dist >= maxdist_C_M) // m's min_dist and max_dist are distances from L_f to m
	{
		//minus m to c's minInf
		pEnC->minInf += pEnM->amount_of_nodes;
		nIsUpdate = 1;
		it_c2 = pEnM->L_c.begin();
		nSizeofLc = pEnM->L_c.size();
		for(int j = 0; j < nSizeofLc; j++)
		{
			if((*it_c2) == pEnC)
			{
				pEnM->L_c.erase(it_c2);
				break;
			}
			else
			{
				it_c2++;
			}
		}
	}
	//m is only affected by L_f
	if(pEnM->minExist_dist < mindist_C_M)
	{
		//remove m from Lm and minus c's maxInf by m
		//minus m to c's maxInf
		pEnC->maxInf -= pEnM->amount_of_nodes;
		nIsUpdate = 2;
		it_c2 = pEnM->L_c.begin();
		nSizeofLc = pEnM->L_c.size();
		for(int j = 0; j < nSizeofLc; j++)
		{
			if((*it_c2) == pEnC)
			{
				pEnM->L_c.erase(it_c2);
				break;
			}
			else
			{
				it_c2++;
			}
		}
	}
	return nIsUpdate;
}

//update c caused by expanding imp(f)
//Return value: not updated return 0, update minInf return 1, update maxInf return 2
int UpdateC_M2(EntryC*& pEnC, EntryM*& pEnM)
{
	int nIsUpdate = 0, nRemove_m;
	float mindist_C_M = 0, maxdist_C_M = 0;
	DistCompute help;

	list<EntryC*>::iterator it_c2;
	list<EntryM*>::iterator tempIt_m2;
	int nSizeofLc;
	//initialize L_m, minInf and maxInf
	mindist_C_M = help.minDist(pEnC->pC->m_rect, pEnM->pM->m_rect);
	maxdist_C_M = help.maxDist(pEnC->pC->m_rect, pEnM->pM->m_rect);
	//m is only affected by C
	if(pEnM->min_dist >= maxdist_C_M) // m's min_dist and max_dist are distances from L_f to m
	{
		//minus m to c's minInf
		pEnC->minInf += pEnM->amount_of_nodes;
		nIsUpdate = 1;


		nRemove_m = pEnC->L_m.size();
		tempIt_m2 = pEnC->L_m.begin();
		for(int j = 0; j < nRemove_m; j++)
		{
			if((*tempIt_m2) == pEnM)
			{
				pEnC->L_m.erase(tempIt_m2);
				break;
			}
			else
				tempIt_m2++;
		}

		it_c2 = pEnM->L_c.begin();
		nSizeofLc = pEnM->L_c.size();
		for(int j = 0; j < nSizeofLc; j++)
		{
			if((*it_c2) == pEnC)
			{
				pEnM->L_c.erase(it_c2);
				break;
			}
			else
			{
				it_c2++;
			}
		}
	}
	//m is only affected by L_f
	if(pEnM->minExist_dist < mindist_C_M)
	{
		//remove m from Lm and minus c's maxInf by m
		//minus m to c's maxInf
		pEnC->maxInf -= pEnM->amount_of_nodes;
		nIsUpdate = 2;

		nRemove_m = pEnC->L_m.size();
		tempIt_m2 = pEnC->L_m.begin();
		for(int j = 0; j < nRemove_m; j++)
		{
			if((*tempIt_m2) == pEnM)
			{
				pEnC->L_m.erase(tempIt_m2);
				break;
			}
			else
				tempIt_m2++;
		}

		it_c2 = pEnM->L_c.begin();
		nSizeofLc = pEnM->L_c.size();
		for(int j = 0; j < nSizeofLc; j++)
		{
			if((*it_c2) == pEnC)
			{
				pEnM->L_c.erase(it_c2);
				break;
			}
			else
			{
				it_c2++;
			}
		}
	}
	return nIsUpdate;
}


bool getImp_m(list<EntryM*> &Q_m, EntryM*& imp_m2)
{
	if(Q_m.empty())
	{
		return false;
	}
	int number_of_m, nSize_of_Lf, nSinze_of_Lm, nAffect_c;
	int index = -1;
	list<EntryM*>::iterator it_imp_m2, it_m2, tempIt_m2, it_Lm2;
	list<EntryF*>::iterator it_f2;
	it_m2 = Q_m.begin();
	//remove m with empty S_c
	while((*it_m2)->L_c.empty())
	{
		//erase it_m from f's L_m
		nSize_of_Lf = (*it_m2)->L_f.size();
		it_f2 = (*it_m2)->L_f.begin();
		for(int j = 0; j < nSize_of_Lf; j++)
		{
			nSinze_of_Lm = (*it_f2)->L_m.size();
			it_Lm2 = (*it_f2)->L_m.begin();
			for(int k = 0; k < nSinze_of_Lm; k++)
			{
				if((*it_m2) == (*it_Lm2))
				{
					(*it_f2)->L_m.erase(it_Lm2);
					break;
				}
				else
					it_Lm2++;
			}
			it_f2++;
		}
		//end erase it_m from f's L_m
		Q_m.erase(it_m2);
		it_m2 = Q_m.begin();
		if(Q_m.empty())
			return false;
	}
	number_of_m = Q_m.size();
	it_imp_m2 = it_m2;
	double imp_value = 0, imp_temp = 0;
	nSize_of_Lf = (*it_m2)->L_f.size();
	nAffect_c = (*it_m2)->L_c.size();
	imp_value = (nAffect_c * (*it_m2)->amount_of_nodes * nSize_of_Lf * (*it_m2)->square);
	//find imp(m)
	for(int i = 1; i < number_of_m; i++)
	{
		tempIt_m2 = it_m2;
		it_m2++;
		//remove m with empty S_c
		if ((*it_m2)->L_c.empty())
		{
			//erase it_m from f's L_m
			nSize_of_Lf = (*it_m2)->L_f.size();
			it_f2 = (*it_m2)->L_f.begin();
			for(int j = 0; j < nSize_of_Lf; j++)
			{
				nSinze_of_Lm = (*it_f2)->L_m.size();
				it_Lm2 = (*it_f2)->L_m.begin();
				for(int k = 0; k < nSinze_of_Lm; k++)
				{
					if((*it_m2) == (*it_Lm2))
					{
						(*it_f2)->L_m.erase(it_Lm2);
						break;
					}
					else
						it_Lm2++;
				}
				it_f2++;
			}
	
			(*it_m2)->L_f.clear();
			Q_m.erase(it_m2);
			it_m2 = tempIt_m2;
			continue;
		}
		
		nSize_of_Lf = (*it_m2)->L_f.size();
		nAffect_c = (*it_m2)->L_c.size();
		imp_temp = (nAffect_c * (*it_m2)->amount_of_nodes * nSize_of_Lf * (*it_m2)->square);
		if(imp_temp > imp_value)
		{
			imp_value = imp_temp;
			it_imp_m2 = it_m2;
		}
	}
	imp_m2 = *it_imp_m2;
	Q_m.erase(it_imp_m2);
	return true;
}
//get imp(f)
bool getImp_f(list<EntryF*> &Q_f, EntryF*& imp_f2)
{
	if(Q_f.empty())
	{
		return false;
	}
	int number_of_f, number_of_m;
	list<EntryF*>::iterator it2, tempIt_f2, it_f2;
	it2 = Q_f.begin();
	//remove f with empty S_m
	while((*it2)->L_m.empty())
	{
		Q_f.erase(it2);
		if(Q_f.empty())
		{
			return false;
		}
		it2 = Q_f.begin();
	}

	tempIt_f2 = it2;
	number_of_f = Q_f.size();
	number_of_m = (*it2)->L_m.size();

	double imp_value = 0, imp_temp = 0;
	imp_value = (number_of_m * (*it2)->square);
	//find imp(f)
	for(int i = 1; i < number_of_f; i++)
	{
		it2++;
		//remove f with empty S_m
		while((*it2)->L_m.empty())
		{
			it_f2 = it2;
			it_f2--;
			Q_f.erase(it2);
			it2 = it_f2;
			if(Q_f.empty())
			{
				return false;
			}
		}
		
		number_of_m = (*it2)->L_m.size();
		imp_temp = (number_of_m * (*it2)->square);
		if(imp_temp > imp_value)
		{
			imp_value = imp_temp;
			tempIt_f2 = it2;
		}
	}
	if(imp_value == 0)
		return false;
	it2 = tempIt_f2;
	imp_f2 = *it2;
	Q_f.erase(it2);
	return true;
}

//input list L_f, expand its imp(f)
int Expand_F(list<EntryF*> &L_f)
{
	EntryF *imp_f2;
	EntryF *pTempEnF;
	imp_f2 = new EntryF;
	pTempEnF = new EntryF;
	list<EntryF*> tempL_f;
	list<EntryF*>::iterator it_f2, tempIt_f2, tempIt_f_remove2, tempIt_f_for_tempLf;
	list<EntryM*>::iterator it_m2, tempIt_m2;
	list<EntryC*>::iterator tempIt_c2, it_c2;
	float mindist_F_M, mindist_imp_f_M, temp_mindist, minExist_dist/*, exist_dist*/;
	int nNode_count, nAffect_m, nListF, nResult;
	
	DistCompute help;
	Point m, f;
	bool bMinDistUpdated, bMinExistDistUpdated;
	int nReturnResult = 1;
	for(int i = 0; i < GETF; i++)//Can modify GETF to change the number of F to expand, we use 1.
	{
		bool bResult = true;
		if(L_f.empty())
			break;

		do 
		{
			if(L_f.empty())
			{
				break;
			}
			//get imp(f)
			bResult = getImp_f(L_f, imp_f2);
			if(!bResult)
			{
				return 0;
			}
		} while (imp_f2->L_m.empty());

		if(imp_f2->pF->m_child)
		{
			//remove imp(f) from L_f
			it_m2 = imp_f2->L_m.begin();
			nAffect_m = imp_f2->L_m.size();
			for(int i = 0; i < nAffect_m; i++)
			{
				it_f2 = (*it_m2)->L_f.begin();
				nListF = (*it_m2)->L_f.size();
				for(int j = 0; j < nListF; j++)
				{
					if((*it_f2) == imp_f2)
					{
						(*it_m2)->L_f.erase(it_f2);
						break;
					}
					else
						it_f2++;
				}
				it_m2++;
			}//end remove imp(f) from m's L_f

			nNode_count = imp_f2->pF->m_child->m_count;
			for(int i = 0; i < nNode_count; i++)
			{
				//inherit relation information from parent node, and compute their own relation information
				pTempEnF->pF = &(imp_f2->pF->m_child->m_branch[i]);
				pTempEnF->square = ((imp_f2->pF->m_child->m_branch[i].m_rect.m_max[0] - imp_f2->pF->m_child->m_branch[i].m_rect.m_min[0]) 
									+ (imp_f2->pF->m_child->m_branch[i].m_rect.m_max[1] - imp_f2->pF->m_child->m_branch[i].m_rect.m_min[1]));
				pTempEnF->L_m = imp_f2->L_m;
				if(imp_f2->pF->m_child->IsLeaf())
				{
					pTempEnF->pF->m_child = NULL;
				}	
				tempL_f.push_back(pTempEnF);
				pTempEnF = new EntryF;
			}

			//update their own relation information
			it_m2 = imp_f2->L_m.begin();
			nAffect_m = imp_f2->L_m.size();
			for(int j = 0; j < nAffect_m; j++)
			{
				bMinDistUpdated = false, bMinExistDistUpdated = false;
				float fminMinDist, tempMinMinDist;
				tempIt_f_for_tempLf = tempL_f.begin();

				fminMinDist = help.minDist((*tempIt_f_for_tempLf)->pF->m_rect, (*it_m2)->pM->m_rect);
				minExist_dist = help.minExistDist((*it_m2)->pM->m_rect, (*tempIt_f_for_tempLf)->pF->m_rect);
				if(minExist_dist < (*it_m2)->minExist_dist)
				{
					(*it_m2)->minExist_dist = minExist_dist;
					bMinExistDistUpdated = true;
				}
				for(int i = 1; i < nNode_count; i++)
				{
					tempIt_f_for_tempLf++;
					tempMinMinDist = help.minDist((*tempIt_f_for_tempLf)->pF->m_rect, (*it_m2)->pM->m_rect);
					if(tempMinMinDist < fminMinDist)
						fminMinDist = tempMinMinDist;
					//update minExist_dist
					minExist_dist = help.minExistDist((*it_m2)->pM->m_rect, (*tempIt_f_for_tempLf)->pF->m_rect);
					if(minExist_dist < (*it_m2)->minExist_dist)
					{
						(*it_m2)->minExist_dist = minExist_dist;
						bMinExistDistUpdated = true;
					}
				}

				mindist_imp_f_M = help.minDist(imp_f2->pF->m_rect, (*it_m2)->pM->m_rect);
				if((*it_m2)->min_dist == mindist_imp_f_M)
				{
					if((*it_m2)->min_dist < fminMinDist)
					{
						int nAffect_f;
						float fminMinDist_it_m, tempMinMinDist_it_m;
						tempIt_f2 = (*it_m2)->L_f.begin();
						nAffect_f = (*it_m2)->L_f.size();
						fminMinDist_it_m = fminMinDist;
						for(int k = 0; k < nAffect_f; k++)
						{
							tempMinMinDist_it_m = help.minDist((*it_m2)->pM->m_rect, (*tempIt_f2)->pF->m_rect);
							if(tempMinMinDist_it_m < fminMinDist_it_m)
								fminMinDist_it_m = tempMinMinDist_it_m;
							tempIt_f2++;
						}
						(*it_m2)->min_dist = fminMinDist_it_m;
						bMinDistUpdated = true;
					}
				}
					
				if(bMinDistUpdated || bMinExistDistUpdated)
				{
					int nAffect_c;
					nAffect_c = (*it_m2)->L_c.size();
					tempIt_c2 = (*it_m2)->L_c.begin();
					for(int i = 0; i < nAffect_c; i++)
					{
						it_c2 = tempIt_c2;
						it_c2++;
						nResult = UpdateC_M2((*tempIt_c2), (*it_m2));//not updated return 0, update minInf return 1, update maxInf return 2
						if(nResult == 2)
							nReturnResult = 2;
						tempIt_c2 = it_c2;
						
					}
				}
				if(bMinExistDistUpdated)
				{
					int nAffect_f;
					tempIt_f2 = (*it_m2)->L_f.begin();
					nAffect_f = (*it_m2)->L_f.size();
					for(int k = 0; k < nAffect_f; k++)
					{
						tempIt_f_remove2 = tempIt_f2;

						temp_mindist = help.minDist((*it_m2)->pM->m_rect, (*tempIt_f2)->pF->m_rect);
						if(temp_mindist >= (*it_m2)->minExist_dist)
						{
							int nRemove_m;
							tempIt_m2 = (*tempIt_f2)->L_m.begin();
							nRemove_m = (*tempIt_f2)->L_m.size();
							for(int j = 0; j < nRemove_m; j++)
							{
								if((*tempIt_m2) == (*it_m2))
								{
									(*tempIt_f2)->L_m.erase(tempIt_m2);
									break;
								}
								else
								{
									tempIt_m2++;
								}
							}
							tempIt_f2++;
							(*it_m2)->L_f.erase(tempIt_f_remove2);
						}
						else
							tempIt_f2++;
					}
				}
				tempIt_f2 = tempL_f.begin();
				for(int i = 0; i < nNode_count; i++)
				{
					mindist_F_M = help.minDist((*tempIt_f2)->pF->m_rect, (*it_m2)->pM->m_rect);
					if(mindist_F_M <= (*it_m2)->minExist_dist)
					{
						(*it_m2)->L_f.push_back(*tempIt_f2);
					}
					else// if(mindist_F_M > (*it_m)->minExist_dist)
					{
						int nRemove_m;
						tempIt_m2 = (*tempIt_f2)->L_m.begin();
						nRemove_m = (*tempIt_f2)->L_m.size();
						for(int j = 0; j < nRemove_m; j++)
						{
							if((*tempIt_m2) == (*it_m2))
							{
								(*tempIt_f2)->L_m.erase(tempIt_m2);
								break;
							}
							else
							{
								tempIt_m2++;
							}
						}
					}
					tempIt_f2++;
				}
				it_m2++;
			}//end for each m in imp_f's L_m
			for(int i = 0; i < nNode_count; i++)
			{
				pTempEnF = tempL_f.back();
				tempL_f.pop_back();
				if(pTempEnF->L_m.size())
				{
					L_f.push_back(pTempEnF);
				}
			}
		}//end imp_f->m_child->m_count > 1
		else
		{
			L_f.push_back(pTempEnF);
		}
	}//for imp_f
	return nReturnResult;
}

//input list L_m, expand its imp(m)
bool Expand_M(list<EntryM*> &L_m)
{
	//expand imp(m)
	bool bIsMaxInfUpdate = false;
	int nNode_count = 0, nAffect_c = 0, nAffect_m = 0, nListF = 0;

	EntryM *imp_m2;
	EntryM *pTempM;
	imp_m2 = new EntryM;
	pTempM = new EntryM;
	list<EntryM*> tempQ_m;
	list<EntryC*>::iterator tempIt_c2, it_c2;
	list<EntryF*>::iterator tempIt_f2, it_f2;
	list<EntryM*>::iterator tempIt_m2, it_m2;
	list<EntryC*>::iterator next_it2;

	DistCompute help;
	float mindist_F_M = 0, maxdist_F_M = 0, minMaxDist = 0;
	float temp_maxdist = 0, temp_mindist = 0;
	float mindist_C_M = 0, maxdist_C_M = 0, fTempMaxDist = 0;
	Point m, f;
	float minExist_dist = 0, exist_dist = 0;
	for(int i = 0; i < GETM; i++) //can change GETM to choose how many entries to expand, we use 1.
	{
		bool bGetImp_m = false;
		if(L_m.empty())
			break;
		bGetImp_m = getImp_m(L_m, imp_m2);
		if(!bGetImp_m)
		{
			break;
		}
		if(imp_m2->amount_of_nodes > 1)
		{
			nNode_count = imp_m2->pM->m_child->m_count;
			//remove imp(m) from c's S_m
			nAffect_c = imp_m2->L_c.size();
			tempIt_c2 = imp_m2->L_c.begin();
			for(int i = 0; i < nAffect_c; i++)
			{
				tempIt_m2 = (*tempIt_c2)->L_m.begin();
				nAffect_m = (*tempIt_c2)->L_m.size();
				for(int j = 0; j < nAffect_m; j++)
				{
					if((*tempIt_m2) == imp_m2)
					{
						//remove tempIt_m
						(*tempIt_c2)->L_m.erase(tempIt_m2);
						break;
					}
					else
						tempIt_m2++;
				}//end-for each L_m
				tempIt_c2++;
			}//end-for each c

			//remove imp(m) from f's S_m
			nListF = imp_m2->L_f.size();
			it_f2 = imp_m2->L_f.begin();
			for(int i = 0; i < nListF; i++)
			{
				tempIt_m2 = (*it_f2)->L_m.begin();
				nAffect_m = (*it_f2)->L_m.size();
				for(int j = 0; j < nAffect_m; j++)
				{
					if((*tempIt_m2) == imp_m2)
					{
						//remove tempIt_m
						(*it_f2)->L_m.erase(tempIt_m2);
						break;
					}
					else
						tempIt_m2++;
				}//end-for each S_m
				it_f2++;
			}//end-for each f

			for(int i = 0; i < nNode_count; i++)
			{
				pTempM->pM = &(imp_m2->pM->m_child->m_branch[i]);
				pTempM->square = ((imp_m2->pM->m_child->m_branch[i].m_rect.m_max[0] - imp_m2->pM->m_child->m_branch[i].m_rect.m_min[0])
								+ (imp_m2->pM->m_child->m_branch[i].m_rect.m_max[1] - imp_m2->pM->m_child->m_branch[i].m_rect.m_min[1]));
				pTempM->L_c = imp_m2->L_c;
				pTempM->L_f = imp_m2->L_f;
				pTempM->amount_of_nodes = imp_m2->pM->m_child->m_branch[i].amount_of_nodes;
				tempQ_m.push_back(pTempM);
				pTempM = new EntryM;
			}
			tempIt_m2 = tempQ_m.begin();
			for(int i = 0; i < nNode_count; i++)
			{
 				//update child nodes' relation information by computing min_dist and minExist_dist for each m
				it_f2 = (*tempIt_m2)->L_f.begin();
				nListF = (*tempIt_m2)->L_f.size();
				mindist_F_M = help.minDist((*it_f2)->pF->m_rect, (*tempIt_m2)->pM->m_rect);
				minExist_dist = help.minExistDist((*tempIt_m2)->pM->m_rect, (*it_f2)->pF->m_rect);
				it_f2++;
				for(int j = 1; j < nListF; j++)
				{
					temp_mindist = help.minDist((*it_f2)->pF->m_rect, (*tempIt_m2)->pM->m_rect);
					exist_dist = help.minExistDist((*tempIt_m2)->pM->m_rect, (*it_f2)->pF->m_rect);
					if(exist_dist < minExist_dist)
					{
						minExist_dist = exist_dist;
					}
					if(temp_mindist < mindist_F_M)
					{
						mindist_F_M = temp_mindist;
					}
					it_f2++;		
				}
				(*tempIt_m2)->min_dist = mindist_F_M;
				(*tempIt_m2)->minExist_dist = minExist_dist;
				
				it_f2 = (*tempIt_m2)->L_f.begin();
				for(int j = 0; j < nListF; j++)
				{
					temp_mindist = help.minDist((*it_f2)->pF->m_rect, (*tempIt_m2)->pM->m_rect);
					if(temp_mindist > (*tempIt_m2)->minExist_dist)
					{
						tempIt_f2 = it_f2;
						tempIt_f2++;
						(*tempIt_m2)->L_f.erase(it_f2);
						it_f2 = tempIt_f2;
					}
					else
						it_f2++;
				}
				
				nAffect_c = (*tempIt_m2)->L_c.size();
				tempIt_c2 = (*tempIt_m2)->L_c.begin();
				
				next_it2 = tempIt_c2;
				int nIsUpdate = 0;
				for(int j = 0; j < nAffect_c; j++)
				{
					next_it2++;
					nIsUpdate = UpdateC_M(*tempIt_c2, *tempIt_m2);//return 2 when minInf is updated
					
					if(nIsUpdate == 2)
					{
						bIsMaxInfUpdate = true;
					}
					if(nIsUpdate == 0)
					{
						(*tempIt_c2)->L_m.push_back(*tempIt_m2);
						tempIt_c2++;
					}
					else
					{
						tempIt_c2 = next_it2;
					}
				}
				tempIt_m2++;
			}
			//add sub entries of imp_m to f's S_m
			tempIt_m2 = tempQ_m.begin();
			for(int i = 0; i < nNode_count; i++)
			{
				if(!(*tempIt_m2)->L_c.empty())
				{
					nListF = (*tempIt_m2)->L_f.size();
					it_f2 = (*tempIt_m2)->L_f.begin();
					for(int j = 0; j < nListF; j++)
					{
						(*it_f2)->L_m.push_back(*tempIt_m2);
						it_f2++;
					}
					L_m.push_back((*tempIt_m2));
				}
				tempIt_m2++;
			}
		}//end-if imp_m->amount_of_nodes > 1
		else
		{
			L_m.insert(L_m.end(), imp_m2);
		}//end else, imp_m->amount_of_nodes == 1
	}//for imp_m
	return bIsMaxInfUpdate;
}

//EEP algoirthm
bool EEP(aNode* root_m, Node* root_f, Node* root_c, int k, long nMa, long nFa, long nCa)
{
	int nSize, nSize_of_top = 0;
	int nNode_count = 0, count_c = 1, nAffect_c = 0, nAffect_m = 0, nListF = 0;
	int sort_count = 0, nSize_of_lm = 0, nSize_of_lf = 0;
	int nOutput = 1;
	int nMaxInf = 0, nTempMaxInf = 0;
	double imp_value_f = 0, imp_value_m = 0;
	DistCompute help;
	Sort sort_help;

	list<EntryC*> Q_c2;
	list<EntryC*> tempQ_c;
	list<EntryF*> Q_f;
	list<EntryM*> Q_m;
	list<EntryC*> output_c2, temp_result_c2;

	EntryC *pEnC, *pImpEnC, *secondC2, *firstC2;
	EntryF *pEnF;
	EntryM *pEnM;
	pEnC = new EntryC;
	pImpEnC = new EntryC;
	secondC2 = new EntryC;
	firstC2 = new EntryC;
	pEnF = new EntryF;
	pEnM = new EntryM;

	list<EntryC*>::iterator it_c2, tempIt_c2;
	list<EntryF*>::iterator it_f2, tempIt_f2;
	list<EntryM*>::iterator it_m2, tempIt_m2;
	
	float mindist_C_M = 0, maxdist_C_M = 0, mindist_F_M = 0, maxdist_F_M = 0, minMaxDist = 0, temp_mindist = 0, temp_maxdist = 0;
	float minExist_dist_for_MF = 0;
	float minExist_dist = 0, temp_dist = 0, fTempMaxDist = 0;
	Point m;
	int notSort = 0;
	//push root nodes into three lists, denoted by Q_m, Q_f, Q_c
	for(int i = 0; i < root_m->m_count; i++)
	{
		pEnM->pM = &(root_m->m_branch[i]);
		pEnM->amount_of_nodes = root_m->m_branch[i].amount_of_nodes;
		pEnM->square = ((root_m->m_branch[i].m_rect.m_max[0] - root_m->m_branch[i].m_rect.m_min[0])
			+ (root_m->m_branch[i].m_rect.m_max[1] - root_m->m_branch[i].m_rect.m_min[1]));
		Q_m.push_back(pEnM);
		pEnM = new EntryM;
	}
	for(int i = 0; i < root_f->m_count; i++)
	{
		pEnF->pF = &(root_f->m_branch[i]);
		pEnF->square = ((root_f->m_branch[i].m_rect.m_max[0] - root_f->m_branch[i].m_rect.m_min[0]) 
						+ (root_f->m_branch[i].m_rect.m_max[1] - root_f->m_branch[i].m_rect.m_min[1]));
		Q_f.push_back(pEnF);
		pEnF = new EntryF;
	}
	
	for(int i = 0; i < root_c->m_count; i++)
	{
		pEnC->pC = &(root_c->m_branch[i]);
		if(root_c->IsLeaf())
		{
			pEnC->pC->m_child = NULL;
		}
		Q_c2.push_back(pEnC);
		pEnC = new EntryC;
	}
	//expand Q_c to at least k entries
	nSize = Q_c2.size();
	while(nSize < k)
	{
		Expand(Q_c2, k);
		nSize = Q_c2.size();
	}
	nSize = Q_c2.size();
	count_c = nSize + 1;
	//update minInf and maxInf for each c
	UpdateC(Q_c2, Q_f, Q_m, nMa);
	long start_sort, end_sort, time_for_sort = 0;
	bool bNeedSort = false, bNewEntriesAdded = false;
	bool bExpandSort = false;
	int nExpandF = 1;
	list<EntryC*>::iterator it_imp_c2, it_second_c2;

	while(1)
	{
		
		bExpandSort = Expand_M(Q_m);
		if(nExpandF)
		{
			nExpandF = Expand_F(Q_f);
		}

		if(bExpandSort || nExpandF == 2 || bNeedSort)//in some case, Q_c is already soted, so we no need to sort it again
		{
			//sort list Q_c, so as to get imp(c)
			sort_help.InsertionSort(Q_c2);
			bNeedSort = false;
			bExpandSort = false;
			it_imp_c2 = Q_c2.begin();
			pImpEnC = (*it_imp_c2);
		}
		else
		{
			notSort++;
			it_imp_c2 = Q_c2.begin();
			pImpEnC = (*it_imp_c2);
		}
		it_second_c2 = it_imp_c2;
		it_second_c2++;
		secondC2 = (*it_second_c2);
		if(pImpEnC->pC->m_child != NULL)
		{
			Q_c2.pop_front();
			nMaxInf = secondC2->maxInf;
			if(pImpEnC->minInf >= nMaxInf)
			{
				nNode_count = pImpEnC->pC->m_child->m_count;
				for(int i = 0; i < nNode_count; i++)
				{
					pEnC->pC = &(pImpEnC->pC->m_child->m_branch[i]);
					pEnC->maxInf = pImpEnC->maxInf;
					pEnC->minInf = pImpEnC->minInf;
					result_c2.push_back(pEnC);
					pEnC = new EntryC;
				}
				continue;
			}
			//Expand imp_c
			nNode_count = pImpEnC->pC->m_child->m_count;
			for(int i = 0; i < nNode_count; i++)
			{
				//child nodes of list C inherit relation information from parent node
				pEnC->pC = &(pImpEnC->pC->m_child->m_branch[i]);
				if(pImpEnC->pC->m_child->m_level == 0)
				{
					pEnC->pC->m_child = NULL;
				}
				pEnC->maxInf = pImpEnC->maxInf;
				pEnC->minInf = pImpEnC->minInf;
				pEnC->L_m = pImpEnC->L_m;
				tempQ_c.push_front(pEnC);
				pEnC = new EntryC;
			}

			nAffect_m = pImpEnC->L_m.size();
			it_m2 = pImpEnC->L_m.begin();
			for(int i = 0; i < nAffect_m; i++)
			{
				//remove imp_c from m's S_c set
 				tempIt_c2 = (*it_m2)->L_c.begin();
				nSize = (*it_m2)->L_c.size();
				for (int j = 0; j < nSize; j++)
				{			
					if(pImpEnC == (*tempIt_c2))
					{
						(*it_m2)->L_c.erase(tempIt_c2);
						break;
					}
					else
						tempIt_c2++;
				}//end-for remove
				//add new entries to L_m's L_c
				tempIt_c2 = tempQ_c.begin();
				for(int j = 0; j < nNode_count; j++)
				{
					mindist_C_M = help.minDist((*tempIt_c2)->pC->m_rect, (*it_m2)->pM->m_rect);
					maxdist_C_M = help.maxDist((*tempIt_c2)->pC->m_rect, (*it_m2)->pM->m_rect);
					//m is only affected by C
					if((*it_m2)->min_dist >= maxdist_C_M) // m's min_dist and max_dist are distances from L_f to m
					{
						//minus m to c's minInf
						(*tempIt_c2)->minInf += (*it_m2)->amount_of_nodes;
						tempIt_m2 = (*tempIt_c2)->L_m.begin();
						nSize_of_lm = (*tempIt_c2)->L_m.size();
						for(int j = 0; j < nSize_of_lm; j++)
						{
							if((*tempIt_m2) == (*it_m2))
							{
								(*tempIt_c2)->L_m.erase(tempIt_m2);
								break;
							}
							else
							{
								tempIt_m2++;
							}
						}
					}
					//m is only affected by L_f
					else if((*it_m2)->minExist_dist < mindist_C_M)
					{
						//remove m from Lm and minus c's maxInf by m
						//minus m to c's maxInf
						(*tempIt_c2)->maxInf -= (*it_m2)->amount_of_nodes;
						if(secondC2->maxInf > (*tempIt_c2)->maxInf)
							bNeedSort = true;
						tempIt_m2 = (*tempIt_c2)->L_m.begin();
						nSize_of_lm = (*tempIt_c2)->L_m.size();
						for(int j = 0; j < nSize_of_lm; j++)
						{
							if((*tempIt_m2) == (*it_m2))
							{
								(*tempIt_c2)->L_m.erase(tempIt_m2);
								break;
							}
							else
							{
								tempIt_m2++;
							}
						}
					}
					else
					{
						//notice L_m's entries, new Cs are added
						(*it_m2)->L_c.push_back((*tempIt_c2));
					}
					tempIt_c2++;
				}
				it_m2++;
			}//end-for each m
			tempIt_c2 = tempQ_c.begin();
			for(int i = 0; i < nNode_count; i++)
			{
				Q_c2.push_front(*tempIt_c2);
				tempIt_c2++;
			}
		}//end if imp_c->amount_of_nodes > 1
		else
		{
			nMaxInf = secondC2->maxInf;
			if(pImpEnC->minInf >= nMaxInf)
			{
				do
				{
					result_c2.push_back(pImpEnC);
					Q_c2.pop_front();
					it_imp_c2 = it_second_c2;
					pImpEnC = (*it_imp_c2);
					it_second_c2++;
					if(it_second_c2 == Q_c2.end())
					{
						result_c2.push_back(pImpEnC);
						if(result_c2.size() == K)
							break;
					}
					secondC2 = (*it_second_c2);
					nMaxInf = secondC2->maxInf;
				}while((result_c2.size() < K && pImpEnC->minInf >= nMaxInf));
			}
		}

		if(result_c2.size() == K)
		{
			return true;
		}
	}
	return true;
}
#endif INFLUENCE_H