#ifndef DISTCOMPUTE_H
#define DISTCOMPUTE_H

#include "Utility.h"

#define HEAP_MAX_SIZE 400
long long counter = 0;
long long ExactCounter = 0;

class DistCompute
{
public:
	DistCompute(){}
	~DistCompute(){}
	float distance (float x1, float y1, float x2, float y2)
	{
		float tempdist;
		tempdist = ((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
		//distance counter
		counter++;
		return tempdist;
	}

	float distance(Point A, Point B)
	{
		float tempdist;
		tempdist = ((A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y));
		//distance counter
		counter++;
		tempdist = sqrt(tempdist);
		return tempdist;
	}
public:
	float minDist(Point p, Rect r)
	{
		float minDist_p_r = 0, tempDist = 0;
		Point A, B, C, D;
		Point temp;
		DistCompute help;

		A.x = r.m_min[0];
		A.y = r.m_min[1];
		B.x = r.m_max[0];
		B.y = r.m_max[1];

		temp.x = r.m_min[0];
		temp.y = r.m_max[1];
		C.x = temp.x;
		C.y = temp.y;
		temp.x = r.m_max[0];
		temp.y = r.m_min[1];
		D.x = temp.x;
		D.y = temp.y;

		//////////////////////////////////////////////////////////////////////////
		//exact distance counter
		if(A.x == B.x && A.y == B.y)
			ExactCounter++;
		//////////////////////////////////////////////////////////////////////////

		//determine which side of the rectangle p is located.
		if(p.x < A.x)
		{
			if(p.y <= B.y && p.y >= A.y)
			{
				//compute distance from p to AC line segment
				minDist_p_r = A.x - p.x;
			}
			else
			{
				tempDist = help.distance(p, C);
				minDist_p_r = help.distance(p, A);
				minDist_p_r = (tempDist < minDist_p_r ? tempDist : minDist_p_r);
			}
		}
		else if(p.x > B.x)
		{
			if(p.y <= B.y && p.y >= A.y)
			{
				//compute distance from p to BD line segment
				minDist_p_r = p.x - B.x;
			}
			else
			{
				tempDist = help.distance(p, B);
				minDist_p_r = help.distance(p, D);
				minDist_p_r = (tempDist < minDist_p_r ? tempDist : minDist_p_r);
			}
		}
		else if(p.y < A.y)
		{
			//compute distance from p to AD
			minDist_p_r = A.y - p.y;
		}
		else if(p.y > B.y)
		{
			//compute distance from p to CB
			minDist_p_r = p.y - B.y;
		}
		else
		{
			minDist_p_r = 0;
		}
		return minDist_p_r;
	}
//Minimum distance between two rectangles.
	float minDist(Rect R_A, Rect R_B)
	{
		float mindist;

		if(R_A.m_min[0] <= R_B.m_max[0] && R_B.m_min[0] <= R_A.m_max[0] && R_A.m_min[1] <= R_B.m_max[1] && R_B.m_min[1] <= R_A.m_max[1])
		{
			return 0;
		}
		Point R_Aa, R_Ab, R_Ac, R_Ad;  // c------------------b
		R_Aa.x = R_A.m_min[0];         // |                  |
		R_Aa.y = R_A.m_min[1];         // |                  |
									   // |                  |
		R_Ab.x = R_A.m_max[0];         // |                  |
		R_Ab.y = R_A.m_max[1];	       // a------------------d

		R_Ac.x = R_A.m_min[0];
		R_Ac.y = R_A.m_max[1];

		R_Ad.x = R_A.m_max[0];
		R_Ad.y = R_A.m_min[1];

		Point R_Ba, R_Bb, R_Bc, R_Bd;
		R_Ba.x = R_B.m_min[0];
		R_Ba.y = R_B.m_min[1];

		R_Bb.x = R_B.m_max[0];
		R_Bb.y = R_B.m_max[1];

		R_Bc.x = R_B.m_min[0];
		R_Bc.y = R_B.m_max[1];

		R_Bd.x = R_B.m_max[0];
		R_Bd.y = R_B.m_min[1];

		//////////////////////////////////////////////////////////////////////////
		//exact distance counter
		if(R_Aa.x == R_Ab.x && R_Aa.y == R_Ab.y && R_Ba.x == R_Bb.x && R_Ba.y == R_Bb.y)
			ExactCounter++;
		//////////////////////////////////////////////////////////////////////////

		if(R_A.m_max[0] < R_B.m_min[0])//R_A is to the left side of R_B
		{
			if((R_A.m_min[1] < R_B.m_max[1] && R_A.m_max[1] > R_B.m_max[1]) || (R_A.m_max[1] > R_B.m_min[1] && R_A.m_min[1] < R_B.m_min[1])
				|| (R_B.m_min[1] < R_A.m_min[1] && R_A.m_min[1] < R_B.m_max[1]) || (R_B.m_min[1] < R_A.m_max[1] && R_A.m_max[1] < R_B.m_max[1]))//R_A and R_B have a vertical line
			{
				//X axis distance
				mindist = R_B.m_min[0] - R_A.m_max[0];
			}
			else //the minimum value of 4 points distance
			{
				if(R_Ad.y > R_Bc.y)
					mindist = distance(R_Ad, R_Bc);
				else
					mindist = distance(R_Ab, R_Ba);
			}
			return mindist;
		}
		else if(R_A.m_min[0] > R_B.m_max[0])//R_A is to the right side of R_B
		{
			if((R_A.m_min[1] < R_B.m_max[1] && R_A.m_max[1] > R_B.m_max[1]) || (R_A.m_max[1] > R_B.m_min[1] && R_A.m_min[1] < R_B.m_min[1])
				|| (R_B.m_min[1] < R_A.m_min[1]) && (R_A.m_min[1] < R_B.m_max[1]) || (R_B.m_min[1] < R_A.m_max[1]) && (R_A.m_max[1] < R_B.m_max[1]))//R_A and R_B have a vertical line
			{
				//X axis distance
				mindist = R_A.m_min[0] - R_B.m_max[0];
			}
			else //the minimum value of 4 points distance
			{
				if(R_Aa.y > R_Bb.y)
					mindist = distance(R_Aa, R_Bb);
				else
					mindist = distance(R_Ac, R_Bd);
			}
			return mindist;
		}
		else if(R_A.m_min[1] > R_B.m_max[1])//R_A is on the top of R_B
		{
			mindist = R_Aa.y - R_Bc.y;
			return mindist;
		}
		else if(R_A.m_max[1] < R_B.m_min[1])//R_A is on the bottom of R_B
		{
			mindist = R_Ba.y - R_Ac.y;
			return mindist;
		}
		else
			return -1;
	}

//maximum distance from point p to rectangle R
	float maxDist(Point p, Rect R)
	{
		float fTempDist1 = 0, fTempDist2 = 0, fTempDist3 = 0,fMaxDist;
		Point Ra, Rb, Rc, Rd;      // c------------------b
		Ra.x = R.m_min[0];         // |                  |
		Ra.y = R.m_min[1];         // |                  |
								   // |                  |
		Rb.x = R.m_max[0];         // |                  |
		Rb.y = R.m_max[1];	       // a------------------d

		Rc.x = R.m_min[0];
		Rc.y = R.m_max[1];

		Rd.x = R.m_max[0];
		Rd.y = R.m_min[1];
		//////////////////////////////////////////////////////////////////////////
		//exact distance counter
		if(Ra.x == Rb.x && Ra.y == Rb.y)
			ExactCounter++;
		//////////////////////////////////////////////////////////////////////////
		fTempDist1 = distance(p, Ra);
		fTempDist2 = distance(p, Rb);
		fTempDist3 = distance(p, Rc);
		fMaxDist = distance(p, Rd);
		assert(fMaxDist >= 0 && fTempDist1 >= 0 && fTempDist2 >= 0 && fTempDist1 >= 0);
		if(fMaxDist < fTempDist1)
			fMaxDist = fTempDist1;
		if(fMaxDist < fTempDist2)
			fMaxDist = fTempDist2;
		if(fMaxDist < fTempDist3)
			fMaxDist = fTempDist3;
		return fMaxDist;
	}

	//maximum distance between two rectangles
	float maxDist(Rect R_A, Rect R_B)
	{
		float fTempDist1 = 0, fTempDist2 = 0, fTempDist3 = 0,fMaxDist;
		
		Point R_Aa, R_Ab, R_Ac, R_Ad;  // c------------------b
		R_Aa.x = R_A.m_min[0];         // |                  |
		R_Aa.y = R_A.m_min[1];         // |                  |
									   // |                  |
		R_Ab.x = R_A.m_max[0];         // |                  |
		R_Ab.y = R_A.m_max[1];	       // a------------------d

		R_Ac.x = R_A.m_min[0];
		R_Ac.y = R_A.m_max[1];

		R_Ad.x = R_A.m_max[0];
		R_Ad.y = R_A.m_min[1];
		
		Point R_Ba, R_Bb, R_Bc, R_Bd;
		R_Ba.x = R_B.m_min[0];
		R_Ba.y = R_B.m_min[1];

		R_Bb.x = R_B.m_max[0];
		R_Bb.y = R_B.m_max[1];

		R_Bc.x = R_B.m_min[0];
		R_Bc.y = R_B.m_max[1];

		R_Bd.x = R_B.m_max[0];
		R_Bd.y = R_B.m_min[1];
		//////////////////////////////////////////////////////////////////////////
		//exact distance counter
		if(R_Aa.x == R_Ab.x && R_Aa.y == R_Ab.y && R_Ba.x == R_Bb.x && R_Ba.y == R_Bb.y)
			ExactCounter++;
		//////////////////////////////////////////////////////////////////////////
		fTempDist1 = distance(R_Aa, R_Bb);
		fTempDist2 = distance(R_Ab, R_Ba);
		fTempDist3 = distance(R_Ac, R_Bd);
		fMaxDist = distance(R_Ad, R_Bc);
		assert(fMaxDist > 0 && fTempDist1 > 0 && fTempDist2 > 0 && fTempDist1 > 0);
		if(fMaxDist < fTempDist1)
			fMaxDist = fTempDist1;
		if(fMaxDist < fTempDist2)
			fMaxDist = fTempDist2;
		if(fMaxDist < fTempDist3)
			fMaxDist = fTempDist3;
		return fMaxDist;
	}
	//compute minMaxDist
	float minMaxDist(Point p, Rect R)
	{
		float fTempDist1 = 0, fTempDist2 = 0, fTempDist3 = 0,fminMaxDist;
		Point Ra, Rb, Rc, Rd;      // c------------------b
		Ra.x = R.m_min[0];         // |                  |
		Ra.y = R.m_min[1];         // |                  |
								   // |                  |
		Rb.x = R.m_max[0];         // |                  |
		Rb.y = R.m_max[1];	       // a------------------d
		
		//////////////////////////////////////////////////////////////////////////
		//exact distance counter
		if(Ra.x == Rb.x && Ra.y == Rb.y)
			ExactCounter++;
		//////////////////////////////////////////////////////////////////////////
		//p in the left of R
		if(p.x < Ra.x)
		{
			fTempDist1 = distance(p, Ra);
			fTempDist2 = distance(p, Rc);
		}
		//p in the right of R
		else if(p.x > Rb.x)
		{
			fTempDist1 = distance(p, Rb);
			fTempDist2 = distance(p, Rd);
		}
		//p in the top of R
		else if (p.y > Rb.y)
		{
			fTempDist1 = distance(p, Rc);
			fTempDist2 = distance(p, Rb);
		}
		//p in the bottom of R
		else if(p.y < Ra.y)
		{
			fTempDist1 = distance(p, Ra);
			fTempDist2 = distance(p, Rd);
		}
		//p inside R
		else
			return 0;

		if(fTempDist1 < fTempDist2)
		{
			fminMaxDist = fTempDist1;
		}
		else 
			fminMaxDist = fTempDist2;

		return fminMaxDist;
	}
	float minExistDist(Rect &Rm, Rect &Rf)
	{
		float fminExistDist;
								// c--p1-----------p3-b    c-------------------b 
								// |   \          /   |    |                   |
								// |     \ Rm   /     |    |      Rf           |
								// |       \  /       |    |                   |
								// |        /\        |    |                   |
								// a------p4---p2-----d    a-------------------d

		Point Rm_a, Rm_b, Rm_c, Rm_d;
		Point Rf_a, Rf_b, Rf_c, Rf_d;

		Rm_a.x = Rm.m_min[0];
		Rm_a.y = Rm.m_min[1];
		Rm_b.x = Rm.m_max[0];
		Rm_b.y = Rm.m_max[1];

		Rm_c.x = Rm.m_min[0];
		Rm_c.y = Rm.m_max[1];
		Rm_d.x = Rm.m_max[0];
		Rm_d.y = Rm.m_min[1];

		Rf_a.x = Rf.m_min[0];
		Rf_a.y = Rf.m_min[1];
		Rf_b.x = Rf.m_max[0];
		Rf_b.y = Rf.m_max[1];

		if(Rm_a.x == Rm_b.x && Rm_a.y == Rm_b.y)
		{
			fminExistDist = maxDist(Rm_a, Rf);
			return fminExistDist;
		}

		if(Rf_a.x == Rf_b.x && Rf_a.y == Rf_b.y)
		{
			fminExistDist = maxDist(Rf_a, Rm);
			return fminExistDist;
		}
		
		Rf_c.x = Rf.m_min[0];
		Rf_c.y = Rf.m_max[1];
		Rf_d.x = Rf.m_max[0];
		Rf_d.y = Rf.m_min[1];

		Point centre;
		Point p1, p2, p3, p4;
		bool bExist[4] = {false, false, false, false};
		float slope_f_ab, slope_f_cd;
		float side_for_ab_mb, side_for_ab_ma, side_for_ab_mc, side_for_ab_md;
		float side_for_cd_ma, side_for_cd_mb, side_for_cd_mc, side_for_cd_md;
		float dist1 = 0, dist2 = 0, dist3 = 0, dist4 = 0, dist_ma = 0, dist_mb = 0, dist_mc = 0, dist_md = 0;

		centre.x = (Rf_a.x + Rf_b.x) / 2;
		centre.y = (Rf_a.y + Rf_b.y) / 2;
		slope_f_ab = -(Rf_a.x - Rf_b.x) / (Rf_a.y - Rf_b.y);
		slope_f_cd = -(Rf_c.x - Rf_d.x) / (Rf_c.y - Rf_d.y);

		side_for_ab_mb = Rm_b.y - centre.y - slope_f_ab * (Rm_b.x - centre.x);
		side_for_ab_ma = Rm_a.y - centre.y - slope_f_ab * (Rm_a.x - centre.x);
		side_for_ab_md = Rm_d.y - centre.y - slope_f_ab * (Rm_d.x - centre.x);
		side_for_ab_mc = Rm_c.y - centre.y - slope_f_ab * (Rm_c.x - centre.x);

		side_for_cd_md = Rm_d.y - centre.y - slope_f_cd * (Rm_d.x - centre.x);
		side_for_cd_mc = Rm_c.y - centre.y - slope_f_cd * (Rm_c.x - centre.x);
		side_for_cd_ma = Rm_a.y - centre.y - slope_f_cd * (Rm_a.x - centre.x);
		side_for_cd_mb = Rm_b.y - centre.y - slope_f_cd * (Rm_b.x - centre.x);

		if(side_for_ab_mb > 0 && side_for_ab_ma < 0)// line cross Rm with two interaction points
		{
			//p1 may be on the left edge
			p1.y = Rm_b.y;
			p1.x = (Rm_b.y - centre.y) / slope_f_ab + centre.x;
			if(p1.x < Rm_a.x)
			{
				p1.x = Rm_a.x;
				p1.y = slope_f_ab * (Rm_a.x - centre.x) + centre.y;
			}
			//p2 may be on the right edge
			p2.y = Rm_a.y;
			p2.x = (Rm_a.y - centre.y) / slope_f_ab + centre.x;
			if(p2.x > Rm_b.x)
			{
				p2.x = Rm_b.x;
				p2.y = slope_f_ab * (Rm_b.x - centre.x) + centre.y;
			}
			bExist[0] = true;
			bExist[1] = true;
		}
		else if(side_for_ab_mb <= 0) // none interaction point, located in left
		{
			//
		}
		else if(side_for_ab_ma >= 0) // none interaction point, located in right
		{
			//
		}
		else
			cout << "error in minExistDist" << endl;

		if(side_for_cd_md < 0 && side_for_cd_mc > 0) //line interaction with Rm at lease two points
		{
			p3.y = Rm_d.y;
			p3.x = (Rm_d.y - centre.y) / slope_f_cd + centre.x;
			if(p3.x < Rm_c.x)
			{
				p3.x = Rm_c.x;
				p3.y = slope_f_cd * (Rm_c.x - centre.x) + centre.y;
			}
			p4.y = Rm_c.y;
			p4.x = (Rm_c.y - centre.y) / slope_f_cd + centre.x;
			if(p4.x > Rm_d.x)
			{
				p4.x = Rm_d.x;
				p4.y = slope_f_cd * (Rm_d.x - centre.x) + centre.y;
			}
			bExist[2] = true;
			bExist[3] = true;
		}
		else if(side_for_cd_md >= 0) //none interaction point
		{
			//
		}
		else if(side_for_cd_mc <= 0) //none interaction point
		{
			//
		}
		else
			cout << "error in minExistDist" << endl;
		//compute minExistDist
		if(bExist[0])
		{
			dist1 = distance(p1, Rf_b);
			dist2 = distance(p2, Rf_b);
			if(Rm_a.y > centre.y)//on the top of the region
			{
				dist_ma = distance(Rm_a, Rf_a);//a is always in the second 2nd sector
				if(side_for_ab_mc >= 0)//c is in second 1st sector
				{
					dist_mc = distance(Rm_c, Rf_b);
				}
				else//c is in second 2nd sector
				{
					dist_mc = distance(Rm_c, Rf_a);
				}

				if(Rm_b.x <= centre.x)//b in second 1st sector
				{
					dist_mb = distance(Rm_b, Rf_b);
					if(side_for_ab_md >= 0)//d in second 1st sector
					{
						dist_md = distance(Rm_d, Rf_b);
					}
					else //d in second 2nd sector
					{
						dist_md = distance(Rm_d, Rf_a);
					}
				}
				else if(side_for_cd_mb >= 0)//b in first 2nd sector
				{
					dist_mb = distance(Rm_b, Rf_c);
					if(side_for_cd_md >= 0)//d in first 2nd sector
					{
						dist_md = distance(Rm_d, Rf_c);
					}
					else //d in first 1st sector
					{
						dist_md = distance(Rm_d, Rf_d);
					}
				}
				else //located in first 1st sector
				{
					dist_mb = distance(Rm_b, Rf_d);
					dist_md = distance(Rm_d, Rf_d);
				}
			}
			else if(Rm_b.y < centre.y)//on the bottom of the region
			{
				//b is always on the fourth 2nd sector
				dist_mb = distance(Rm_b, Rf_b);
				if(side_for_cd_md >= 0)//d is on the fourth 2nd sector
				{
					dist_md = distance(Rm_d, Rf_b);
				}
				else//d in the fourth 1st sector
				{
					dist_md = distance(Rm_d, Rf_a);
				}
				if(Rm_a.x <= centre.x)
				{
					if(side_for_cd_ma >= 0)//a and c are in third 1st sector
					{
						dist_ma = distance(Rm_a, Rf_c);
						dist_mc = distance(Rm_c, Rf_c);
					}
					else //a is in third 2nd sector
					{
						dist_ma = distance(Rm_a, Rf_d);
						if(side_for_cd_mc >= 0)//c is in third 1st sector
						{
							dist_mc = distance(Rm_c, Rf_c);
						}
						else//c is in third 2nd sector
						{
							dist_mc = distance(Rm_c, Rf_d);
						}
					}
				}
				else//a and c are in fourth sector
				{
					dist_ma = distance(Rm_a, Rf_a);
					if(side_for_ab_mc >= 0)//c in fourth 2nd sector
					{
						dist_mc = distance(Rm_c, Rf_b);
					}
					else //c in fourth 1st sector
					{
						dist_mc = distance(Rm_c, Rf_a);
					}
				}
			}
			else //between centre
			{
				if(Rm_b.x <= centre.x)//Rm is in the left of the region
				{
					dist_mb = distance(Rm_b, Rf_b);
					if(side_for_ab_md <= 0 && side_for_cd_md > 0)//md is in the third 1st sector
					{
						dist_md = distance(Rm_d, Rf_c);
					}
					else if(side_for_cd_md <=0)//md is in the third 2nd sector
					{
						dist_md = distance(Rm_d, Rf_d);
					}
					else
						cout << "error in minExistDist" << endl;

					if(side_for_ab_mc >= 0)//c is in the second 1st sector
					{
						dist_mc = distance(Rm_c, Rf_b);
					}
					else//c is in the second 2nd sector
					{
						dist_mc = distance(Rm_c, Rf_a);
					}
					if(side_for_cd_ma >= 0)//a is in the third 1st sector
					{
						dist_ma = distance(Rm_a, Rf_c);
					}
					else//a is in the third 2nd sector
					{
						dist_ma = distance(Rm_a, Rf_d);
					}
				}
				else if(Rm_c.x < centre.x)//Rm is in the middle of the region
				{
					if(side_for_cd_mb >= 0)//b is in first 2nd sector
					{
						dist_mb = distance(Rm_b, Rf_c);
						if(side_for_ab_md >= 0)//d is in the fourth 2nd sector
						{
							dist_md = distance(Rm_d, Rf_b);
						}
						else//d is in the fourth 1st sector
						{
							dist_md = distance(Rm_d, Rf_a);
						}
					}
					else//b is in the first 1st sector
					{
						dist_mb = distance(Rm_b, Rf_d);
						if(side_for_ab_md >= 0)//d is in the fourth 2nd sector
						{
							dist_md = distance(Rm_d, Rf_b);
						}
						else//d is in the fourth 1st sector
						{
							dist_md = distance(Rm_d, Rf_a);
						}
					}

					if(side_for_ab_mc >= 0)//c is in the second 1st sector
					{
						dist_mc = distance(Rm_c, Rf_b);
					}
					else//c is in the second 2nd sector
					{
						dist_mc = distance(Rm_c, Rf_a);
					}
					if(side_for_cd_ma >= 0)//a is in the third 1st sector
					{
						dist_ma = distance(Rm_a, Rf_c);
					}
					else//a is in the third 2nd sector
					{
						dist_ma = distance(Rm_a, Rf_d);
					}
				}
				else//Rm is in the right of the region
				{
					if(side_for_cd_mb >=0)//b is in the first 2nd sector
					{
						dist_mb = distance(Rm_b, Rf_c);
						if(side_for_ab_md >= 0)//d is in the fourth 2nd sector
						{
							dist_md = distance(Rm_d, Rf_b);
						}
						else//d is in the fourth 1st sector
						{
							dist_md = distance(Rm_d, Rf_a);
						}
					}
					else//b is in the first 1st sector
					{
						dist_mb = distance(Rm_b, Rf_d);
						if(side_for_ab_md >= 0)//d is in fourth 2nd sector
						{
							dist_md = distance(Rm_d, Rf_b);
						}
						else//d is in the fourth 1st sector
						{
							dist_md = distance(Rm_d, Rf_a);
						}
					}
					//////////////////////////////////////////////////////////////////////////
					if(side_for_cd_mc >= 0)//c is in the first 2nd sector
					{
						dist_mc = distance(Rm_c, Rf_c);
					}
					else//c is in the first 1st sector
					{
						dist_mc = distance(Rm_c, Rf_d);
					}
					if(side_for_ab_ma >= 0)//a is in the fourth 2nd sector
					{
						dist_ma = distance(Rm_a, Rf_b);
					}
					else//a is in the fourth 1st sector
					{
						dist_ma = distance(Rm_a, Rf_a);
					}
				}
			}
			if(bExist[2])
			{
				dist3 = distance(p3, Rf_d);
				dist4 = distance(p4, Rf_d);
			}
		}//end if(bExist[0])
		else//p1 and p2 are not exist
		{
			if(side_for_ab_mb <= 0)//Rm_b is in the left of ab
			{
				if(Rm_d.y > centre.y)//top of region
				{
					dist_ma = distance(Rm_a, Rf_a);
					dist_mb = distance(Rm_b, Rf_a);
					dist_mc = distance(Rm_c, Rf_a);
					dist_md = distance(Rm_d, Rf_a);
				}
				else if(Rm_b.y < centre.y)//bottom of region
				{
					if(Rm_b.x < centre.x)
					{
						if(side_for_cd_md >= 0)//abcd are in third 1st sector
						{
							dist_ma = distance(Rm_a, Rf_c);
							dist_mb = distance(Rm_b, Rf_c);
							dist_mc = distance(Rm_c, Rf_c);
							dist_md = distance(Rm_d, Rf_c);
						}
						else//d is in the third 2nd sector
						{
							dist_md = distance(Rm_d, Rf_d);
							if(side_for_cd_mb >= 0)//b is in the third 1st sector
							{
								dist_mb = distance(Rm_b, Rf_c);
							}
							else
							{
								dist_mb = distance(Rm_b, Rf_d);
							}

							if(side_for_cd_ma >= 0)//a is in the third 1st sector
							{
								dist_ma = distance(Rm_a, Rf_c);
							}
							else//d is in the third 2nd sector
							{
								dist_ma = distance(Rm_a, Rf_d);
							}

							if(side_for_cd_mc >= 0)//c is in the third 1st sector
							{
								dist_mc = distance(Rm_c, Rf_c);
							}
							else//c is in the third 2nd sector
							{
								dist_mc = distance(Rm_c, Rf_d);
							}
						}
					}
					else//b,d are in the fourth 1st sector
					{
						dist_mb = distance(Rm_b, Rf_a);
						dist_md = distance(Rm_d, Rf_a);
						if(Rm_a.x >= centre.x)
						{
							dist_ma = distance(Rm_a, Rf_a);
							dist_mc = distance(Rm_c, Rf_a);
						}
						else
						{
							if(side_for_cd_ma >= 0)//a and c are in the third 1st sector
							{
								dist_ma = distance(Rm_a, Rf_c);
								dist_mc = distance(Rm_c, Rf_c);
							}
							else//a is in the third 2nd sector
							{
								dist_ma = distance(Rm_a, Rf_d);
								if(side_for_cd_mc >= 0)//c is in the third 1st sector
								{
									dist_mc = distance(Rm_c, Rf_c);
								}
								else
								{
									dist_mc = distance(Rm_c, Rf_d);
								}
							}
						}
					}
				}//end bottom
				else//between centre
				{
					if(Rm_b.y >= centre.y)
					{
						dist_mc = distance(Rm_c, Rf_a);
						dist_mb = distance(Rm_b, Rf_a);
					}
					else
						cout << "error in minExistDist" << endl;
					if(side_for_cd_md >= 0)
					{
						dist_ma = distance(Rm_a, Rf_c);
						dist_md = distance(Rm_d, Rf_c);
					}
					else//d is in third 2nd sector
					{
						dist_md = distance(Rm_d, Rf_d);
						if(side_for_cd_ma >= 0)
						{
							dist_ma = distance(Rm_a, Rf_c);
						}
						else
							dist_ma = distance(Rm_a, Rf_d);
					}
				}//end between centre
			}
			else if(side_for_ab_ma >= 0)//Rm_a is in the right of ab
			{
				if(Rm_b.y < centre.y)//on the bottom of region, abcd are in the fourth 2nd sector
				{
					dist_ma = distance(Rm_a, Rf_b);
					dist_mb = distance(Rm_b, Rf_b);
					dist_mc = distance(Rm_c, Rf_b);
					dist_md = distance(Rm_d, Rf_b);
				}
				else if(Rm_a.y > centre.y)//on the top region
				{
					if(Rm_b.x <= centre.x)//abcd are in the second 1st sector
					{
						dist_ma = distance(Rm_a, Rf_b);
						dist_mb = distance(Rm_b, Rf_b);
						dist_mc = distance(Rm_c, Rf_b);
						dist_md = distance(Rm_d, Rf_b);
					}
					else
					{
						if(side_for_cd_md >= 0)//b and d are in the first 2nd sector
						{
							dist_mb = distance(Rm_b, Rf_c);
							dist_md = distance(Rm_d, Rf_c);
						}
						else//d is in the first 1st sector
						{
							dist_md = distance(Rm_d, Rf_d);
							if(side_for_cd_mb >= 0)//b is in the first 2nd sector
							{
								dist_mb = distance(Rm_b, Rf_c);
							}
							else//b is in the first 1st sector
							{
								dist_mb = distance(Rm_b, Rf_d);
							}
						}

						if(Rm_c.x <= centre.x)//a and c are in the second 1st sector
						{
							dist_ma = distance(Rm_a, Rm_b);
							dist_mc = distance(Rm_c ,Rm_b);
						}
						else
						{
							if(side_for_cd_ma >= 0)//a and c are in the first 2nd sector
							{
								dist_ma = distance(Rm_a, Rf_c);
								dist_mc = distance(Rm_c, Rf_c);
							}
							else//a is in the first 1st sector
							{
								dist_ma = distance(Rm_a, Rf_d);
								if(side_for_cd_mc >= 0)//c is in the first 2nd sector
								{
									dist_mc = distance(Rm_c, Rf_c);
								}
								else//c is in the first 1st sector
								{
									dist_mc = distance(Rm_c, Rf_d);
								}
							}
						}
					}//end else abcd are not in 2nd quadrant
				}
				else//between centre
				{
					//a and d are always in fourth 2nd sector
					dist_ma = distance(Rm_a, Rf_b);
					dist_md = distance(Rm_d, Rf_b);
					if(side_for_cd_mc <= 0)//c and b are in the first 1st sector
					{
						dist_mc = distance(Rm_c, Rf_d);
						dist_mb = distance(Rm_b, Rf_d);
					}
					else//c is in the first 2nd sector
					{
						dist_mc = distance(Rm_c, Rf_c);
						if(side_for_cd_mb >= 0)//b is in the first 2nd sector
						{
							dist_mb = distance(Rm_b, Rf_c);
						}
						else//b is in the first 1st sector
						{
							dist_mb = distance(Rm_b, Rf_d);
						}
					}
				}
			}
			else
				cout << "error" << endl;
		}
		fminExistDist = dist1;

		if(fminExistDist < dist2)
			fminExistDist = dist2;
		if(fminExistDist < dist3)
			fminExistDist = dist3;
		if(fminExistDist < dist4)
			fminExistDist = dist4;
		
		if(fminExistDist < dist_ma)
			fminExistDist = dist_ma;
		if(fminExistDist < dist_mb)
			fminExistDist = dist_mb;
		if(fminExistDist < dist_mc)
			fminExistDist = dist_mc;
		if(fminExistDist < dist_md)
			fminExistDist = dist_md;

		return fminExistDist;
	}
};
#endif DISTCOMPUTE_H