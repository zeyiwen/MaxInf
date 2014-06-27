#ifndef UTILITY_H
#define UTILITY_H

#include "struct.h"
#include <algorithm>
#include <functional>
#include <list>
#include <fstream>
#include <iostream>
#include <time.h>
#include "RTree_EEP.h"
#include "aRTree.h"
#include "distcompute.h"
#include "sort.h"
#include "stopwatch.h"
#include "constant.h"

using namespace std;

//Constant for array in EEP
#define MAX_SIZE_DATAM 2000000
#define MAX_SIZE_DATAF 1000000
#define MAX_SIZE_DATAC 500000

#define SORT_TIME 1
#define GETM 1
#define GETF 1

list<EntryC*> result_c2;

#endif UTILITY_H