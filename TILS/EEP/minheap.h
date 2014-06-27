/*£¨MinHeap£©V1.0Coded by Arthur Yoo*/

#ifndef HEAP_H
#define HEAP_H
template <class type> class MinHeap{
public:
	MinHeap(int maxSize);
	MinHeap(type a[], int n);
	void Insert(type d);
	void DeleteTop();
	~MinHeap();
	void Show();
	bool Top(type*& top);
private:
	type *heapP;
	int CurSize;
	int heapMaxSize;
	void FilterDown(int start);
	void FilterUp(int start);
};

template<class type>
MinHeap<type> :: MinHeap(int maxSize){
	if(maxSize <= 0){
		cerr << "The capacity of the heap is not reasonable!" << endl;
		exit(1);
	}
	heapMaxSize = maxSize;
	heapP = new type [heapMaxSize];
	CurSize = 0;
}


template<class type>
MinHeap<type> :: MinHeap(type a[], int n){
	if(n <= 0){
		cerr << "The capacity of the heap is not reasonable!" << endl;
		exit(1);
	}
	heapMaxSize = n;
	CurSize = n;
	heapP = new type [heapMaxSize];
	pKey = new float[heapMaxSize];
	heapP = a;
	int i = (heapMaxSize - 2) / 2;
	while(i >= 0){
		FilterDown(i);
		i--;
	}
}

template<class type>
MinHeap<type> :: ~MinHeap(){
	if(CurSize != 0)
		delete [] heapP; 
}

template<class type>
void MinHeap<type> :: FilterDown(int start){
	int i, j;
	i = start;
	type temp = heapP[i];
	j = 2 * i + 1; 
	while(j <= CurSize - 1){
		if(j < CurSize - 1 && heapP[j] > heapP[j + 1]){
			j++;
		}
		if(temp <= heapP[j])
			break;
		else{
			heapP[i] = heapP[j];
			i = j;
			j = 2 * j + 1;
		}
		heapP[i] = temp;
	}
}

template<class type>
void MinHeap<type> :: FilterUp(int start){
	int i,j;
	j = start;
	type temp = heapP[j];
	i = (j - 1) / 2;
	while(j > 0){
		if(heapP[i] <= temp){
			break;
		}
		else{
			heapP[j] = heapP[i];
			j = i;
			i = (j - 1) / 2;
		}
	}
	heapP[j] = temp;
}

template<class type>
void MinHeap<type> :: Insert(type d){
	if(CurSize == heapMaxSize + 1){
		cout << "The heap is full" <<endl;
	}
	heapP[CurSize] = d;
	FilterUp(CurSize);
	CurSize++;
}

template<class type>
void MinHeap<type> :: DeleteTop(){
	if(CurSize == 0){
		return;
	}
	heapP[0] = heapP[CurSize - 1];
	CurSize--;
	FilterDown(0);
}

template<class type>
void MinHeap<type> :: Show(){
	if(CurSize == 0){
		return;
	}
	else{
		cout << heapP[0];
		for(int i = 1; i < CurSize; i++){
			cout << " " << heapP[i];
		}  
		cout << endl;
	}
}

template<class type>
bool MinHeap<type> :: Top(type*& top)
{
	if(CurSize == 0)
	{
		return false;
	}
	else
	{
		top = heapP;
	}
	return true;
}

#endif



