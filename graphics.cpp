#include<iostream>
#include<graphics.h>
#include<bits/stdc++.h>
using namespace std;
void swap(int *xp, int *yp)  
{  
    int temp = *xp;  
    *xp = *yp;  
    *yp = temp;  
}  
void selectionSort(int arr[], int n)  
{
	int min_idx;    
    for (int i = 0; i < n-1; i++)  
    {  
        // Find the minimum element in unsorted array 
        min_idx = i;
        setcolor(BLACK);
		line(6*min_idx,1000,6*min_idx,1000-arr[min_idx]);
		setcolor(GREEN);
		line(6*i,1000,6*i,1000-arr[min_idx]);
		delay(50);
        for (int j = i+1; j < n; j++)  
        {setcolor(RED);
		line(6*j,1000,6*j,1000-arr[j]);
		setcolor(WHITE);
		line(6*j,1000,6*j,1000-arr[j]);
		if (arr[j] < arr[min_idx])  
            {
			min_idx = j; }
		 }
  
        // Swap the found minimum element with the first element  
		setcolor(BLACK);
		line(6*i,1000,6*i,1000-arr[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-arr[min_idx]);
		setcolor(BLACK);
		line(6*min_idx,1000,6*min_idx,1000-arr[min_idx]);
		setcolor(WHITE);
		line(6*(min_idx),1000,6*(min_idx),1000-arr[i]);
		
		swap(&arr[min_idx], &arr[i]);
	}
}
void insertionSort(int arr[], int n)  
{  
    int i, key, j;  
    for (i = 1; i < n; i++) 
    {  
        key = arr[i];
        setcolor(BLACK);
		line(6*i,1000,6*i,1000-key);
		setcolor(GREEN);
		line(6*i,1000,6*i,1000-key);  
		delay(50);
        j = i - 1;  
  
        /* Move elements of arr[0..i-1], that are  
        greater than key, to one position ahead  
        of their current position */
        while (j >= 0 && arr[j] > key) 
        {  
            
			setcolor(BLACK);
			line(6*(j+1),1000,6*(j+1),1000-arr[j+1]);
			setcolor(WHITE);
			line(6*(j+1),1000,6*(j+1),1000-arr[j]);  
            arr[j + 1] = arr[j];
			j = j - 1;  
        }  
        
		setcolor(BLACK);
		line(6*(j+1),1000,6*(j+1),1000-arr[j+1]);
		setcolor(WHITE);
		line(6*(j+1),1000,6*(j+1),1000-key);  
        arr[j + 1] = key;      
    }	  
}   
int partition (int arr[], int low, int high) 
{ 
    int pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        if (arr[j] <= pivot) 
        { 
            i++;    // increment index of smaller element 
        setcolor(BLACK);
		line(6*i,1000,6*i,1000-arr[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-arr[j]);
		setcolor(BLACK);
		line(6*j,1000,6*j,1000-arr[j]);
		setcolor(WHITE);
		line(6*(j),1000,6*j,1000-arr[i]);
		delay(10);
            swap(&arr[i], &arr[j]); 
        } 
    }
		setcolor(BLACK);
		line(6*high,1000,6*high,1000-arr[high]);
		setcolor(WHITE);
		line(6*(high),1000,6*(high),1000-arr[i+1]);
		setcolor(BLACK);
		line(6*(i+1),1000,6*(i+1),1000-arr[i+1]);
		setcolor(WHITE);
		line(6*(i+1),1000,6*(i+1),1000-arr[high]); 
    swap(&arr[i + 1], &arr[high]);
	delay(10); 
    return (i + 1); 
} 
void quickSort(int arr[], int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition(arr, low, high);
		 setcolor(GREEN);
		line(6*(pi),1000,6*(pi),1000-arr[pi]);
		setcolor(RED);
		line(6*(low),1000,6*(low),1000-arr[low]);
		setcolor(RED);
		line(6*(high),1000,6*(high),1000-arr[high]);
		delay(50); 
		
		 setcolor(WHITE);
		line(6*(pi),1000,6*(pi),1000-arr[pi]);
		line(6*(low),1000,6*(low),1000-arr[low]);
		line(6*(high),1000,6*(high),1000-arr[high]);
		
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
}
void merge(int arr[], int l, int m, int r) 
{ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
  
    /* create temp arrays */
    int L[n1], R[n2]; 
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) 
        L[i] = arr[l + i]; 
    for (j = 0; j < n2; j++) 
        R[j] = arr[m + 1+ j]; 
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2) 
    { 
        setcolor(GREEN);
		line(6*k,1000,6*k,1000-arr[k]);
		delay(10);
		if (L[i] <= R[j]) 
        { 
        	setcolor(BLACK);
			line(6*k,1000,6*k,1000-arr[k]);
			setcolor(WHITE);
			line(6*k,1000,6*k,1000-L[i]);
            arr[k] = L[i]; 
            i++; 
        } 
        else
        { 
            setcolor(BLACK);
			line(6*k,1000,6*k,1000-arr[k]);
			setcolor(WHITE);
			line(6*k,1000,6*k,1000-R[j]);
			arr[k] = R[j]; 
            j++; 
        } 
        k++; 
    } 
  
    /* Copy the remaining elements of L[], if there 
       are any */
    while (i < n1) 
    { 	setcolor(GREEN);
		line(6*k,1000,6*k,1000-arr[k]);
		delay(10);
    	setcolor(BLACK);
		line(6*k,1000,6*k,1000-arr[k]);
		setcolor(WHITE);
		line(6*k,1000,6*k,1000-L[i]);
        arr[k] = L[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
       are any */
    while (j < n2) 
    { 	setcolor(GREEN);
		line(6*k,1000,6*k,1000-arr[k]);
		delay(10);
        setcolor(BLACK);
		line(6*k,1000,6*k,1000-arr[k]);
		setcolor(WHITE);
		line(6*k,1000,6*k,1000-R[j]);
			
		arr[k] = R[j]; 
        j++; 
        k++; 
    } 
} 
  
/* l is for left index and r is right index of the 
   sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        // large l and h 
        int m = l+(r-l)/2; 
  		
        // Sort first and second halves 
        mergeSort(arr, l, m); 
        mergeSort(arr, m+1, r); 
  	delay(50);
        merge(arr, l, m, r); 
    } 
}
void heapify(int arr[], int n, int i) 
{ 
    int largest = i; // Initialize largest as root 
    int l = 2*i + 1; // left = 2*i + 1 
    int r = 2*i + 2; // right = 2*i + 2 
  
    // If left child is larger than root 
    if (l < n && arr[l] > arr[largest]) 
        largest = l; 
  
    // If right child is larger than largest so far 
    if (r < n && arr[r] > arr[largest]) 
        largest = r; 
  
    // If largest is not root 
    if (largest != i) 
    { 
    	setcolor(BLACK);
		line(6*i,1000,6*i,1000-arr[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-arr[largest]);
		setcolor(BLACK);
		line(6*largest,1000,6*largest,1000-arr[largest]);
		setcolor(WHITE);
		line(6*(largest),1000,6*(largest),1000-arr[i]);
        swap(arr[i], arr[largest]); 
  		delay(50);
        // Recursively heapify the affected sub-tree 
        heapify(arr, n, largest); 
    } 
} 
  
// main function to do heap sort 
void heapSort(int arr[], int n) 
{ 
    // Build heap (rearrange array) 
    for (int i = n / 2 - 1; i >= 0; i--) 
        heapify(arr, n, i); 
  
    // One by one extract an element from heap 
    for (int i=n-1; i>=0; i--) 
    { 
    	setcolor(BLACK);
		line(0,1000,0,1000-arr[0]);
		setcolor(WHITE);
		line(0,1000,0,1000-arr[i]);
		setcolor(BLACK);
		line(6*i,1000,6*i,1000-arr[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*i,1000-arr[0]);
        // Move current root to end 
        swap(arr[0], arr[i]); 
  
        // call max heapify on the reduced heap 
        heapify(arr, i, 0); 
    } 
}  
int getMax(int arr[], int n) 
{ 
    int mx = arr[0]; 
    for (int i = 1; i < n; i++) 
        if (arr[i] > mx) 
            mx = arr[i]; 
    return mx; 
} 
  
// A function to do counting sort of arr[] according to 
// the digit represented by exp. 
void countSort(int arr[], int n, int exp) 
{ 
    int output[n]; // output array 
    int i, count[10] = {0}; 
  
    // Store count of occurrences in count[] 
    for (i = 0; i < n; i++) 
        count[ (arr[i]/exp)%10 ]++; 
  
    // Change count[i] so that count[i] now contains actual 
    //  position of this digit in output[] 
    for (i = 1; i < 10; i++) 
        count[i] += count[i - 1]; 
  
    // Build the output array 
    for (i = n - 1; i >= 0; i--) 
    { 
        output[count[ (arr[i]/exp)%10 ] - 1] = arr[i]; 
        count[ (arr[i]/exp)%10 ]--; 
    } 
  
    // Copy the output array to arr[], so that arr[] now 
    // contains sorted numbers according to current digit 
    for (i = 0; i < n; i++) 
        {setcolor(GREEN);
		line(6*i,1000,6*i,1000-arr[i]);
		delay(50);
		setcolor(BLACK);
		line(6*i,1000,6*i,1000-arr[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*i,1000-output[i]);
		
		arr[i] = output[i];} 
} 
  
// The main function to that sorts arr[] of size n using  
// Radix Sort 
void radixsort(int arr[], int n) 
{ 
    // Find the maximum number to know number of digits 
    int m = getMax(arr, n); 
  
    // Do counting sort for every digit. Note that instead 
    // of passing digit number, exp is passed. exp is 10^i 
    // where i is current digit number 
    for (int exp = 1; m/exp > 0; exp *= 10) 
        countSort(arr, n, exp); 
}
/* Implementation of introsort*/
void swapValue(int *a, int *b) 
{ 
    int *temp = a; 
    a = b; 
    b = temp; 
    return; 
} 
/* Function to sort an array using insertion sort*/
void InsertionSort(int arr[], int *begin, int *end) 
{ 
    // Get the left and the right index of the subarray 
    // to be sorted 
    int left = begin - arr; 
    int right = end - arr; 
  
    for (int i = left+1; i <= right; i++) 
    { 
        int key = arr[i]; 
        int j = i-1; 
  
       /* Move elements of arr[0..i-1], that are 
          greater than key, to one position ahead 
          of their current position */
        while (j >= left && arr[j] > key) 
        { 
        	setcolor(BLACK);
			line(6*(j+1),1000,6*(j+1),1000-arr[j+1]);
			setcolor(WHITE);
			line(6*(j+1),1000,6*(j+1),1000-arr[j]); 
            arr[j+1] = arr[j]; 
            j = j-1; 
        } 
        setcolor(BLACK);
		line(6*(j+1),1000,6*(j+1),1000-arr[j+1]);
		setcolor(WHITE);
		line(6*(j+1),1000,6*(j+1),1000-key);  
        arr[j+1] = key; 
   } 
  
   return; 
} 
int* Partition(int arr[], int low, int high) 
{ 
    int pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        if (arr[j] <= pivot) 
        { 
            // increment index of smaller element 
            i++; 
  			setcolor(BLACK);
		line(6*i,1000,6*i,1000-arr[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-arr[j]);
		setcolor(BLACK);
		line(6*j,1000,6*j,1000-arr[j]);
		setcolor(WHITE);
		line(6*(j),1000,6*j,1000-arr[i]);
		delay(10);
            swap(arr[i], arr[j]); 
        } 
    } 
    setcolor(BLACK);
		line(6*high,1000,6*high,1000-arr[high]);
		setcolor(WHITE);
		line(6*(high),1000,6*(high),1000-arr[i+1]);
		setcolor(BLACK);
		line(6*(i+1),1000,6*(i+1),1000-arr[i+1]);
		setcolor(WHITE);
		line(6*(i+1),1000,6*(i+1),1000-arr[high]); 
    swap(arr[i + 1], arr[high]); 
    return (arr + i + 1); 
} 
// A function that find the middle of the 
// values pointed by the pointers a, b, c 
// and return that pointer 
int *MedianOfThree(int * a, int * b, int * c) 
{ 
    if (*a < *b && *b < *c) 
        return (b); 
  
    if (*a < *c && *c <= *b) 
        return (c); 
  
    if (*b <= *a && *a < *c) 
        return (a); 
  
    if (*b < *c && *c <= *a) 
        return (c); 
  
    if (*c <= *a && *a < *b) 
        return (a); 
  
    if (*c <= *b && *b <= *a) 
        return (b); 
} 
  
// A Utility function to perform intro sort 
void IntrosortUtil(int arr[], int * begin, 
                  int * end, int depthLimit) 
{ 
    // Count the number of elements 
    int size = end - begin; 
  
      // If partition size is low then do insertion sort 
    if (size < 16) 
    { 
        InsertionSort(arr, begin, end); 
        return; 
    } 
  
    // If the depth is zero use heapsort 
    if (depthLimit == 0) 
    { 
        make_heap(begin, end+1); 
        sort_heap(begin, end+1); 
        return; 
    } 
  
    // Else use a median-of-three concept to 
    // find a good pivot 
    int * pivot = MedianOfThree(begin, begin+size/2, end); 
  
    // Swap the values pointed by the two pointers 
    swapValue(pivot, end); 
  
   // Perform Quick Sort 
    int * partitionPoint = Partition(arr, begin-arr, end-arr); 
    IntrosortUtil(arr, begin, partitionPoint-1, depthLimit - 1); 
    IntrosortUtil(arr, partitionPoint + 1, end, depthLimit - 1); 
  
    return; 
} 
void Introsort(int arr[], int *begin, int *end) 
{ 
    int depthLimit = 2 * log(end-begin); 
  
    // Perform a recursive Introsort 
    IntrosortUtil(arr, begin, end, depthLimit); 
  
      return; 
}  

void shellSort(int arr[], int n) 
{ 
    // Start with a big gap, then reduce the gap 
    for (int gap = n/2; gap > 0; gap /= 2) 
    { 
        // Do a gapped insertion sort for this gap size. 
        // The first gap elements a[0..gap-1] are already in gapped order 
        // keep adding one more element until the entire array is 
        // gap sorted  
        for (int i = gap; i < n; i += 1) 
        { 
            // add a[i] to the elements that have been gap sorted 
            // save a[i] in temp and make a hole at position i 
            int temp = arr[i]; 
  
            // shift earlier gap-sorted elements up until the correct  
            // location for a[i] is found 
            int j;             
            for (j = i; j >= gap && arr[j - gap] > temp; j -= gap) 
                {
		setcolor(BLACK);
		line(6*j,1000,6*j,1000-arr[j]);
		setcolor(WHITE);
		line(6*(j),1000,6*j,1000-arr[j-gap]);
				arr[j] = arr[j - gap]; }
              
            //  put temp (the original a[i]) in its correct location 
            		setcolor(BLACK);
		line(6*j,1000,6*j,1000-arr[j]);
		setcolor(WHITE);
		line(6*(j),1000,6*j,1000-temp);
            arr[j] = temp; 
            delay(10);
        } 
    } 
} 

void bubbleSort(int arr[], int n) 
{ 
   int i, j; 
   for (i = 0; i < n-1; i++)       
  
       // Last i elements are already in place    
       for (j = 0; j < n-i-1; j++)  
           if (arr[j] > arr[j+1]) 
              {setcolor(BLACK);
		line(6*j,1000,6*j,1000-arr[j]);
		setcolor(WHITE);
		line(6*(j),1000,6*(j),1000-arr[j+1]);
		setcolor(BLACK);
		line(6*(j+1),1000,6*(j+1),1000-arr[j+1]);
		setcolor(WHITE);
		line(6*(j+1),1000,6*(j+1),1000-arr[j]); 
			  swap(&arr[j], &arr[j+1]);} 
} 

void CocktailSort(int a[], int n) 
{ 
    bool swapped = true; 
    int start = 0; 
    int end = n - 1; 
  
    while (swapped) { 
        // reset the swapped flag on entering 
        // the loop, because it might be true from 
        // a previous iteration. 
        swapped = false; 
  
        // loop from left to right same as 
        // the bubble sort 
        for (int i = start; i < end; ++i) { 
            if (a[i] > a[i + 1]) {
				 setcolor(BLACK);
		line(6*i,1000,6*i,1000-a[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-a[i+1]);
		setcolor(BLACK);
		line(6*(i+1),1000,6*(i+1),1000-a[i+1]);
		setcolor(WHITE);
		line(6*(i+1),1000,6*(i+1),1000-a[i]);
                swap(a[i], a[i + 1]); 
                swapped = true; 
            } 
        } 
  
        // if nothing moved, then array is sorted. 
        if (!swapped) 
            break; 
  
        // otherwise, reset the swapped flag so that it 
        // can be used in the next stage 
        swapped = false; 
  
        // move the end point back by one, because 
        // item at the end is in its rightful spot 
        --end; 
  
        // from right to left, doing the 
        // same comparison as in the previous stage 
        for (int i = end - 1; i >= start; --i) { 
            if (a[i] > a[i + 1]) { 
            setcolor(BLACK);
		line(6*i,1000,6*i,1000-a[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-a[i+1]);
		setcolor(BLACK);
		line(6*(i+1),1000,6*(i+1),1000-a[i+1]);
		setcolor(WHITE);
		line(6*(i+1),1000,6*(i+1),1000-a[i]);
                swap(a[i], a[i + 1]);
                swapped = true; 
            } 
        } 
  
        // increase the starting point, because 
        // the last stage would have moved the next 
        // smallest number to its rightful spot. 
        ++start; 
    } 
} 
void gnomeSort(int arr[], int n) 
{ 
    int index = 0; 
  
    while (index < n) { 
        if (index == 0) 
            index++; 
        if (arr[index] >= arr[index - 1]) 
            index++; 
        else { 
        	setcolor(BLACK);
		line(6*index,1000,6*index,1000-arr[index]);
		setcolor(WHITE);
		line(6*(index),1000,6*(index),1000-arr[index-1]);
		setcolor(BLACK);
		line(6*(index-1),1000,6*(index-1),1000-arr[index-1]);
		setcolor(WHITE);
		line(6*(index-1),1000,6*(index-1),1000-arr[index]);
            swap(arr[index], arr[index - 1]); 
            index--; 
        } 
    } 
    return; 
}
/*The parameter dir indicates the sorting direction, ASCENDING 
   or DESCENDING; if (a[i] > a[j]) agrees with the direction, 
   then a[i] and a[j] are interchanged.*/
void compAndSwap(int a[], int i, int j, int dir) 
{ 
    if (dir==(a[i]>a[j])) 
        {setcolor(BLACK);
		line(6*i,1000,6*i,1000-a[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-a[j]);
		setcolor(BLACK);
		line(6*(j),1000,6*(j),1000-a[j]);
		setcolor(WHITE);
		line(6*(j),1000,6*(j),1000-a[i]);
		delay(20);
		swap(a[i],a[j]); }
} 
  
/*It recursively sorts a bitonic sequence in ascending order, 
  if dir = 1, and in descending order otherwise (means dir=0). 
  The sequence to be sorted starts at index position low, 
  the parameter cnt is the number of elements to be sorted.*/
void bitonicMerge(int a[], int low, int cnt, int dir) 
{ 
    if (cnt>1) 
    { 
        int k = cnt/2; 
        for (int i=low; i<low+k; i++) 
            compAndSwap(a, i, i+k, dir); 
        bitonicMerge(a, low, k, dir); 
        bitonicMerge(a, low+k, k, dir); 
    } 
} 
  
/* This function first produces a bitonic sequence by recursively 
    sorting its two halves in opposite sorting orders, and then 
    calls bitonicMerge to make them in the same order */
void bitonicSort(int a[],int low, int cnt, int dir) 
{ 
    if (cnt>1) 
    { 
        int k = cnt/2; 
  
        // sort in ascending order since dir here is 1 
        bitonicSort(a, low, k, 1); 
  
        // sort in descending order since dir here is 0 
        bitonicSort(a, low+k, k, 0); 
  
        // Will merge wole sequence in ascending order 
        // since dir=1. 
        bitonicMerge(a,low, cnt, dir); 
    } 
}
// To check if array is sorted or not 
bool isSorted(int a[], int n) 
{ 
    while ( --n > 1 ) 
        if (a[n] < a[n-1]) 
            return false; 
    return true; 
} 
  
// To generate permuatation of the array 
void shuffle(int a[], int n) 
{ 
    for (int i=0; i < n; i++) 
        {int x=rand()%n;
		setcolor(BLACK);
		line(6*i,1000,6*i,1000-a[i]);
		setcolor(WHITE);
		line(6*(i),1000,6*(i),1000-a[x]);
		setcolor(BLACK);
		line(6*(x),1000,6*(x),1000-a[x]);
		setcolor(WHITE);
		line(6*(x),1000,6*(x),1000-a[i]);
		delay(20);
		swap(a[i], a[x]);} 
} 
  
// Sorts array a[0..n-1] using Bogo sort 
void bogosort(int a[], int n) 
{ 
    // if array is not sorted then shuffle 
    // the array again 
    while ( !isSorted(a, n) ) 
        shuffle(a, n); 
}   
int main()
{
	int gd=DETECT, gm;
	initgraph(&gd,&gm,"");
	setlinestyle(0,0,5);
	int length[200],arr[200];
	for(int i=0;i<200;i++){
		length[i]=rand()%1000;arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	settextstyle(8,0,2);
	outtextxy(1200,0,"Selection Sort");
	selectionSort(arr,200);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Insertion Sort");
	insertionSort(arr,200);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Quick Sort");
	quickSort(arr,0,199);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Merge Sort");
	mergeSort(arr,0,199);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Heap Sort");
	heapSort(arr,200);
	
	getch();   
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Radix Sort");
	radixsort(arr,200);
	
	getch();   
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Intro Sort");
	Introsort(arr,arr,arr+199);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Shell Sort");
	shellSort(arr,200);
	
	getch();   
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Bubble Sort");
	bubbleSort(arr,200); 
	
	getch();   
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Cocktail Sort");
	CocktailSort(arr,200);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Gnome Sort");
	gnomeSort(arr,200);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<128;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Bitonic Sort");
	bitonicSort(arr,0,128,1);
	
	getch();
	cleardevice();
	setcolor(WHITE);
	for(int i=0;i<200;i++){
		arr[i]=length[i];
		line(6*i,1000,6*i,1000-length[i]);
	}
	getch();
	outtextxy(1200,0,"Bogo Sort");
	bogosort(arr,200);
	
	getch();
	closegraph();
}

