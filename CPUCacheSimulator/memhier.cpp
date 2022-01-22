#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <string_view>

using namespace std;

int TLBNumSets;
int TLBSetSize;
int TLBIndexSize;

int numVPages;
int numPPages;
int pageSize;
int pageTableIndexSize,pageOffset;

int l1Numsets;
int l1SetSize;
int l1LineSize;
bool l1WT;
int l1IndexSize,l1Offset;

int l2Numsets;
int l2SetSize;
int l2LineSize;
bool l2WT;
int l2IndexSize,l2Offset;

bool virtAdd;
bool TLBenabled;
bool l2;

int dtlbhits = 0;
int dtlbmisses = 0;
int pthits = 0;
int ptmisses = 0;
int l1hits = 0;
int l1misses = 0;
int l2hits = 0;
int l2misses = 0;

int totalReads = 0;
int totalWrites = 0;

int mainmemoryrefs = 0;
int pagetablerefs = 0;
int diskrefs = 0;



typedef struct cacheEntry{
	unsigned long tag;
	unsigned long LRU;
	bool dirty;
	unsigned long address;
};

typedef struct pageTableEntry{
	unsigned long virt;
	unsigned long LRU;
};

pageTableEntry *PAGETABLE;
cacheEntry **L1CACHE;
cacheEntry **L2CACHE;
cacheEntry **TLB;

//Prints a formated line that contains all the information of a data acquisition
void printLine(unsigned long *VA, unsigned long *VPNum, unsigned long *PO, unsigned long *TLBTag, unsigned long *TLBIndex, bool *TLBHit, bool *PTHit, unsigned long *PhysPageNum, unsigned long *DCTag, unsigned long *DCIndex, bool DCHit, unsigned long *L2Tag, unsigned long *L2Index, bool *L2Hit, bool L2Relevent){
	printf("%08x ", *VA);

	if(VPNum != NULL){
		printf("%6x ", *VPNum);
	}else{
		printf("       ");
	}

	printf("%4x ", *PO);

	if(TLBTag != NULL){
		printf("%6x ", *TLBTag);
	}else{
		printf("       ");
	}

	if(TLBIndex != NULL){
		printf("%3x ", *TLBIndex);
	}else{
		printf("    ");
	}
	
	if(TLBHit == NULL){
		printf("     ");
	}else if(*TLBHit){
		printf("hit  ");
	}else{
		printf("miss ");
	}

	if(PTHit == NULL){
		printf("     ");
	}else if(*PTHit){
		printf("hit  ");
	}else{
		printf("miss ");
	}

	printf("%4x %6x %3x ",*PhysPageNum,*DCTag,*DCIndex);
	if(DCHit){
		printf("hit  ");
	}else{
		printf("miss ");
	}
	
	if(l2 && L2Relevent){
		if(L2Tag != NULL){
			printf("%6x ", *L2Tag);
		}else{
			printf("       ");
		}

		if(L2Index != NULL){
			printf("%3x ", *L2Index);
		}else{
			printf("    ");
		}

		if(L2Hit == NULL){
			printf("    ");
		}else if(*L2Hit){
			printf("hit ");
		}else{
			printf("miss");
		}
		
	}else{
		
	}
	printf("\n");
	
}

//Prints Statistics after the simulation is done
void printStatistics(){
	printf("\nSimulation statistics\n\n");
	printf("dtlb hits        : %d\n",dtlbhits);
	printf("dtlb misses      : %d\n",dtlbmisses);

	if(dtlbmisses+dtlbhits != 0){
		printf("dtlb hit ratio   : %7.6f\n",((double)dtlbhits)/((double)(dtlbhits+dtlbmisses)));
	}else{
		printf("dtlb hit ratio   : N/A\n");
	}
	
	printf("\npt hits          : %d\n",pthits);
	printf("pt faults        : %d\n",ptmisses);

	if(ptmisses+pthits != 0){
		printf("pt hit ratio     : %7.6f\n",((double)pthits)/((double)(pthits+ptmisses)));
	}else{
		printf("pt hit ratio     : N/A\n");
	}

	printf("\ndc hits          : %d\n",l1hits);
	printf("dc misses        : %d\n",l1misses);

	if(l1misses+l1hits != 0){
		printf("dc hit ratio     : %7.6f\n",((double)l1hits)/((double)(l1hits+l1misses)));
	}else{
		printf("dc hit ratio     : N/A\n");
	}

	printf("\nL2 hits          : %d\n",l2hits);
	printf("L2 misses        : %d\n",l2misses);

	if(l2misses+l2hits != 0){
		printf("L2 hit ratio     : %7.6f\n",((double)l2hits)/((double)(l2hits+l2misses)));
	}else{
		printf("L2 hit ratio     : N/A\n");
	}

	printf("\nTotal reads      : %d\n",totalReads);
	printf("Total writes     : %d\n",totalWrites);
	printf("Ratio of reads   : %7.6f\n\n",((double)totalReads)/((double)(totalReads+totalWrites)));

	printf("main memory refs : %d\n",mainmemoryrefs);
	printf("page table refs  : %d\n",pagetablerefs);
	printf("disk refs        : %d\n",diskrefs);

}

//Reads the config
void readConfig(){
	string tempS;
	int tempI;
	ifstream config("trace.config");

	config >> tempS >> tempS >> tempS >> tempS >> tempS >> tempS >> tempI;
	TLBNumSets = tempI;
	config >> tempS >> tempS >> tempI;
	TLBSetSize = tempI;

	config >> tempS >> tempS >> tempS >> tempS >> tempS >> tempS >> tempS >> tempI;
	numVPages = tempI;
	config >> tempS >> tempS >> tempS >> tempS >> tempI;
	numPPages = tempI;
	config >> tempS >> tempS >> tempI;
	pageSize = tempI;

	config >> tempS >> tempS >> tempS >> tempS >> tempS >> tempS >> tempI;
	l1Numsets = tempI;
	config >> tempS >> tempS >> tempI;
	l1SetSize = tempI;
	config >> tempS >> tempS >> tempI;
	l1LineSize = tempI;
	config >> tempS >> tempS >> tempS >> tempS >> tempS;
	if(tempS == "n"){
		l1WT = false;
	}else{
		l1WT = true;
	}

	config >> tempS >> tempS >> tempS >> tempS >> tempS >> tempS >> tempI;
	l2Numsets = tempI;
	config >> tempS >> tempS >> tempI;
	l2SetSize = tempI;
	config >> tempS >> tempS >> tempI;
	l2LineSize = tempI;
	config >> tempS >> tempS >> tempS >> tempS >> tempS;
	if(tempS == "n"){
		l2WT = false;
	}else{
		l2WT = true;
	}

	config >> tempS >> tempS >> tempS;
	if(tempS == "n"){
		virtAdd = false;
	}else{
		virtAdd = true;
	}
	config >> tempS >> tempS;
	if(tempS == "n"){
		TLBenabled = false;
	}else{
		TLBenabled = true;
	}
	config >> tempS >> tempS >> tempS;
	if(tempS == "n"){
		l2 = false;
	}else{
		l2 = true;
	}

	TLBIndexSize = log2(TLBNumSets);
	pageTableIndexSize = log2(numVPages);
	pageOffset = log2(pageSize);
	l1IndexSize = log2(l1Numsets);
	l1Offset = log2(l1LineSize);
	l2IndexSize = log2(l2Numsets);
	l2Offset = log2(l2LineSize);



	printf("Data TLB contains %d sets.\nEach set contains %d entries.\nNumber of bits used for the index is %d.\n\n",TLBNumSets,TLBSetSize,TLBIndexSize);
	printf("Number of virtual pages is %d.\nNumber of physical pages is %d.\nEach page contains %d bytes.\nNumber of bits used for the page table index is %d.\nNumber of bits used for the page offset is %d.\n\n",numVPages,numPPages,pageSize,pageTableIndexSize,pageOffset);
	
	printf("D-cache contains %d sets.\nEach set contains %d entries.\nEach line is %d bytes.\n", l1Numsets, l1SetSize, l1LineSize); 
	if(l1WT){
		printf("The cache uses a no write-allocate and write-through policy.\n");
	}else{
		printf("The cache uses a write-allocate and write-back policy.\n");
	}
	printf("Number of bits used for the index is %d.\nNumber of bits used for the offset is %d.\n\n",l1IndexSize,l1Offset);

	printf("L2-cache contains %d sets.\nEach set contains %d entries.\nEach line is %d bytes.\n", l2Numsets, l2SetSize, l2LineSize); 
	if(l2WT){
		printf("The cache uses a no write-allocate and write-through policy.\n");
	}else{
		printf("The cache uses a write-allocate and write-back policy.\n");
	}
	printf("Number of bits used for the index is %d.\nNumber of bits used for the offset is %d.\n\n",l2IndexSize,l2Offset);

	if(virtAdd){
		printf("The addresses read in are virtual addresses.\n");
	}else{
		printf("The addresses read in are physical addresses.\n"); //mabye
	}

	if(!virtAdd && TLBenabled){
		printf("hierarchy: TLB cannot be enabled when virtual addresses are disabled\n");
		exit(1);
	}
	if(!TLBenabled){
		printf("TLB is disabled in this configuration.\n");
	}
	if(!l2){
		printf("L2 cache is disabled in this configuration.\n");
	}
	printf("\n");
}

//Gets tag and index for DC information from address 
void l1Info(unsigned long *address,unsigned long *l1tag, unsigned long *l1index){
	unsigned long mask = 0xFFFFFFFF;
	*l1index  = (*address >> l1Offset )& (mask >> (32 - l1IndexSize));
	*l1tag = (*address >> (l1Offset + l1IndexSize));
}

//Gets tag and index for L2 information from address 
void l2Info(unsigned long *address,unsigned long *l2tag, unsigned long *l2index){
	unsigned long mask = 0xFFFFFFFF;
	*l2index  = (*address >> l2Offset )& (mask >> (32 - l2IndexSize));
	*l2tag = (*address >> (l2Offset + l2IndexSize));
}

//Allocates memory for level 1 cache and initializes tag with -1
void makeL1Cache(){
	L1CACHE = (cacheEntry **)calloc(l1Numsets,sizeof(cacheEntry *));
	for(int i = 0; i < l1Numsets; i++){
		L1CACHE[i] = (cacheEntry *)calloc(l1SetSize, sizeof(cacheEntry));
		for(int j = 0; j < l1SetSize; j++){
			L1CACHE[i][j].tag = -1;
		}
	}
}

//Debugging Function
void printL1Cache(){
	printf("\nL1\n");
	printf("--index-|--tag---|---LRU--|--adr---|-dirty--\n");
	printf("--------|--------|--------|--------|--------\n");
	for(int i = 0; i < l1Numsets; i++){
		for(int j = 0; j < l1SetSize; j++){
			printf("--------|--------|--------|--------|--------\n");
			printf("%8x|%8x|%8x|%8x|%8x\n",i,L1CACHE[i][j].tag,L1CACHE[i][j].LRU,L1CACHE[i][j].address,L1CACHE[i][j].dirty);
		}
		printf("--------|--------|--------|--------|--------\n");
		printf("--------|--------|--------|--------|--------\n");
	}
}

//Allocates memory for level 2 cache and initializes tag with -1
void makeL2Cache(){
	L2CACHE = (cacheEntry **)calloc(l2Numsets,sizeof(cacheEntry *));
	for(int i = 0; i < l2Numsets; i++){
		L2CACHE[i] = (cacheEntry *)calloc(l2SetSize, sizeof(cacheEntry));
		for(int j = 0; j < l2SetSize; j++){
			L2CACHE[i][j].tag = -1;
		}
	}
}

//Debugging Function
void printL2Cache(){
	printf("\nL2\n");
	printf("\t\t\t\t--index-|--tag---|---LRU--|--adr---|-dirty--\n");
	printf("\t\t\t\t--------|--------|--------|--------|--------\n");
	for(int i = 0; i < l2Numsets; i++){
		for(int j = 0; j < l2SetSize; j++){
			printf("\t\t\t\t--------|--------|--------|--------|--------\n");
			printf("\t\t\t\t%8x|%8x|%8x|%8x|%8x\n",i,L2CACHE[i][j].tag,L2CACHE[i][j].LRU,L2CACHE[i][j].address,L2CACHE[i][j].dirty);
		}
		printf("\t\t\t\t--------|--------|--------|--------|--------\n");
		printf("\t\t\t\t--------|--------|--------|--------|--------\n");
	}
}

//Allocates memory for the tlb, initializes tag with -1 and ensures the least recently used bit starts at 0
void makeTLB(){
	TLB = (cacheEntry **)calloc(TLBNumSets,sizeof(cacheEntry *));
	for(int i = 0; i < TLBNumSets; i++){
		TLB[i] = (cacheEntry *)calloc(TLBSetSize, sizeof(cacheEntry));
		for(int j = 0; j < TLBSetSize; j++){
			TLB[i][j].tag = -1;
			TLB[i][j].LRU = 0;
		}
	}
}

//Debugging Function
void printTLB(){
	printf("\nTLB\n");
	printf("\t\t\t\t--index-|--tag---|---LRU--|--adr---|-dirty--\n");
	printf("\t\t\t\t--------|--------|--------|--------|--------\n");
	for(int i = 0; i < TLBNumSets; i++){
		for(int j = 0; j < TLBSetSize; j++){
			printf("\t\t\t\t--------|--------|--------|--------|--------\n");
			printf("\t\t\t\t%8x|%8x|%8x|%8x|%8x\n",i,TLB[i][j].tag,TLB[i][j].LRU,TLB[i][j].address,TLB[i][j].dirty);
		}
		printf("\t\t\t\t--------|--------|--------|--------|--------\n");
		printf("\t\t\t\t--------|--------|--------|--------|--------\n");
	}
}

//Invalidates a DC entry based off of L2 being replaced
void invalidateL1(unsigned long address){
	for(int i = 0; i < l1Numsets; i++){
		for(int j=0; j < l1SetSize; j++){
			//Find cache that we are invalidating
			if(L1CACHE[i][j].address == address){
				if(L1CACHE[i][j].dirty){
					//If its dirty and the l2 isn't enabled then memory would need to be accessed
					if(!l2){
						mainmemoryrefs++;
					}
				}
				//printf("invalidating DC line with tag %x and index %x and adr %x since l2 being replaced \n", L1CACHE[i][j].tag, i, L1CACHE[i][j].address);
				L1CACHE[i][j].tag = -1;
				L1CACHE[i][j].LRU = 0;
				L1CACHE[i][j].address = 0;
				L1CACHE[i][j].dirty = 0;
			}
		}
	}
}

//Puts physical address in corrosponding tlb entry
void updateTLB(unsigned long address,unsigned long ppagenum){
	
	unsigned long mask = 0xFFFFFFFF;
	unsigned long temp = (address >> pageOffset);
	unsigned long lookIndex = temp & (mask >> (32 - TLBIndexSize)); //Indicates which tlb set to look in
	unsigned long tag = temp >> TLBIndexSize;
	int maxLRU = TLB[lookIndex][0].LRU;
	int minLRU = TLB[lookIndex][0].LRU;
	
	//printf("Updating TLB with address %x, tag %x, at index %x\n",address,tag,lookIndex);
	//printf("current maxLRU = %d\n",maxLRU);
	int maxLRUIndex = 0;
	int minLRUIndex = 0;

	//Search through tlb set
	for(int i=0; i < TLBSetSize; i++){
		//printf("setoffset + i*2 :%d\n",setoffset + i*2);	 
		if(TLB[lookIndex][i].tag == tag){
			//
		}
		if(TLB[lookIndex][i].LRU > maxLRU){
			//printf("Found lru %d\n",TLB[lookIndex][i].LRU);
			maxLRU = TLB[lookIndex][i].LRU;
			maxLRUIndex = i;
		}
		if(TLB[lookIndex][i].LRU < minLRU){
			minLRU = TLB[lookIndex][i].LRU;
			minLRUIndex = i;
		}
	}
	//Update the least recently used tlb entry with the new phyiscal address
	TLB[lookIndex][minLRUIndex].tag = tag;
	TLB[lookIndex][minLRUIndex].address = ppagenum;
	TLB[lookIndex][minLRUIndex].LRU = maxLRU+1;
}

//Handles reads and writes to L2 cache
bool accessL2Cache(unsigned long address, char RW, unsigned long index, unsigned long tag){
	int i;
	int setoffset = 2*l2SetSize*index;
	//printf("Setoffset :%d\n",setoffset);
	
	bool hit = false;
	int hitIndex = 0;
	int maxLRU = L2CACHE[index][0].LRU;
	int minLRU = L2CACHE[index][0].LRU;
	
	int maxLRUIndex = 0;
	int minLRUIndex = 0;

	//Search through l2 set
	for(i=0; i < l2SetSize; i++){
		//printf("setoffset + i*2 :%d\n",setoffset + i*2);	 
		if(L2CACHE[index][i].tag == tag){
			hit = true;
			hitIndex = i;
		}
		if(L2CACHE[index][i].LRU > maxLRU){
			maxLRU = L2CACHE[index][i].LRU;
			
			maxLRUIndex = i;
		}
		if(L2CACHE[index][i].LRU < minLRU){
			minLRU = L2CACHE[index][i].LRU;
			minLRUIndex = i;
		}
	}
	

	if(hit){ //On cache hit update the least recently used 
		L2CACHE[index][hitIndex].LRU = maxLRU+1;
		l2hits++;
		if(RW=='W'){ //If it was a write hit then make that block dirty
			L2CACHE[index][hitIndex].dirty = true;
		}
		return true;
	}else{
		l2misses++;
		if(!l2WT){//If in write back mode
			if(L2CACHE[index][minLRUIndex].dirty){ //write back to memory
				mainmemoryrefs++;
				//printf("writing back L2 %x with tag %x and index %x\n",L2CACHE[index][minLRUIndex].address,L2CACHE[index][minLRUIndex].tag,index);
				invalidateL1(L2CACHE[index][minLRUIndex].address);
			}
			mainmemoryrefs++;
		}

		//Since we missed we need to update a block in the cache with the new data
		L2CACHE[index][minLRUIndex].tag = tag;
		L2CACHE[index][minLRUIndex].address = address;
		L2CACHE[index][minLRUIndex].LRU = maxLRU+1;
		if(RW=='W'){
			L2CACHE[index][minLRUIndex].dirty = true;
		}
		if(RW=='R'){
			L2CACHE[index][minLRUIndex].dirty = false;
		}
		return false;
	}
	
	
	
}

//Handles accesses to DC cache and will pass off to l2 on a miss.
bool accessL1Cache(unsigned long address, char RW,unsigned long index, unsigned long tag, bool *didL2Hit, bool *wasL2Accesses){
	
	//printf("Accessing l1 with address %x\n",address);
	int i;
	int setoffset = 2*l1SetSize*index;
	
	
	bool hit = false;
	int hitIndex = 0;
	int maxLRU = L1CACHE[index][0].LRU;
	int minLRU = L1CACHE[index][0].LRU;
	
	int maxLRUIndex = 0;
	int minLRUIndex = 0;

	//Search through l1 cache
	for(i=0; i < l1SetSize; i++){
		if(L1CACHE[index][i].tag == tag){
			hit = true;
			hitIndex = i;
		}
		if(L1CACHE[index][i].LRU > maxLRU){
			maxLRU = L1CACHE[index][i].LRU;
			
			maxLRUIndex = i;
		}
		if(L1CACHE[index][i].LRU < minLRU){
			minLRU = L1CACHE[index][i].LRU;
			minLRUIndex = i;
		}
	}
	

	if(hit){ //On cache hit update the least recently used 
		L1CACHE[index][hitIndex].LRU = maxLRU+1;
		L1CACHE[index][hitIndex].address = address; //not sure
		l1hits++;
		if(RW=='W'){
			L1CACHE[index][hitIndex].dirty = true;
		}
		return true;
	}else{
		l1misses++;
		if(!l1WT){
			if(L1CACHE[index][minLRUIndex].dirty){//Write back to memory or L2
				//printf("writing back tag %x from index %d\n",L1CACHE[index][minLRUIndex].tag , index);
				if(l2){
					unsigned long debugtag,debugindex;
					l1Info(&(L1CACHE[index][minLRUIndex].address),&debugtag,&debugindex);
					unsigned long l2tag,l2index;
					l2Info(&(L1CACHE[index][minLRUIndex].address),&l2tag,&l2index);
					//printf("writing back DC %x with tag %x and index %x\n",L1CACHE[index][minLRUIndex].address,debugtag,debugindex);
					accessL2Cache(L1CACHE[index][minLRUIndex].address,'W',l2index,l2tag);
				}else{
					mainmemoryrefs++;
				}
			}

			if(l2){//If level 2 cache is enabled 
					//read to l2
					unsigned long l2tag,l2index;
					l2Info(&address,&l2tag,&l2index);
					*didL2Hit = accessL2Cache(address,RW,l2index,l2tag);
					*wasL2Accesses = true;
			}else{
				mainmemoryrefs++;
			}
		}
			
		//Since we missed we need to update a block in the cache with the new data
		L1CACHE[index][minLRUIndex].tag = tag;
		L1CACHE[index][minLRUIndex].address = address;
		L1CACHE[index][minLRUIndex].LRU = maxLRU+1;
		if(RW=='W'){
			L1CACHE[index][minLRUIndex].dirty = true;
		}
		if(RW=='R'){
			L1CACHE[index][minLRUIndex].dirty = false;
		}
		return false;
	}
}

//Allocates memory for the page table, initializes the virtual address with -1 and ensures the least recently used bit starts at 0
void makePageTable(){
	PAGETABLE = (pageTableEntry *)malloc(sizeof(pageTableEntry) * numPPages);
	for(int i = 0; i < numPPages; i++){
		PAGETABLE[i].virt = -1;
		PAGETABLE[i].LRU = 0;
	}
}

//Debugging Function
void printPageTable(){
	printf("\nPAGETABLE\n");
	printf("index|--virt--|LRU\n");
	for(int i = 0; i < numPPages; i++){
		printf("%5x|%8x|%3x\n",i,PAGETABLE[i].virt,PAGETABLE[i].LRU);
	}
}


//Invalidates entries after a page fault
void invalidate(unsigned long badPPageNum){
	//Invalidates tlb entries
	if(TLBenabled){
		for(int i = 0; i < TLBNumSets; i++){
			for(int j=0; j < TLBSetSize; j++){
			
				if(TLB[i][j].address == badPPageNum){
					//printf("removing bad page %x from TbL\n",badPPageNum);
					
					TLB[i][j].tag = -1;
					TLB[i][j].LRU = 0;
					TLB[i][j].address = 0;
					TLB[i][j].dirty = 0;
				}
			}
		}
	}
	//Invalidates l1 entries 
	for(int i = 0; i < l1Numsets; i++){
		for(int j=0; j < l1SetSize; j++){
			if((L1CACHE[i][j].address >> pageOffset) == badPPageNum){
				if(L1CACHE[i][j].dirty){
					if(l2){
						//l2hits++;
					}else{
						mainmemoryrefs++;
					}
				}
				//printf("invalidating DC line with tag %x and index %x and adr %x since page %x is gone \n", L1CACHE[i][j].tag, i, L1CACHE[i][j].address,badPPageNum);
				L1CACHE[i][j].tag = -1;
				L1CACHE[i][j].LRU = 0;
				L1CACHE[i][j].address = 0;
				L1CACHE[i][j].dirty = 0;
				//printf("removing bad l1 entry\n");
			}
		}
	}
	
	//Invalidates l2 entries 
	if(l2){
		for(int i = 0; i < l2Numsets; i++){
			for(int j=0; j < l2SetSize; j++){
				//printf("%x,%x address is %x\n",i,j,L2CACHE[i][j].address);
				if((L2CACHE[i][j].address >> pageOffset) == badPPageNum){
					//printf("s bad l2 with tag %x address %x entry\n",L2CACHE[i][j].tag,L2CACHE[i][j].address);
					if(L2CACHE[i][j].dirty){
						mainmemoryrefs++;
						diskrefs++;	
					}
					
					L2CACHE[i][j].tag = -1;
					L2CACHE[i][j].LRU = 0;
					L2CACHE[i][j].address = 0;
					L2CACHE[i][j].dirty = 0;
				}
			}
		}
	}
}

//returns the physical page given virtual adress
unsigned long pageTableLookup(unsigned long address,bool* didPageTableHit){
	pagetablerefs++;
	unsigned long pagenum = (address >> pageOffset);
	unsigned long index = -1;
	unsigned long LRU = PAGETABLE[0].LRU;
	unsigned long LRUindex = 0;
	unsigned long maxLRU = PAGETABLE[0].LRU;
	unsigned long maxLRUindex = 0;

	//Search through the page table
	for(int i = 0; i < numPPages; i++){
		if(PAGETABLE[i].virt == pagenum){
			index = i;
		}
		if(PAGETABLE[i].LRU < LRU){
			LRU = PAGETABLE[i].LRU;
			LRUindex=i;
		}
		if(PAGETABLE[i].LRU > maxLRU){
			maxLRU = PAGETABLE[i].LRU;
			maxLRUindex=i;
		}
	}

	if(index != -1){//We got a hit on the page table
		pthits++;
		*didPageTableHit = true;
		PAGETABLE[index].LRU = maxLRU+1;
		if(TLBenabled){
			updateTLB(address,index);
		}
		return index;
	}else{ //If index is still -1 then we missed in the page table
		if(LRU != 0){ //If the least recently used is not zero then that means we are replacing a page table entry and need to invalidate values in the cache corresponding to the replaced entry
			//printf("Evicting a physical block %x\n", LRUindex);
			invalidate(LRUindex);
		}
		ptmisses++;
		diskrefs++;
		*didPageTableHit = false;
		PAGETABLE[LRUindex].virt = pagenum;
		PAGETABLE[LRUindex].LRU = maxLRU+1;
		if(TLBenabled){
			updateTLB(address,LRUindex);
		}
		return LRUindex;
	}

}

//Look in the TLB for a hit
bool TLBLookup(unsigned long address,unsigned long *tlbindex, unsigned long *tlbtag, unsigned long *paddress){
	unsigned long mask = 0xFFFFFFFF;
	
	unsigned long pagenum = (address >> pageOffset);
	*tlbindex = pagenum & (mask >> (32 - TLBIndexSize));
	*tlbtag = pagenum >> TLBIndexSize;
	int lookIndex = *tlbindex;
	int tag = *tlbtag;
	bool hit = false;
	int hitIndex = 0;
	int maxLRU = TLB[lookIndex][0].LRU;
	int minLRU = TLB[lookIndex][0].LRU;
	
	int maxLRUIndex = 0;
	int minLRUIndex = 0;

	//Search through the tlb
	for(int i=0; i < TLBSetSize; i++){ 
		if(TLB[lookIndex][i].tag == tag){
			//printf("Found hit in lookIndex %x index %x for tag %x\n",lookIndex,	i, tag);
			hit = true;
			hitIndex = i;
		}
		if(TLB[lookIndex][i].LRU > maxLRU){
			maxLRU = TLB[lookIndex][i].LRU;
			
			maxLRUIndex = i;
		}
		if(TLB[lookIndex][i].LRU < minLRU){
			minLRU = TLB[lookIndex][i].LRU;
			minLRUIndex = i;
		}
	}
	

	if(hit){//On hit we can use the address to look in the pagetable
		TLB[lookIndex][hitIndex].LRU = maxLRU+1;
		*paddress = TLB[lookIndex][hitIndex].address;
		unsigned long maxLRU = 0;
		for(int i = 0; i < numPPages; i++){
			
			if(PAGETABLE[i].LRU > maxLRU){
				maxLRU = PAGETABLE[i].LRU;
			}
		}
		
		PAGETABLE[TLB[lookIndex][hitIndex].address].LRU = maxLRU+1;

		dtlbhits++;
		return true;
	}else{
		dtlbmisses++;
		return false;
	}

}




//Gets information about a physical address
void physicalInfo(unsigned long *address, unsigned long *po, unsigned long *pagenum){
	unsigned long mask = 0xFFFFFFFF;
	*po = *address & (mask >> (32 - pageOffset));
	*pagenum = (*address >> pageOffset);
}



int main(int argc, char *argv[]){
	readConfig();
	if(virtAdd){
		printf("Virtual ");
	}else{
		printf("Physical");
	}
	
	printf(" Virt.  Page TLB    TLB TLB  PT   Phys        DC  DC          L2  L2\nAddress  Page # Off  Tag    Ind Res. Res. Pg # DC Tag Ind Res. L2 Tag Ind Res.\n-------- ------ ---- ------ --- ---- ---- ---- ------ --- ---- ------ --- ----\n");
	if(TLBenabled){
		makeTLB();
	}
	//printTLB();
	if(virtAdd){
		makePageTable();
	}
	makeL1Cache();
	if(l2){
			makeL2Cache();
	}
	//printL1Cache();
	stringstream ss;
	string temp;
	char RW;
	unsigned long mask = 0xFFFFFFFF;
	unsigned long address;
	unsigned long paddress;
	unsigned long po;
	unsigned long pagenum;
	unsigned long l1tag;
	unsigned long l1index;
	unsigned long *l2tag = NULL;
	unsigned long *l2index = NULL;
	unsigned long *tlbtag = NULL;
	unsigned long *tlbindex = NULL;
	bool *didPageTableHit = NULL;
	bool didL1Hit;
	bool *didL2Hit;
	bool *wasL2Accesses;
	bool *didTLBHit = NULL;
	unsigned long *vpnum = NULL;
	unsigned long pindex;

	while(cin >> temp){ //Read standard in until end
		
		RW = temp[0];
		if(temp.size() == 2){
			cin >> temp;
		}else{
			temp = temp.substr(2,temp.size());
		}
		

		if(RW=='R'){
			totalReads++;
		}else if(RW=='W'){
			totalWrites++;
		}else{
			cout << "BAD input" << endl;
		}

		bool L2Relevent = false;
		
		
		ss << hex << temp;
		ss >> address;
		ss = stringstream();

		if(virtAdd){//If virtual addresses are enabled then we need the page table and or the tlb to translate the address
			
			if(TLBenabled){
				didTLBHit = (bool *)malloc(sizeof(bool));
				tlbtag = (unsigned long*)malloc(sizeof(unsigned long));
				tlbindex = (unsigned long*)malloc(sizeof(unsigned long));
				*didTLBHit = TLBLookup(address,tlbindex,tlbtag,&pindex);
				paddress = ((pindex) << pageOffset) | (address & (mask >> (32 - pageOffset)));
				if(*didTLBHit){
					didPageTableHit = NULL;
				}else{
					didPageTableHit = (bool *)malloc(sizeof(bool));
					int index = pageTableLookup(address,didPageTableHit);
					address & (mask >> (32 - pageOffset));
					//printf("ppn: %x and offset: %x\n",((index) << pageOffset), (address & (mask >> (32 - pageOffset))));
					paddress = ((index) << pageOffset) | (address & (mask >> (32 - pageOffset)));
					
				}

			}else{
				didPageTableHit = (bool *)malloc(sizeof(bool));
				int index = pageTableLookup(address,didPageTableHit);
				address & (mask >> (32 - pageOffset));
				//printf("ppn: %x and offset: %x\n",((index) << pageOffset), (address & (mask >> (32 - pageOffset))));
				paddress = ((index) << pageOffset) | (address & (mask >> (32 - pageOffset)));
				
			}
			vpnum = (unsigned long*)malloc(sizeof(unsigned long));
			*vpnum = (address >> pageOffset);
		}else{
			paddress = address;
		}

		//Get info about the address
		physicalInfo(&paddress,&po,&pagenum);
		l1Info(&paddress,&l1tag,&l1index);
		if(l2){
			L2Relevent = true;
			l2tag = (unsigned long*)malloc(sizeof(unsigned long));
			l2index = (unsigned long*)malloc(sizeof(unsigned long));
			didL2Hit = (bool *)malloc(sizeof(bool));
			*didL2Hit = false;
			wasL2Accesses = (bool *)malloc(sizeof(bool));
			*wasL2Accesses = false;
			l2Info(&paddress,l2tag,l2index);
		}
		
		didL1Hit = accessL1Cache(paddress, RW, l1index, l1tag, didL2Hit, wasL2Accesses);
		if(!*wasL2Accesses){
			l2tag = NULL;
			l2index = NULL;
			didL2Hit = NULL;
			L2Relevent = false;
		}
		
		printLine(&address, vpnum, &po, tlbtag, tlbindex, didTLBHit, didPageTableHit, &pagenum,&l1tag,&l1index,didL1Hit, l2tag, l2index, didL2Hit, L2Relevent);
		//printL1Cache();
		//printL2Cache();
		//printTLB();
		//printPageTable();
		//break;
		
	}

	//Output statistics at the end of the simulation
	printStatistics();
	

}
