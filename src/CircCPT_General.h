/* 

#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdio.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdlib.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\math.h>
#include <C:\Program Files\R\R-3.4.3\include\R.h>
#include <C:\Program Files\R\R-3.4.3\include\Rmath.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\time.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdint.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\string.h>

 */

#define TRUE 1
#define FALSE 0

#define DoMEM 0
int memoryID = 0;
char profFile[] = "Profile.txt";
int MEMMAX = 1000000000;
void **MEMaddress;
int *MEMalloc;

#define DoPROFILE 0
char ptFile[] = "pointerInformation.txt";
int ProfID = 0;
int *INOUT;
clock_t *DURATION;
int *FNID;
int PROFMAX = 10000000;

int INDEX = 0; //file index for printing cpt cache & mcmc chain

typedef struct MCMCitem{
  int *tau;
  int j;
  int m;
  double value;
  double prob;
  struct MCMCitem *prev;
  struct MCMCitem *next;
} MCMCitem_t;


typedef struct MCMCchain{
  MCMCitem_t *first;
  MCMCitem_t *last;
  int length;
} chain_t;


typedef struct PropItem{
  int m;
  int *tau;
  int j;
  double value;
  struct MCMCchain *chain;
  struct PropItem *prev;
  struct PropItem *next;
} Prop_t;

typedef struct PropList{
  int maxlen;
  int length;
  struct PropItem *first;
  struct PropItem *last;
} list_t;

struct G2
{
  int cpt;
  int seglen;
  double value;
};


void *my_calloc(size_t count, size_t size){ 
  void *pt;
  pt = calloc(count, size);
  if(DoMEM == TRUE){
    if(memoryID == 0){
      MEMaddress = (void **)calloc(MEMMAX, sizeof(void*));
      MEMalloc = (int*)calloc(MEMMAX, sizeof(int));
    }
    if(memoryID<MEMMAX){
      MEMalloc[memoryID] = TRUE;
      MEMaddress[memoryID] = pt;
      //Rprintf("SET:  %d, %p\n",MEMalloc[memoryID],MEMaddress[memoryID]);
    }
    memoryID++;
    if(memoryID==MEMMAX){
      FILE *file;
      file = fopen(ptFile, "w");
      int i;
      fprintf(file, "Allocate, Address\n");
      for(i=0; i<memoryID; i++) fprintf(file, "%d, %p\n",MEMalloc[i],MEMaddress[i]);
      fclose(file);
      free(MEMalloc);
      free(MEMaddress);
    }
  }
  return pt;
}

void my_free(void *pt){
  if(DoMEM == TRUE){
    if(memoryID<MEMMAX){
      MEMalloc[memoryID] = FALSE;
      MEMaddress[memoryID] = pt;
      //Rprintf("FREE: %d, %p\n",MEMalloc[memoryID],MEMaddress[memoryID]);
    }
    memoryID++;
    if((memoryID==MEMMAX) | ((pt==MEMaddress[0])&(memoryID<=MEMMAX))){  //nb first alloc MUST be last free (CACHEMAX!)
      FILE *file;
      file = fopen(ptFile, "w");
      int i;
      fprintf(file, "Allocate, Address\n");
      for(i=0; i<memoryID; i++) fprintf(file, "%d, %p\n",MEMalloc[i],MEMaddress[i]);
      fclose(file);
      free(MEMalloc);
      free(MEMaddress);
    }
  }
  free(pt);
  return;
}



void profile(int in, int fn){
  if((DoPROFILE != TRUE) | (ProfID>=PROFMAX)) return;
//  Rprintf("%d(%d)\n",fn,in);
  
  if((fn==0) & (in==TRUE)){
    //Started function, allocate memory for profiling
    INOUT = (int*)calloc(PROFMAX,sizeof(int));
    DURATION = (clock_t*)calloc(PROFMAX,sizeof(clock_t));
    FNID = (int*)calloc(PROFMAX,sizeof(int));
  }

  FNID[ProfID] = fn;
  INOUT[ProfID] = in;
  DURATION[ProfID] = clock();
  ProfID++;

  if(((fn==0) & (in==FALSE)) | (ProfID == PROFMAX)){
    //Finished function or hit cache max
    FILE *file;
    file = fopen(profFile, "w");
    fprintf(file, "ID, Time, In\n");
    int i;
    for(i=0;i<ProfID;i++) fprintf(file,"%d, %ld, %d\n",FNID[i], DURATION[i], INOUT[i]);
    fclose(file);
    free(INOUT);free(DURATION);free(FNID);
  }

  return;
}


