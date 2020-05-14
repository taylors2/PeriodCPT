#ifndef FILE_CACHE
#define FILE_CACHE

#include "PeriodCPT_General.h"
#include "PeriodCPT_MCMCgeneric.h"

cache_t *Make_Cache_Item(MCMCitem_t *generator, chain_t *chain){
  //Rprintf("Make_Cache_Item\n");
  cache_t *p   = (cache_t *)my_calloc(1, sizeof(cache_t));
  p->generator = generator;
  p->chain     = chain;
  p->count     = 5;  //initialise
  return p;
}

void Delete_Cache_Item(cache_t *item){
  //Rprintf("Delete_Cache_Item\n");
  if(item != NULL){
    Delete_MCMCitem(item->generator);
    Delete_Chain(item->chain);
    my_free(item);
  }
  return;
}

cache_t **Make_Cache_List(int maxlen){
  //Rprintf("Make_Cache_List\n");
  cache_t **p = (cache_t **)my_calloc(maxlen,sizeof(cache_t *));
  return p;
}

void Delete_Cache_List(cache_t **list, int maxlen){
  //Rprintf("Delete_Cache_List\n");
  for(int i = 0; i < maxlen; i++){
    Delete_Cache_Item(list[i]);
    list[i] = NULL;
  }
  my_free(list);
}


int Compare_Cache_Item(cache_t *test, MCMCitem_t *instance){
  //Rprintf("Compare_Cache_Item\n");
  return Compare_MCMCitem(test->generator, instance, TRUE);
}

cache_t *Find_in_Cache(cache_t **cache, int *n, MCMCitem_t *mcmc){
  //Rprintf("Find_in_Cache\n");
  int i = 0;
  int compare = 0;
  while((i < *n) & (compare == 0)){
    compare = Compare_Cache_Item(cache[i], mcmc);
    if(compare == 0) i++;
  }
  if(compare == 1){
    //cache[i]->count++;
    int done = 0;
    while((i>0) & (done != 1)){
      if(cache[i]->count > cache[i-1]->count){
        cache_t *tmp;
        tmp        = cache[i];
        cache[i]   = cache[i-1];
        cache[i-1] = tmp;
        i--;
      }else{
        done = 1;
      }
    }
    return cache[i];
  }else{
    return NULL;
  }
}

void Push_To_Cache(cache_t **cache, int *last, int *length, cache_t *item){
  //Rprintf("Push_To_Cache\n");
  //Assume item does not already exist in cache
  int n = *last;
  if(n == *length){
    Delete_Cache_Item(cache[n-1]);
    cache[n-1] = NULL;
    n--;
  }

  cache[n] = item;
  n++;
  *last = n;

  //Bubble up
  int i = n-1;
  int done = 0;
  while((i>0) & (done != 1)){
    if(cache[i]->count > cache[i-1]->count){
      cache_t *tmp;
      tmp        = cache[i];
      cache[i]   = cache[i-1];
      cache[i-1] = tmp;
      i--;
    }else{
      done = 1;
    }
  }


  return;
}


void Manage_Cache(cache_t **cache, int *ncache, int *maxcache){
  //Rprintf("Manage_Cache\n");
  //discount all counts
  for(int i = 0; i<*ncache; i++){
    cache[i]->count--;
  }
  int j = *ncache - 1;
  while(cache[j]->count <= 0){
    Delete_Cache_Item(cache[j]);
    cache[j] = NULL;
    *ncache = j;
    j--;
    if(*ncache == 0){
      return;
    }
  }
  return;
}


void Bubble_Chace_Item(cache_t **cache, int *n, int *max, cache_t *psetj){
  //Rprintf("Bubble_Chace_Item\n");
  int i = 0;
  while(cache[i] != psetj) i++;

  //Bubble up
  int done = 0;
  while((i>0) & (done != 1)){
    if(cache[i]->count > cache[i-1]->count){
      cache_t *tmp;
      tmp        = cache[i];
      cache[i]   = cache[i-1];
      cache[i-1] = tmp;
      i--;
    }else{
      done = 1;
    }
  }
  return;
}



#endif //FILE_CACHE
