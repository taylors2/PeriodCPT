#include "PeriodCPT_General.h"
#include "PeriodCPT_MCMCgeneric.h"

cache_t *Make_Cache_Item(MCMCitem_t *generator, chain_t *chain){
  cache_t *p   = (cache_t *)calloc(1, sizeof(cache_t));
  p->generator = generator;
  p->chain     = chain;
  p->count     = 0;
  return p;
}

void Delete_Cache_Item(cache_t *item){
  if(item != 0){
    Delete_MCMCitem(item->generator);
    Delete_Chain(item->chain);
  }
  return;
}

cache_t **Make_Cache_List(int maxlen){
  cache_t **p = (cache_t **)calloc(maxlen,sizeof(cache_t *));
  return p;
}

void Delete_Cache_List(cache_t **list, int maxlen){
  for(int i = 0; i < maxlen; i++){
    Delete_Cache_Item(list[i]);
    list[i] = 0;
  }
  free(list);
}


int Compare_Cache_Item(cache_t *test, MCMCitem_t *instance){
  return Compare_MCMCitem(test->generator, instance, TRUE);
}

chain_t *Find_in_Cache(cache_t **cache, int *n, MCMCitem_t *mcmc){
  int compare, i;
  compare = 0;
  for(int i = 0; (i < *n) && (compare == 0); i++){
    compare = Compare_Cache_Item(cache[i], mcmc);
  }
  if(compare == 1){
    cache[i]->count++;
    while((i>0) && (cache[i]->count > cache[i-1]->count)){
      cache_t *tmp;
      tmp        = cache[i-1];
      cache[i-1] = cache[i];
      cache[i]   = tmp;
      i--;
    }
    return cache[i]->chain;
  }else{
    return NULL;
  }
}

void Push_To_Cache(cache_t **cache, int *last, int *length, cache_t *item){
  //Assume item does not already exist in cache
  int n = *last;
  if(n == *length){
    Delete_Cache_Item(cache[*last-1]);
    n--;
  }
  cache[n] = item;
  n++;
  *last = n;
  return;
}


