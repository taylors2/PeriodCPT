/* 

#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdio.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdlib.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\math.h>
#include <C:\Program Files\R\R-3.4.3\include\R.h>
#include <C:\Program Files\R\R-3.4.3\include\Rmath.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\time.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\stdint.h>
#include <C:\Rtools\gcc-4.6.3\i686-w64-mingw32\include\string.h>
#include"CircCPT_General.c"
#include"CircCPT_MCMCcache.c"


 */



list_t * makeList(int maxlen){
  int fnname = 400;
  profile(TRUE, fnname);
  list_t *p = (list_t *)my_calloc(1,sizeof(list_t));
  p->first = NULL;
  p->last = NULL;
  p->length = 0;
  p->maxlen = maxlen;
  profile(FALSE, fnname);
  return p;
}

Prop_t * makePropItem(MCMCitem_t *mcmc, chain_t *chain){
  int fnname = 401;
  profile(TRUE, fnname);
  Prop_t *this;

  this = (Prop_t *)my_calloc(1,sizeof(Prop_t));
  this->m = mcmc->m;
  this->j = mcmc->j;
  this->tau = (int*)my_calloc(this->m,sizeof(int));
  int i;
  for(i=0; i<this->m; i++) this->tau[i] = mcmc->tau[i];
  this->value = mcmc->value;
  this->chain = chain;
  this->prev = NULL;
  this->next = NULL;
  profile(FALSE, fnname);
  return this;
}

void DeletePropItem(Prop_t *this){
  int fnname = 402;
  profile(TRUE, fnname);
  my_free(this->tau);
  DeleteMCMCchain(this->chain);
  my_free(this);
  profile(FALSE, fnname);
  return;
}

void PopLastListItem(list_t *list){
  int fnname = 403;
  profile(TRUE, fnname);
  Prop_t *this;
  this = list->last;
  if(this!=NULL){
    list->last = this->prev;
    list->length--;
    if(list->last == NULL){
      list->first = NULL;
    }else{
      list->last->next = NULL;
    }
    DeletePropItem(this);
  }
  profile(FALSE, fnname);
  return;
}

void DeleteList(list_t *list){
  int fnname = 404;
  profile(TRUE, fnname);
  Prop_t *this, *next;
  this = list->first;
  while(this!=NULL){
    next = this->next;
    DeletePropItem(this);
    this = next;
  }
  //while(list->first != NULL) PopLastListItem(list);
  my_free(list);
  profile(FALSE, fnname);
  return;
}

void PushToList(list_t *list, MCMCitem_t *current, chain_t *chain){
  int fnname = 405;
  profile(TRUE, fnname);

  if(list->length == list->maxlen){
    PopLastListItem(list);  //!!!should not remove last item if it contains current!!!
  }
  
  Prop_t *this,*place;
  int point = FALSE;
  this = makePropItem(current, chain);
  if(list->first == NULL){
    list->first = this;
    list->last = this;
  }else{
    place = list->first;
    while(point == FALSE){
      if(place==NULL){
        point = TRUE;
      }else if(place->value <= this->value){
        point = TRUE;
      }else{
        place = place->next;
      }
    }
    //put this before place
    if(place==NULL){ //push to the end
      this->prev = list->last;
      this->next = NULL;
      this->prev->next = this;
      list->last = this;
    }else if(place->prev == NULL){ //push to the front
      place->prev = this;
      this->next = place;
      this->prev = NULL;
      list->first = this;
    }else{
      place->prev->next = this;
      this->prev = place->prev;
      place->prev = this;
      this->next = place;
    }
  }
  list->length++;

  profile(FALSE, fnname);
  return;
}

int CompareCase(MCMCitem_t *current, Prop_t *this){
  int fnname = 406;
  profile(TRUE, fnname);

  if(this->value != current->value){
    profile(FALSE, fnname);
    return FALSE;
  }
  if(this->m != current->m){
    profile(FALSE, fnname);
    return FALSE;
  }
  if(this->j != current->j){
    profile(FALSE, fnname);
    return FALSE;
  }
  int i;
  for(i=0;i<this->m;i++){
    if(current->tau[i] != this->tau[i]){
      profile(FALSE, fnname);
      return FALSE;
    }
  }
  profile(FALSE, fnname);
  return TRUE;
}

chain_t *FindCase(list_t *list, MCMCitem_t *current){
  int fnname = 407;
  profile(TRUE, fnname);
  
  Prop_t *this;
  this = list->first;
  int search;

  while(this != NULL){
    search = CompareCase(current, this);
    if(search==TRUE){
      profile(FALSE, fnname);
      return this->chain;
    }
    this = this->next;
  }
  //If here then not found case in list
  profile(FALSE, fnname);
  return NULL;
}

void PrintPropCache(list_t * list, int index){
  int fnname = 408;
  profile(TRUE, fnname);
  
  Prop_t *this;
  this = list->first;
  int maxm;
  maxm = 0;
  while(this != NULL){
    if(maxm < this->m) maxm = this->m;
    this = this->next;
  }
  
  FILE *f;
  int j;
  char fname[100]; 
  sprintf(fname,"%d_PROPlist.txt", index);

  f = fopen(fname, "w");
  fprintf(f,"me, prev, next, id, m, j, value, chain");
  for(j=0; j<maxm; j++) fprintf(f,", tau%d", j+1);
  fprintf(f,"\n");
  
  this = list->first;
  int i = 0;
  while(i < list->length){
    fprintf(f,"%p, ", this);
    if(this->prev == NULL){
      fprintf(f,"-1, ");
    }else{
      fprintf(f,"%p, ",this->prev);
    }
    if(this->next == NULL){
      fprintf(f,"-1, ");
    }else{
      fprintf(f,"%p, ",this->next);
    }
    fprintf(f,"%d, %d, %d, %f, %p", i, this->m, this->j, this->value, this->chain);
    for(j=0;j<maxm;j++){
      if(j<this->m){
        fprintf(f,", %d", this->tau[j]);
      }else{
        fprintf(f,", -1");
      }
    }
    fprintf(f,"\n");
    this = this->next;
    i++;
  }
  fclose(f);
  
  profile(FALSE, fnname);
  return;
}






