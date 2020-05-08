#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

typedef struct segval {
  int cpt;
  int seglen;
  double value;
} segval_t;
