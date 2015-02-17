#include "gfcm.h"
#include "mem.h"
#include <math.h>

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t CalcMult(uint32_t c, uint32_t s){
  uint32_t n;
  uint64_t x[c], p = 1;
  for(n = 0 ; n < c ; ++n){ x[n] = p; p *= s; }
  return x[c-1];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitGArrayTab(GFCM *M){
  M->A.cnts = (GACC *) Calloc(M->nPMod*M->nSym, sizeof(GACC));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FreeGFCM(GFCM *M){
  Free(M->A.cnts);
  Free(M->freqs);
  Free(M);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

GFCM *CreateGFCM(uint32_t c, uint32_t a, uint8_t n){
  GFCM   *M = (GFCM *) Calloc(1, sizeof(GFCM));
  M->nSym   = n;
  M->ctx    = c;
  M->aDen   = a;
  M->mult   = CalcMult(c, n);
  M->idx    = 0;
  M->freqs  = (uint32_t *) Calloc(M->nSym+1, sizeof(uint32_t));
  M->nPMod  = (uint64_t) pow(n,c);

  InitGArrayTab(M);

  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateGFCM(GFCM *M, uint32_t c){
  uint32_t n;
  GACC *ac = &M->A.cnts[M->idx*M->nSym];
  if(++ac[c] == MAXGACC_C){
    for(n = 0 ; n < M->nSym ; ++n)
      ac[n]>>=1;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputeGFCM(GFCM *M){
  uint32_t x;
  GACC *a = &M->A.cnts[M->idx*M->nSym];

  for(x = 0 ; x < M->nSym ; ++x)
    M->freqs[x] = 1 + M->aDen * a[x];

  M->freqs[M->nSym] = M->freqs[0];
  for(x = 1 ; x < M->nSym ; ++x)
    M->freqs[M->nSym] += M->freqs[x];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void GetIdx(uint8_t *p, GFCM *M){
  M->idx = ((M->idx-*(p-M->ctx)*M->mult)*M->nSym)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


