#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "defs.h"
#include "mem.h"
#include "common.h"
#include "context.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static HCCounters zeroCounters = {0x00, 0x00, 0x00, 0x00};
static HCCounters auxCounters;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t XHASH(uint64_t x){
  return (x * 786433 + 196613) % 68719476735;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitHashTable(CModel *M, U32 c){ 
  uint32_t k;
  M->hTable.maxC    = c;
  M->hTable.index   = (ENTMAX *) Calloc(HASH_SIZE, sizeof(ENTMAX));
  M->hTable.entries = (Entry **) Calloc(HASH_SIZE, sizeof(Entry *));
  for(k = 0 ; k < HASH_SIZE ; ++k)
    M->hTable.entries[k] = (Entry *) Calloc(M->hTable.maxC, sizeof(Entry));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FreeCModel(CModel *M){
  U32 k;
  if(M->mode == HASH_TABLE_MODE){
    for(k = 0 ; k < HASH_SIZE ; ++k)
      Free(M->hTable.entries[k]);
    Free(M->hTable.entries);
    Free(M->hTable.index);
    }
  else // TABLE_MODE
    Free(M->array.counters);
  Free(M);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitArray(CModel *M){
  M->array.counters = (ACC *) Calloc(M->nPModels<<2, sizeof(ACC));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InsertKey(HashTable *H, U32 hi, U64 idx, U8 s){
  if(++H->index[hi] == H->maxC)
    H->index[hi] = 0;

  #ifdef PREC32B
  H->entries[hi][H->index[hi]].key = (U32)(idx&0xffffffff);
  #else
  H->entries[hi][H->index[hi]].key = (U16)(idx&0xffff);
  #endif  
  H->entries[hi][H->index[hi]].counters = (0x01<<(s<<1));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static HCC *GetHCCounters(HashTable *H, U64 key){
  U32 n, hIndex = key % HASH_SIZE;
  #ifdef PREC32B
  U32 b = key & 0xffffffff;
  #else
  U16 b = key & 0xffff;
  #endif

  #ifdef FSEARCHMODE
  U32 pos = H->index[hIndex];
  // FROM INDEX-1 TO 0
  for(n = pos+1 ; n-- ; ){
    if(H->entries[hIndex][n].key == b){
      auxCounters[0] =  H->entries[hIndex][n].counters &  0x03;
      auxCounters[1] = (H->entries[hIndex][n].counters & (0x03<<2))>>2;
      auxCounters[2] = (H->entries[hIndex][n].counters & (0x03<<4))>>4;
      auxCounters[3] = (H->entries[hIndex][n].counters & (0x03<<6))>>6;
      return auxCounters;
      }
    }

  // FROM MAX_COLISIONS TO INDEX
  for(n = (H->maxC-1) ; n > pos ; --n){
    if(H->entries[hIndex][n].key == b){
      auxCounters[0] =  H->entries[hIndex][n].counters &  0x03;
      auxCounters[1] = (H->entries[hIndex][n].counters & (0x03<<2))>>2;
      auxCounters[2] = (H->entries[hIndex][n].counters & (0x03<<4))>>4;
      auxCounters[3] = (H->entries[hIndex][n].counters & (0x03<<6))>>6;
      return auxCounters;
      }
    }
  #else
  // FROM 0 TO MAX
  for(n = 0 ; n < H->maxC ; ++n){
    if(H->entries[hIndex][n].key == b){
      auxCounters[0] =  H->entries[hIndex][n].counters &  0x03;
      auxCounters[1] = (H->entries[hIndex][n].counters & (0x03<<2))>>2;
      auxCounters[2] = (H->entries[hIndex][n].counters & (0x03<<4))>>4;
      auxCounters[3] = (H->entries[hIndex][n].counters & (0x03<<6))>>6;
      return auxCounters;
      }
    }
  #endif

  return NULL;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PModel *CreatePModel(U32 n){
  PModel *P = (PModel *) Malloc(sizeof(PModel));
  P->freqs  = (U32    *) Malloc(n* sizeof(U32));
  return P;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FloatPModel *CreateFloatPModel(U32 n){
  FloatPModel *F = (FloatPModel *) Malloc(sizeof(FloatPModel));
  F->freqs = (double *) Malloc(n * sizeof(double));
  return F;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef SWAP
void SwapPos(Entry *A, Entry *B){
  Entry *tmp = B;
  *B = *A;
  A = tmp;
  }
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateCModelCounter(CModel *M, U32 sym, U64 im){
  U32 n;
  ACC *AC;
  U64 idx = im;

  if(M->mode == HASH_TABLE_MODE){
    U8   counter;
    U32  s, hIndex = (idx = XHASH(idx)) % HASH_SIZE;
    #ifdef PREC32B
    U32 b = idx & 0xffffffff;
    #else
    U16 b = idx & 0xffff;
    #endif

    for(n = 0 ; n < M->hTable.maxC ; n++){
      if(M->hTable.entries[hIndex][n].key == b){
        counter = (M->hTable.entries[hIndex][n].counters>>(sym<<1))&0x03;
        if(counter == 3){
          for(s = 0 ; s < 4 ; ++s){
            if(s != sym){
              counter = 
              ((M->hTable.entries[hIndex][n].counters>>(s<<1))&0x03)>>1;
              M->hTable.entries[hIndex][n].counters &= ~(0x03<<(s<<1));
              M->hTable.entries[hIndex][n].counters |= (counter<<(s<<1));
              }
            }
          #ifdef SWAP
          SwapPos(&M->hTable.entries[hIndex][n],
          &M->hTable.entries[hIndex][M->hTable.index[hIndex]]);
          #endif
          return;
          }
        else{ // THERE IS STILL SPACE FOR INCREMENT COUNTER
          ++counter;
          M->hTable.entries[hIndex][n].counters &= ~(0x03<<(sym<<1));
          M->hTable.entries[hIndex][n].counters |= (counter<<(sym<<1));
          #ifdef SWAP
          SwapPos(&M->hTable.entries[hIndex][n],
          &M->hTable.entries[hIndex][M->hTable.index[hIndex]]);
          #endif
          return;
          }
        }
      }
    InsertKey(&M->hTable, hIndex, b, sym); // KEY NOT FOUND: WRITE ON OLDEST
    }
  else{
    AC = &M->array.counters[idx << 2];
    if(++AC[sym] == M->maxCount){    
      AC[0] >>= 1;
      AC[1] >>= 1;
      AC[2] >>= 1;
      AC[3] >>= 1;
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CModel *CreateCModel(U32 ctx, U32 aDen, U32 ir, U8 ref, U32 col, U32 am){
  CModel *M = (CModel *) Calloc(1, sizeof(CModel));
  U64    prod = 1, *mult;
  U32    n;

  if(ctx > MAX_HASH_CTX){
    fprintf(stderr, "Error: context size cannot be greater than %d\n", 
    MAX_HASH_CTX);
    exit(1);
    }

  mult           = (U64 *) Calloc(ctx, sizeof(U64));
  M->nPModels    = (U64) pow(ALPHABET_SIZE, ctx);
  M->ctx         = ctx;
  M->alphaDen    = aDen;
  M->am          = am;
  M->pModelIdx   = 0;
  M->pModelIdxIR = M->nPModels - 1;
  M->ir          = ir  == 0 ? 0 : 1;
  M->ref         = ref == 0 ? 0 : 1;

  if(ctx >= HASH_TABLE_BEGIN_CTX){
    M->mode     = HASH_TABLE_MODE;
    M->maxCount = DEFAULT_MAX_COUNT >> 8;
    InitHashTable(M, col);
    }
  else{
    M->mode     = ARRAY_MODE;
    M->maxCount = DEFAULT_MAX_COUNT;
    InitArray(M);
    }

  for(n = 0 ; n < M->ctx ; ++n){
    mult[n] = prod;
    prod <<= 2;
    }

  M->multiplier = mult[M->ctx-1];

  #ifdef CORRECT
  if(am != 0){
    M->correct.seq       = CreateCBuffer(BUFFER_SIZE, BGUARD);
    M->correct.in        = 0;
    M->correct.idx       = 0;
    M->correct.idxRev    = M->nPModels-1;
    M->correct.mask      = (uint8_t *) Calloc(BGUARD, sizeof(uint8_t));
    M->correct.threshold = am;
    }
  #endif

  Free(mult);
  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef CORRECT
int32_t BestId(uint32_t *f, uint32_t sum){
  if(sum == 4) return -2; // ZERO COUNTERS

  uint32_t x, best = 0, max = f[0];
  for(x = 1 ; x < 4 ; ++x)
    if(f[x] > max){
      max = f[x];
      best = x;
      }

  for(x = 0 ; x < 4 ; ++x) if(best != x && max == f[x]) return -1;

  return best;
  }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ResetCModelIdx(CModel *M){
  M->pModelIdx   = 0;
  M->pModelIdxIR = M->nPModels - 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline U8 GetPModelIdxIR(U8 *p, CModel *M){
  M->pModelIdxIR = (M->pModelIdxIR>>2)+GetCompNum(*p)*M->multiplier;
  return GetCompNum(*(p-M->ctx));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline U8 GetPModelIdxIRCorr(U8 *p, CModel *M){
  M->correct.idxRev = (M->correct.idxRev>>2)+GetCompNum(*p)*M->multiplier;
  return GetCompNum(*(p-M->ctx));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void GetPModelIdx(U8 *p, CModel *M){
  M->pModelIdx = ((M->pModelIdx-*(p-M->ctx)*M->multiplier)<<2)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void GetPModelIdxCorr(U8 *p, CModel *M){
  M->correct.idx = ((M->correct.idx-*(p-M->ctx)*M->multiplier)<<2)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef CORRECT
void Fail(CModel *M){
  uint32_t x, fails = 0;
  for(x = 0 ; x < M->ctx ; ++x)
    if(M->correct.mask[x] != 0)
      ++fails;
  if(fails <= M->correct.threshold)
    ShiftBuffer(M->correct.mask, M->ctx, 1);
  else 
    M->correct.in = 0;
  }
#endif
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef CORRECT
void Hit(CModel *M){
  ShiftBuffer(M->correct.mask, M->ctx, 0);
  }
#endif

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputePModel(CModel *M, PModel *P){
  HCC *hCounters;
  ACC *aCounters;

  if(M->mode == HASH_TABLE_MODE){
    if(!(hCounters = GetHCCounters(&M->hTable, XHASH(M->pModelIdx))))
      hCounters = zeroCounters;
    P->freqs[0] = 1 + M->alphaDen * hCounters[0];
    P->freqs[1] = 1 + M->alphaDen * hCounters[1];
    P->freqs[2] = 1 + M->alphaDen * hCounters[2];
    P->freqs[3] = 1 + M->alphaDen * hCounters[3];
    }
  else{
    aCounters = &M->array.counters[M->pModelIdx << 2];
    P->freqs[0] = 1 + M->alphaDen * aCounters[0];
    P->freqs[1] = 1 + M->alphaDen * aCounters[1];
    P->freqs[2] = 1 + M->alphaDen * aCounters[2];
    P->freqs[3] = 1 + M->alphaDen * aCounters[3];
    }

  P->sum = P->freqs[0] + P->freqs[1] + P->freqs[2] + P->freqs[3];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double PModelSymbolNats(PModel *P, U32 s){
  return log((double) P->sum / P->freqs[s]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
