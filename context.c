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

static uint64_t XHASH(uint64_t x){
  return (x * 786433 + 196613) % 68719476735;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t YHASH(uint64_t y){
  return (y * 786491 + 216617) % 66719476787;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static uint64_t ZHASH(uint64_t z){
  z = (~z) + (z << 21);
  z = z    ^ (z >> 24);
  z = (z   + (z << 3)) + (z << 8);
  z = z    ^ (z >> 14);
  z = (z   + (z << 2)) + (z << 4);
  z = z    ^ (z >> 28);
  z = z    + (z << 31);
  return z;
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

  #if defined(PREC32B)
  H->entries[hi][H->index[hi]].key = (U32)(idx&0xffffffff);
  #elif defined(PREC16B)
  H->entries[hi][H->index[hi]].key = (U16)(idx&0xffff);
  #else
  H->entries[hi][H->index[hi]].key = (U8)(idx&0xff);
  #endif  
  H->entries[hi][H->index[hi]].counters = (0x01<<(s<<2));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void GetFreqsFromHCC(HCC c, uint32_t a, PModel *P){
   P->sum  = (P->freqs[0] = a * ( c &  0x0f) + 1);
   P->sum += (P->freqs[1] = a * ((c & (0x0f<<4))>>4) + 1);
   P->sum += (P->freqs[2] = a * ((c & (0x0f<<8))>>8) + 1);
   P->sum += (P->freqs[3] = a * ((c & (0x0f<<12))>>12) + 1);
   }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void GetHCCounters(HashTable *H, U64 key, PModel *P, uint32_t a){
  U32 n, hIndex = key % HASH_SIZE;
  #if defined(PREC32B)
  U32 b = key & 0xffffffff;
  #elif defined(PREC16B)
  U16 b = key & 0xffff;
  #else
  U8  b = key & 0xff;
  #endif

  #ifdef FSEARCHMODE
  U32 pos = H->index[hIndex];
  // FROM INDEX-1 TO 0
  for(n = pos+1 ; n-- ; ){
    if(H->entries[hIndex][n].key == b){
      GetFreqsFromHCC(H->entries[hIndex][n].counters, a, P);
      return;
      }
    }
  // FROM MAX_COLISIONS TO INDEX
  for(n = (H->maxC-1) ; n > pos ; --n){
    if(H->entries[hIndex][n].key == b){
      GetFreqsFromHCC(H->entries[hIndex][n].counters, a, P);
      return;
      }
    }
  #else
  // FROM 0 TO MAX
  for(n = 0 ; n < H->maxC ; ++n){
    if(H->entries[hIndex][n].key == b){
      GetFreqsFromHCC(H->entries[hIndex][n].counters, a, P);
      return;
      }
    }
  #endif

  P->freqs[0] = 1;
  P->freqs[1] = 1;
  P->freqs[2] = 1;
  P->freqs[3] = 1;
  P->sum      = 4; 
  return;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PModel *CreatePModel(U32 n){
  PModel *P = (PModel *) Calloc(1, sizeof(PModel));
  P->freqs  = (U32    *) Calloc(n, sizeof(U32));
  return P;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FloatPModel *CreateFloatPModel(U32 n){
  FloatPModel *F = (FloatPModel *) Calloc(1, sizeof(FloatPModel));
  F->freqs = (double *) Calloc(n, sizeof(double));
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
    U16 counter;
    U32 s, hIndex = (idx = ZHASH(idx)) % HASH_SIZE;
    #if defined(PREC32B)
    U32 b = idx & 0xffffffff;
    #elif defined(PREC16B)
    U16 b = idx & 0xffff;
    #else
    U8  b = idx & 0xff;
    #endif

    for(n = 0 ; n < M->hTable.maxC ; n++){
      if(M->hTable.entries[hIndex][n].key == b){
        counter = (M->hTable.entries[hIndex][n].counters>>(sym<<2))&0x0f;
        if(counter == 15){
          for(s = 0 ; s < 4 ; ++s){
            counter = ((M->hTable.entries[hIndex][n].counters>>(s<<2))&0x0f)>>1;
            M->hTable.entries[hIndex][n].counters &= ~(0x0f<<(s<<2));
            M->hTable.entries[hIndex][n].counters |= (counter<<(s<<2));
            }

          counter = (M->hTable.entries[hIndex][n].counters>>(sym<<2))&0x0f;
          ++counter;
          M->hTable.entries[hIndex][n].counters &= ~(0x0f<<(sym<<2));
          M->hTable.entries[hIndex][n].counters |= (counter<<(sym<<2));
          return;
          }
        else{ // THERE IS STILL SPACE FOR INCREMENT COUNTER
          ++counter;
          M->hTable.entries[hIndex][n].counters &= ~(0x0f<<(sym<<2));
          M->hTable.entries[hIndex][n].counters |= (counter<<(sym<<2));
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

CModel *CreateCModel(U32 ctx, U32 aDen, U32 ir, U8 ref, U32 col, U32 edits, 
U32 eDen){
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
  M->edits       = edits;
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

  if(edits != 0){
    // SUBSTITUTIONS
    M->SUBS.seq       = CreateCBuffer(BUFFER_SIZE, BGUARD);
    M->SUBS.in        = 0;
    M->SUBS.idx       = 0;
    M->SUBS.mask      = (uint8_t *) Calloc(BGUARD, sizeof(uint8_t));
    M->SUBS.threshold = edits;
    M->SUBS.eDen      = eDen;
    // ADDITIONS
    M->ADDS.seq       = CreateCBuffer(BUFFER_SIZE, BGUARD);
    M->ADDS.in        = 0;
    M->ADDS.idx       = 0;
    M->ADDS.idx2      = 0;
    M->ADDS.mask      = (uint8_t *) Calloc(BGUARD, sizeof(uint8_t));
    M->ADDS.threshold = edits;
    M->ADDS.eDen      = eDen;
    // DELETIONS
    M->DELS.seq       = CreateCBuffer(BUFFER_SIZE, BGUARD);
    M->DELS.in        = 0;
    M->DELS.idx       = 0;
    M->DELS.idx2      = 0;
    M->DELS.mask      = (uint8_t *) Calloc(BGUARD, sizeof(uint8_t));
    M->DELS.threshold = edits;
    M->DELS.eDen      = eDen;
    }

  Free(mult);
  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t BestId2(uint32_t *f, uint32_t sum){

  uint32_t x, best = 0, max = f[0];
  for(x = 1 ; x < 4 ; ++x)
    if(f[x] > max){
      max = f[x];
      best = x;
      }

  return best;
  }

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

inline void GetPModelIdx(U8 *p, CModel *M){
  M->pModelIdx = ((M->pModelIdx-*(p-M->ctx)*M->multiplier)<<2)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline uint64_t GetPModelIdxCorr(U8 *p, CModel *M, uint64_t idx){
  return (((idx-*(p-M->ctx)*M->multiplier)<<2)+*p);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FailSUBS(CModel *M){
  uint32_t x, fails = 0;
  for(x = 0 ; x < M->ctx ; ++x)
    if(M->SUBS.mask[x] != 0)
      ++fails;
  if(fails <= M->SUBS.threshold)
    ShiftBuffer(M->SUBS.mask, M->ctx, 1);
  else 
    M->SUBS.in = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FailADDS(CModel *M){
  uint32_t x, fails = 0;
  for(x = 0 ; x < M->ctx ; ++x)
    if(M->ADDS.mask[x] != 0)
      ++fails;
  if(fails <= M->ADDS.threshold)
    ShiftBuffer(M->ADDS.mask, M->ctx, 1);
  else
    M->ADDS.in = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FailDELS(CModel *M){
  uint32_t x, fails = 0;
  for(x = 0 ; x < M->ctx ; ++x)
    if(M->DELS.mask[x] != 0)
      ++fails;
  if(fails <= M->DELS.threshold)
    ShiftBuffer(M->DELS.mask, M->ctx, 1);
  else
    M->DELS.in = 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void HitSUBS(CModel *M){
  ShiftBuffer(M->SUBS.mask, M->ctx, 0);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void HitADDS(CModel *M){
  ShiftBuffer(M->ADDS.mask, M->ctx, 0);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void HitDELS(CModel *M){
  ShiftBuffer(M->DELS.mask, M->ctx, 0);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CorrectCModelSUBS(CModel *M, PModel *P, uint8_t sym){
  int32_t best = BestId(P->freqs, P->sum);
  switch(best){
    case -2:  // IT IS A ZERO COUNTER [NOT SEEN BEFORE]
      if(M->SUBS.in != 0)
        FailSUBS(M);
    break;
    case -1:  // IT HAS AT LEAST TWO MAXIMUM FREQS [CONFUSION FREQS]
      if(M->SUBS.in != 0)
        HitSUBS(M);
    break;
    default:  // IT HAS ONE MAXIMUM FREQ
      if(M->SUBS.in == 0){ // IF IS OUT
        M->SUBS.in = 1;
        memset(M->SUBS.mask, 0, M->ctx);
        }
      else{ // IF IS IN
        if(best == sym) HitSUBS(M);
        else{
          FailSUBS(M);
          M->SUBS.seq->buf[M->SUBS.seq->idx] = best; 
          } // UPDATE BUFFER WITH NEW SYMBOL
        }
    }
  UpdateCBuffer(M->SUBS.seq);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CorrectCModelADDS(CModel *M, PModel *P, uint8_t sym){
  int32_t status = 0, best = BestId(P->freqs, P->sum);
  switch(best){
    case -2:  // IT IS A ZERO COUNTER [NOT SEEN BEFORE]
      //if(M->ADDS.in != 0){ FailADDS(M); status = 1; }
      if(M->ADDS.in != 0){ HitADDS(M); status = 0; }
    break;
    case -1:  // IT HAS AT LEAST TWO MAXIMUM FREQS [CONFUSION FREQS]
      if(M->ADDS.in != 0){ HitADDS(M); status = 0; }
    break;
    default:  // IT HAS ONE MAXIMUM FREQ
      if(M->ADDS.in == 0){ // IF IS OUT
        M->ADDS.in = 1;
        memset(M->ADDS.mask, 0, M->ctx);
        status = 0;
        }
      else{ // IF IS IN
        if(best == sym){ HitADDS(M);  status = 0; }
        else           { FailADDS(M); status = 1; }
        }
    }
  if(status == 0) UpdateCBuffer(M->ADDS.seq); else M->ADDS.idx = M->ADDS.idx2;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CorrectCModelDELS(CModel *M, PModel *P, uint8_t sym){
  int32_t status = 0, best = BestId2(P->freqs, P->sum);

  if(M->DELS.in == 0 && best != -1){ // IF IS OUT
    M->DELS.in = 1;
    memset(M->DELS.mask, 0, M->ctx);
    status = 0;
    UpdateCBuffer(M->DELS.seq); 
    return;    
    }

  if(best == -1)
    best = 0;

  if(M->DELS.in == 1){ // IF IS IN
    if(best == sym){ 
      HitDELS(M);  
      status = 0; 
      }
    else{ 
      FailDELS(M); 
      status = 1; 
      }
    }
  
  if(status == 1){
    M->DELS.idx = M->DELS.idx2;
    M->DELS.seq->buf[M->DELS.seq->idx] = best;
    M->DELS.idx = GetPModelIdxCorr(M->DELS.seq->buf+M->DELS.seq->idx-1, M, 
    M->DELS.idx);
    UpdateCBuffer(M->DELS.seq); 

    M->DELS.seq->buf[M->DELS.seq->idx] = sym;
    M->DELS.idx = GetPModelIdxCorr(M->DELS.seq->buf+M->DELS.seq->idx-1, M,
    M->DELS.idx);
    }

  UpdateCBuffer(M->DELS.seq); 
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void ComputePModel(CModel *M, PModel *P, uint64_t idx, uint32_t aDen){
  ACC *ac;
  switch(M->mode){
    case HASH_TABLE_MODE:
      GetHCCounters(&M->hTable, ZHASH(idx), P, aDen);
    break;
    case ARRAY_MODE:
      ac = &M->array.counters[idx<<2];
      P->freqs[0] = 1 + aDen * ac[0];
      P->freqs[1] = 1 + aDen * ac[1];
      P->freqs[2] = 1 + aDen * ac[2];
      P->freqs[3] = 1 + aDen * ac[3];
      P->sum = P->freqs[0] + P->freqs[1] + P->freqs[2] + P->freqs[3];
    break;
    default:
    fprintf(stderr, "Error: not implemented!\n");
    exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void ComputeWeightedFreqs(double w, PModel *P, FloatPModel *PT){
  double f = w / P->sum;
  PT->freqs[0] += (double) P->freqs[0] * f;
  PT->freqs[1] += (double) P->freqs[1] * f;
  PT->freqs[2] += (double) P->freqs[2] * f;
  PT->freqs[3] += (double) P->freqs[3] * f;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double PModelSymbolNats(PModel *P, U32 s){
  return log((double) P->sum / P->freqs[s]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
