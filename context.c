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

static void InitHashTable(CModel *cModel)
  { 
  cModel->hTable.entries = (Entry **) Calloc(HASH_SIZE, sizeof(Entry *));
  cModel->hTable.counters = (HCCounters **) Calloc(HASH_SIZE,
  sizeof(HCCounters *));
  cModel->hTable.size = (ENTMAX *) Calloc(HASH_SIZE, sizeof(ENTMAX));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FreeCModel(CModel *cModel)
  {
  U32 k;

  if(cModel->mode == HASH_TABLE_MODE)
    {
    for(k = 0 ; k < HASH_SIZE ; ++k)
      {
      if(cModel->hTable.size[k] != 0)
        Free(cModel->hTable.entries[k]);
      if(cModel->hTable.counters[k] != NULL)
        Free(cModel->hTable.counters[k]);
      }
    Free(cModel->hTable.entries);
    Free(cModel->hTable.counters);
    Free(cModel->hTable.size);
    }
  else // TABLE_MODE
    Free(cModel->array.counters);

  Free(cModel);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InitArray(CModel *cModel)
  {
  cModel->array.counters = (ACCounter *) Calloc(cModel->nPModels << 2, 
  sizeof(ACCounter));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InsertKey(HashTable *hTable, U32 hIndex, U64 idx)
  {
  hTable->entries[hIndex] = (Entry *) Realloc(hTable->entries[hIndex],
  (hTable->size[hIndex] + 1) * sizeof(Entry), sizeof(Entry));
  hTable->entries[hIndex][hTable->size[hIndex]++].key = (U32)(idx/HASH_SIZE);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static void InsertCounters(HashTable *hTable, U32 hIndex, U32 nHCCounters, 
U32 k, U32 smallCounters)
  {
  hTable->counters[hIndex] = (HCCounters *) Realloc(hTable->counters[hIndex], 
  (nHCCounters + 1) * sizeof(HCCounters), sizeof(HCCounters));

  if(k < nHCCounters)
    memmove(hTable->counters[hIndex][k + 1], hTable->counters[hIndex][k],
      (nHCCounters - k) * sizeof(HCCounters));

  hTable->counters[hIndex][k][0] =  smallCounters &  0x03;
  hTable->counters[hIndex][k][1] = (smallCounters & (0x03 << 2)) >> 2;
  hTable->counters[hIndex][k][2] = (smallCounters & (0x03 << 4)) >> 4;
  hTable->counters[hIndex][k][3] = (smallCounters & (0x03 << 6)) >> 6;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static HCCounter *GetHCCounters(HashTable *hTable, U64 key)
  {
  U32 k = 0, n, hIndex = key % HASH_SIZE;

  for(n = 0 ; n < hTable->size[hIndex] ; n++)
    {
    if(((U64) hTable->entries[hIndex][n].key*HASH_SIZE) + hIndex == key)
      {
      switch(hTable->entries[hIndex][n].counters)
        {
        case 0:
        return hTable->counters[hIndex][k];

        default:
        auxCounters[0] =  hTable->entries[hIndex][n].counters &  0x03;
        auxCounters[1] = (hTable->entries[hIndex][n].counters & (0x03<<2))>>2;
        auxCounters[2] = (hTable->entries[hIndex][n].counters & (0x03<<4))>>4;
        auxCounters[3] = (hTable->entries[hIndex][n].counters & (0x03<<6))>>6;
        return auxCounters;
        }
      }

    if(hTable->entries[hIndex][n].counters == 0)
      k++;
    }

  return NULL;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

PModel *CreatePModel(U32 nSymbols)
  {
  PModel *pModel;
  pModel = (PModel *) Malloc(sizeof(PModel));
  pModel->freqs = (U32 *) Malloc(nSymbols * sizeof(U32));
  return pModel;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FloatPModel *CreateFloatPModel(U32 nSymbols)
  {
  FloatPModel *floatPModel;
  floatPModel = (FloatPModel *) Malloc(sizeof(FloatPModel));
  floatPModel->freqs = (double *) Malloc(nSymbols * sizeof(double));
  return floatPModel;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateCModelCounterIr(CModel *M, U32 symbol)
  {
  U32       n;
  ACCounter *aCounters;
  U64       idx = M->pModelIdxIR;

  if(M->mode == HASH_TABLE_MODE)
    {
    U8 smallCounter;
    U32 i, k = 0, nHCCounters, hIndex = idx % HASH_SIZE;
    for(n = 0 ; n < M->hTable.size[hIndex] ; n++)
      {
      if(((U64) M->hTable.entries[hIndex][n].key*HASH_SIZE)+hIndex == idx)
        {
        if(M->hTable.entries[hIndex][n].counters == 0)  // Large counters
          {
          if(++M->hTable.counters[hIndex][k][symbol] == 255)
            {
            M->hTable.counters[hIndex][k][0] >>= 1;
            M->hTable.counters[hIndex][k][1] >>= 1;
            M->hTable.counters[hIndex][k][2] >>= 1;
            M->hTable.counters[hIndex][k][3] >>= 1;
            }
          return;
          }
        
        smallCounter = (M->hTable.entries[hIndex][n].counters>>(symbol<<1))&0x03;
         // If "counters" is non-zero, then this is at least the
         // second time that this key is generated. Therefore,
         // if the "small" counter of the symbol if full (i.e.,
         // is equal to 3), then the "large" counters have to be
         // inserted into the right position.
        if(smallCounter == 3)
          {
          nHCCounters = k;
          for(i = n + 1 ; i < M->hTable.size[hIndex] ; ++i)
            if(M->hTable.entries[hIndex][i].counters == 0)
              nHCCounters++;

          InsertCounters(&M->hTable, hIndex, nHCCounters, k,
          M->hTable.entries[hIndex][n].counters);
          M->hTable.entries[hIndex][n].counters = 0;
          M->hTable.counters[hIndex][k][symbol]++;
          return;
          }
        else // There is still room for incrementing the "small" counter.
          {
          smallCounter++;
          M->hTable.entries[hIndex][n].counters &= ~(0x03<<(symbol<<1));
          M->hTable.entries[hIndex][n].counters |= (smallCounter<<(symbol<<1));
          return;
          }
        }

      // Keeps counting the number of HCCounters in this entry
      if(!M->hTable.entries[hIndex][n].counters)
        k++;
      }

    // If key not found
    InsertKey(&M->hTable, hIndex, idx);
    M->hTable.entries[hIndex][M->hTable.size[hIndex]-1].counters = (0x01<<(symbol<<1));
    }
  else
    {
    aCounters = &M->array.counters[idx << 2];
    aCounters[symbol]++;
    if(aCounters[symbol] == M->maxCount && M->maxCount != 0)
      {    
      aCounters[0] >>= 1;
      aCounters[1] >>= 1;
      aCounters[2] >>= 1;
      aCounters[3] >>= 1;
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void UpdateCModelCounter(CModel *cModel, U32 symbol)
  {
  unsigned  n;
  ACCounter *aCounters;
  uint64_t  pModelIdx = cModel->pModelIdx;

  if(cModel->mode == HASH_TABLE_MODE)
    {
    unsigned char smallCounter;
    unsigned i, k = 0;
    unsigned nHCCounters;            // The number of HCCounters in this entry
    unsigned hIndex = pModelIdx % HASH_SIZE;                 // The hash index

    for(n = 0 ; n < cModel->hTable.size[hIndex] ; n++)
      {
      if(((uint64_t) cModel->hTable.entries[hIndex][n].key * HASH_SIZE) + 
      hIndex == pModelIdx)
        {
        // If "counters" is zero, then update the "large" counters.
        if(cModel->hTable.entries[hIndex][n].counters == 0)
          {
          if(++cModel->hTable.counters[hIndex][k][symbol] == 255)
            {
            cModel->hTable.counters[hIndex][k][0] >>= 1;
            cModel->hTable.counters[hIndex][k][1] >>= 1;
            cModel->hTable.counters[hIndex][k][2] >>= 1;
            cModel->hTable.counters[hIndex][k][3] >>= 1;
            }
          return;
          }
        
        smallCounter = (cModel->hTable.entries[hIndex][n].counters >> (symbol 
        << 1)) & 0x03;
         // If "counters" is non-zero, then this is at least the
         // second time that this key is generated. Therefore,
         // if the "small" counter of the symbol if full (i.e.,
         // is equal to 3), then the "large" counters have to be
         // inserted into the right position.
        if(smallCounter == 3)
          {
          nHCCounters = k;
          for(i = n + 1 ; i < cModel->hTable.size[hIndex] ; i++)
            if(cModel->hTable.entries[hIndex][i].counters == 0)
              nHCCounters++;

          InsertCounters(&cModel->hTable, hIndex, nHCCounters, k,
          cModel->hTable.entries[hIndex][n].counters);
          cModel->hTable.entries[hIndex][n].counters = 0;
          cModel->hTable.counters[hIndex][k][symbol]++;
          return;
          }
        else // There is still room for incrementing the "small" counter.
          {
          smallCounter++;
          cModel->hTable.entries[hIndex][n].counters &= ~(0x03<<(symbol<<1));
          cModel->hTable.entries[hIndex][n].counters |= (smallCounter<<(symbol
          <<1));
          return;
          }
        }

      // Keeps counting the number of HCCounters in this entry
      if(!cModel->hTable.entries[hIndex][n].counters)
        k++;
      }

    // If key not found
    InsertKey(&cModel->hTable, hIndex, pModelIdx);
    cModel->hTable.entries[hIndex][cModel->hTable.size[hIndex]-1].
    counters = (0x01 << (symbol << 1));
    }
  else
    {
    aCounters = &cModel->array.counters[pModelIdx << 2];
    aCounters[symbol]++;
    if(aCounters[symbol] == cModel->maxCount && cModel->maxCount != 0)
      {    
      aCounters[0] >>= 1;
      aCounters[1] >>= 1;
      aCounters[2] >>= 1;
      aCounters[3] >>= 1;
      }
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CModel *CreateCModel(U32 ctx, U32 aDen, U32 ir, U8 ref) 
  {
  CModel *M = (CModel *) Calloc(1, sizeof(CModel));
  U64    prod = 1, *multipliers;
  U32    n;

  if(ctx > MAX_HASH_CTX)
    {
    fprintf(stderr, "Error: context size cannot be greater than %d\n", 
    MAX_HASH_CTX);
    exit(1);
    }
  
  multipliers      = (U64 *) Calloc(ctx, sizeof(U64));
  M->nPModels      = (U64) pow(ALPHABET_SIZE, ctx);
  M->ctx           = ctx;
  M->alphaDen      = aDen;
  M->pModelIdx     = 0;
  M->pModelIdxIR   = M->nPModels - 1;
  M->ir            = ir  == 0 ? 0 : 1;
  M->ref           = ref == 0 ? 0 : 1;

  if(ctx >= HASH_TABLE_BEGIN_CTX)
    {
    M->mode     = HASH_TABLE_MODE;
    M->maxCount = DEFAULT_MAX_COUNT >> 8;
    InitHashTable(M);
    }
  else 
    {
    M->mode     = ARRAY_MODE;
    M->maxCount = DEFAULT_MAX_COUNT;
    InitArray(M);
    }

  for(n = 0 ; n < M->ctx ; ++n)
    {
    multipliers[n] = prod;
    prod <<= 2;
    }

  M->multiplier = multipliers[M->ctx-1];

  return M;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ResetCModelIdx(CModel *M)
  {
  M->pModelIdx   = 0;
  M->pModelIdxIR = M->nPModels - 1;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline U8 GetPModelIdxIR(U8 *p, CModel *M)
  {
  M->pModelIdxIR = (M->pModelIdxIR>>2)+GetCompNum(*p)*M->multiplier;
  return GetCompNum(*(p - M->ctx));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inline void GetPModelIdx(U8 *p, CModel *M)
  {
  M->pModelIdx = ((M->pModelIdx-*(p-M->ctx)*M->multiplier)<<2)+*p;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ComputePModel(CModel *cModel, PModel *pModel)
  {
  HCCounter *hCounters;
  ACCounter *aCounters;

  if(cModel->mode == HASH_TABLE_MODE)
    {
    if(!(hCounters = GetHCCounters(&cModel->hTable, cModel->pModelIdx)))
      hCounters = zeroCounters;
    pModel->freqs[0] = 1 + cModel->alphaDen * hCounters[0];
    pModel->freqs[1] = 1 + cModel->alphaDen * hCounters[1];
    pModel->freqs[2] = 1 + cModel->alphaDen * hCounters[2];
    pModel->freqs[3] = 1 + cModel->alphaDen * hCounters[3];
    }
  else
    {
    aCounters = &cModel->array.counters[cModel->pModelIdx << 2];
    pModel->freqs[0] = 1 + cModel->alphaDen * aCounters[0];
    pModel->freqs[1] = 1 + cModel->alphaDen * aCounters[1];
    pModel->freqs[2] = 1 + cModel->alphaDen * aCounters[2];
    pModel->freqs[3] = 1 + cModel->alphaDen * aCounters[3];
    }

  pModel->sum = pModel->freqs[0] + pModel->freqs[1] + pModel->freqs[2] + 
  pModel->freqs[3];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double PModelSymbolNats(PModel *pModel, U32 symbol)
  {
  return log((double) pModel->sum / pModel->freqs[symbol]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
