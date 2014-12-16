#ifndef CONTEXT_H_INCLUDED
#define CONTEXT_H_INCLUDED

#include "defs.h"

#define ARRAY_MODE            0
#define HASH_TABLE_MODE       1
#define HASH_TABLE_BEGIN_CTX  15
#define HASH_SIZE             33554471
#define MAX_HASH_CTX          28 

typedef U16  ACCounter;            // Size of context counters for arrays
typedef U8   HCCounter;       // Size of context counters for hash tables
typedef U32  KEYSMAX;                                  // Keys index bits
typedef U16  ENTMAX;                // Entry size (nKeys for each hIndex)
typedef HCCounter HCCounters[4];

typedef struct
  {
  U32        key;                         // The key stored in this entry
  HCCounter  counters;           // "Small" counters: 2 bits for each one
  }
Entry;

typedef struct
  {
  ENTMAX     *size;                       // Number of keys in this entry
  Entry      **entries;              // The heads of the hash table lists
  HCCounters **counters;                          // The context counters
  }
HashTable;

typedef struct
  {
  ACCounter  *counters;
  }
Array;

typedef struct
  {
  U32        ctx;                    // Current depth of context template
  U64        nPModels;            // Maximum number of probability models
  U32        alphaDen;                            // Denominator of alpha
  U32        maxCount;        // Counters /= 2 if one counter >= maxCount
  U64        multiplier;
  U64        pModelIdx;
  U64        pModelIdxIR;
  U32        mode;
  U8         ir;
  U8         ref;
  HashTable  hTable;
  Array      array;
  }
CModel;

typedef struct
  {
  U32        *freqs;
  U32        sum;
  }
PModel;

typedef struct
  {
  double     *freqs;
  }
FloatPModel;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void         FreeCModel            (CModel *);
inline void  GetPModelIdx          (U8 *, CModel *);
inline U8    GetPModelIdxIR        (U8 *, CModel *);
PModel       *CreatePModel         (U32);
FloatPModel  *CreateFloatPModel    (U32);
void         ResetCModelIdx        (CModel *);
void         UpdateCModelCounter   (CModel *, U32);
void         UpdateCModelCounterIr (CModel *, U32);
void         UpdateCModelCounterRM (CModel *, U32);
CModel       *CreateCModel         (U32, U32, U32, U8);
void         ComputePModel         (CModel *, PModel *);
double       PModelSymbolNats      (PModel *, U32);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
