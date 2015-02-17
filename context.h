#ifndef CONTEXT_H_INCLUDED
#define CONTEXT_H_INCLUDED

#include "defs.h"

//#define PREC32B // UNCOMMENT: CONTEXTS UP TO 28 (WILL USE MORE MEMORY!)
#define FSEARCHMODE
//#define SWAP

#define ARRAY_MODE            0
#define HASH_TABLE_MODE       1
#define HASH_TABLE_BEGIN_CTX  15
#define HASH_SIZE             33554471
#define MAX_COLLISIONS        10

#ifdef PREC32B
  #define MAX_HASH_CTX        28 
#else
  #define MAX_HASH_CTX        20 
#endif

typedef U16  ACC;                  // Size of context counters for arrays
typedef U8   HCC;             // Size of context counters for hash tables
typedef U16  ENTMAX;                // Entry size (nKeys for each hIndex)
typedef HCC  HCCounters[4];

typedef struct{
  #ifdef PREC32B
  U32        key;                         // The key stored in this entry
  #else
  U16        key;
  #endif
  HCC        counters;           // "Small" counters: 2 bits for each one
  }
Entry;

typedef struct{
  ENTMAX     *index;                      // Number of keys in this entry
  Entry      **entries;              // The heads of the hash table lists
  uint32_t   maxC;
  uint32_t   maxH;
  }
HashTable;

typedef struct {
  ACC        *counters;
  }
Array;

typedef struct{
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

typedef struct{
  U32        *freqs;
  U32        sum;
  }
PModel;

typedef struct{
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
CModel       *CreateCModel         (U32, U32, U32, U8, U32);
void         ComputePModel         (CModel *, PModel *);
double       PModelSymbolNats      (PModel *, U32);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
