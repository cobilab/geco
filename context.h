#ifndef CONTEXT_H_INCLUDED
#define CONTEXT_H_INCLUDED

#include "defs.h"

#define ARRAY_MODE            0
#define HASH_TABLE_MODE       1
#define HASH_TABLE_BEGIN_CTX  15
#define HASH_SIZE             33554471
#define MAX_HASH_CTX          28 

typedef unsigned short  ACCounter;      // Size of context counters for arrays
typedef unsigned char   HCCounter; // Size of context counters for hash tables
typedef HCCounter       HCCounters[4];

typedef struct
  {
  uint32_t        key;                         // The key stored in this entry
  unsigned char   counters;           // "Small" counters: 2 bits for each one
  }
Entry;

typedef struct
  {
  unsigned short  *entrySize;                  // Number of keys in this entry
  Entry           **entries;              // The heads of the hash table lists
  HCCounters      **counters;                          // The context counters
  }
HashTable;

typedef struct
  {
  ACCounter       *counters;
  }
Array;

typedef struct
  {
  unsigned        ctx;                    // Current depth of context template
  ULL             nPModels;            // Maximum number of probability models
  unsigned        alphaDen;                            // Denominator of alpha
  unsigned        maxCount;        // Counters /= 2 if one counter >= maxCount
  uint64_t        multiplier;
  uint64_t        pModelIdx;
  uint64_t        pModelIdxIR;
  unsigned        mode;
  uint8_t         ir;
  uint8_t         ref;
  HashTable       hTable;
  Array           array;
  }
CModel;

typedef struct
  {
  unsigned        *freqs;
  unsigned        sum;
  }
PModel;

typedef struct
  {
  double          *freqs;
  }
FloatPModel;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void           FreeCModel            (CModel *);
inline void    GetPModelIdx          (uint8_t *, CModel *);
inline uint8_t GetPModelIdxIR        (uint8_t *, CModel *);
PModel         *CreatePModel         (unsigned);
FloatPModel    *CreateFloatPModel    (unsigned);
void           ResetCModelIdx        (CModel *);
void           UpdateCModelCounter   (CModel *, unsigned);
void           UpdateCModelCounterIr (CModel *, unsigned);
CModel         *CreateCModel         (uint32_t, uint32_t, uint32_t, uint8_t);
void           ComputePModel         (CModel *, PModel *);
double         PModelSymbolNats      (PModel *, unsigned);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
