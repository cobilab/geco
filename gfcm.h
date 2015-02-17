#ifndef GFCM_H_INCLUDED
#define GFCM_H_INCLUDED

#include "defs.h"

typedef uint16_t   GACC;       // SIZE OF COUNTERS FOR COUNTER TABLE [8|16|32]
#define MAXGACC_C  (((uint64_t)1<<(sizeof(GACC)*8))-1)

typedef struct{
  GACC     *cnts;     // TABLE COUNTERS FOR GLOBAL MODEL
  }
GARRAY;

typedef struct{
  uint64_t nPMod;     // MAXIMUM NUMBER OF PROBABILITY MODELS
  uint64_t mult;      // MULTIPLICATOR TO CALCULATE INDEX
  uint64_t idx;       // CURRENT CONTEXT INDEX
  uint32_t aDen;      // ESTIMATOR DENOMINATOR
  uint32_t ctx;       // CONTEXT ORDER FOR THE FCM
  uint32_t *freqs;    // FCM SYMBOL PROBABILITIES
  uint8_t  nSym;      // FCM NUMBER OF SYMBOLS
  GARRAY   A;         // COUNTER TABLE LINK
  }
GFCM;

void        FreeGFCM     (GFCM *);
void        UpdateGFCM   (GFCM *, uint32_t);
GFCM        *CreateGFCM  (uint32_t, uint32_t, uint8_t);
void        ComputeGFCM  (GFCM *);
inline void GetIdx       (uint8_t *, GFCM *);

#endif
