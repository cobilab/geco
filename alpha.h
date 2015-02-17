#ifndef ALPHA_H_INCLUDED
#define ALPHA_H_INCLUDED

#include "defs.h"

#define MAX_ALPHA ((1<<(sizeof(char)*8))-1)
#define BUF_SIZE  1048576

typedef struct{
  uint8_t  *symbolic;
  uint8_t  *numeric;
  uint8_t  *bin;
  uint64_t length;
  uint32_t nSym;
  uint8_t  Ns;
  uint8_t  NL;
  uint8_t  lowBase;
  }
Alpha;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t S2NAlpha        (uint8_t, Alpha *);
Alpha   *CreateAlphabet (void);
void    FreeAlphabet    (Alpha *);
void    PrintStreamInfo (Alpha *);
void    BuildAlphabet   (Alpha *);
void    LoadAlphabet    (FILE *, Alpha *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
