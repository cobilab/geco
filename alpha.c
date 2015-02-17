#include <stdio.h>
#include <stdlib.h>
#include "alpha.h"
#include "common.h"
#include "buffer.h"
#include "mem.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t S2NAlpha(uint8_t s, Alpha *A){
  switch(s){
    case 'A': return 0;
    case 'T': return 3;
    case 'C': return 1;
    case 'G': return 2;
    default : return A->lowBase;
    }
  }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// INITIALIZE ALPHABETS
//
Alpha *CreateAlphabet(void){
  Alpha *A    = (Alpha   *) Calloc(1,         sizeof(Alpha  ));
  A->symbolic = (uint8_t *) Calloc(MAX_ALPHA, sizeof(uint8_t));
  A->numeric  = (uint8_t *) Calloc(MAX_ALPHA, sizeof(uint8_t));
  A->bin      = (uint8_t *) Calloc(MAX_ALPHA, sizeof(uint8_t));
  A->length   = 0;
  A->nSym     = 0;
  return A;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// FREE ALPHABETS
//
void FreeAlphabet(Alpha *A){
  Free(A->symbolic);
  Free(A->numeric);
  Free(A->bin);
  Free(A);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PRINT FILE INFORMATION [PROPERTIES & ALPHABET]
//
void PrintStreamInfo(Alpha *A){
  uint32_t k;
  fprintf(stderr, "Original file has got .............. "); 
  PrintHRBytes(A->length);
  fprintf(stderr, " (%"PRIu64" Bytes)\n", A->length);
  fprintf(stderr, "Genomic bases: A,C,G,T");
  if(A->Ns == 1)
    fprintf(stderr, ",N");
  fprintf(stderr, "\n");
  fprintf(stderr, "Extra:\n");
  fprintf(stderr, "  [-] Cardinality .................. %u\n", A->nSym);
  fprintf(stderr, "  [-] Alphabet ..................... ");
  for(k = 0 ; k < A->nSym ; ++k)
    A->symbolic[k] == '\n' ? fprintf(stderr, "\\n") :
    fprintf(stderr, "%c", A->symbolic[k]);
  fprintf(stderr, "\n");
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// BUILD ALPHABET
//
void BuildAlphabet(Alpha *A){
  uint32_t k;
  for(k = 0 ; k < MAX_ALPHA ; ++k)
    if(A->bin[k] == 1){
      A->symbolic[A->nSym] = k;
      A->numeric[k] = A->nSym++;
      }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PARSE FILE & BUILD ALPHABET
//
void LoadAlphabet(FILE *F, Alpha *A){
  int64_t  k;
  uint64_t i, counts[4] = { 0, 0, 0, 0}, min, n;
  uint8_t  s, *buf = (uint8_t *) Calloc(BUF_SIZE, sizeof(uint8_t));

  while((k = fread(buf, 1, BUF_SIZE, F)))
    for(i = 0 ; i < k ; ++i){
      s = *(buf+i);
      switch(s){
        case 'A': counts[0]++; break;
        case 'C': counts[1]++; break;
        case 'G': counts[2]++; break;
        case 'T': counts[3]++; break;
        }
      A->bin[s] = 1; 
      A->length++; 
      }

  A->lowBase = 0;
  min = counts[0];
  for(n = 1 ; n < 4 ; ++n)
    if(counts[n] < min){
      min = counts[n];
      A->lowBase = n;
      }

  rewind(F); 
  if(A->bin['N' ] == 1) A->Ns = 1; 
  if(A->bin['\n'] == 1) A->NL = 1; 
  // REMOVE "ACGT" SYMBOLS [THEY ARE HANDLED SEPARATELY] AND RE-CALC ALPHABET
  A->bin['A'] = A->bin['C'] = A->bin['G'] = A->bin['T'] = 0;
  A->bin['N'] = A->bin['\n'] = 0;
  A->nSym  = 0;
  BuildAlphabet(A);
  PrintStreamInfo(A); 

  Free(buf);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

