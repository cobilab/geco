#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "mem.h"
#include "defs.h"
#include "common.h"
#include "context.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - -

uint64_t Compress(Parameters *P, CModel **cModels, uint8_t id, uint32_t 
refNModels)
  {
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".co");
  FILE        *Writter = Fopen(name, "w");
  uint64_t    nSymbols = 0;
  uint32_t    n, s, k, idxPos;
  int32_t     idx = 0;
  double      *cModelWeight, cModelTotalWeight = 0;
  uint8_t     *readerBuffer, *symbolBuffer, sym, irSym;
  PModel      **pModel, *MX;
  FloatPModel *floatPModel;
  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stderr, "Compressing target sequence %d ...\n", id + 1);

  _bytes_output = 0;
  nSymbols      = NDNASyminFile(Reader);
  pModel        = (PModel  **) Calloc(P->nModels, sizeof(PModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  floatPModel   = CreateFloatPModel(ALPHABET_SIZE);
  readerBuffer  = (uint8_t  *) Calloc(BUFFER_SIZE,          sizeof(uint8_t));
  symbolBuffer  = (uint8_t  *) Calloc(BUFFER_SIZE + BGUARD, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModelWeight  = (double   *) Calloc(P->nModels,           sizeof(double  ));

  for(n = 0 ; n < P->nModels ; ++n)
    {
    cModelWeight[n] = 1.0 / P->nModels;
    if(P->model[n].type == TARGET)
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, 
      P->model[n].ir, TARGET);
    }

  startoutputtingbits();
  start_encode();

  WriteNBits(WATERMARK,                32, Writter);
  WriteNBits(P->checksum,              46, Writter);
  WriteNBits(nSymbols,                 46, Writter);
  WriteNBits((int) (P->gamma * 65536), 32, Writter);
  WriteNBits(P->nModels,               16, Writter);
  for(n = 0 ; n < P->nModels ; ++n)
    {
    WriteNBits(cModels[n]->ctx,        16, Writter);
    WriteNBits(cModels[n]->alphaDen,   16, Writter);
    WriteNBits(cModels[n]->ir,          1, Writter);
    WriteNBits(P->model[n].type,        1, Writter);
    }

  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = DNASymToNum(readerBuffer[idxPos]);
      if(sym == 4) continue;
      symbolBuffer[idx] = sym;

      memset((void *) floatPModel->freqs, 0, ALPHABET_SIZE * sizeof(double));

      for(n = 0 ; n < P->nModels ; ++n)
        {
        GetPModelIdx(symbolBuffer+idx-1, cModels[n]);
        ComputePModel(cModels[n], pModel[n]);

        // The probabilities estimated by each cModel are weighted
        // according to the set of current weights.
        for(s = 0 ; s < ALPHABET_SIZE ; ++s)
          floatPModel->freqs[s] += (double) pModel[n]->freqs[s] /
          pModel[n]->sum * cModelWeight[n];
        }

      MX->sum = 0;
      for(s = 0 ; s < ALPHABET_SIZE ; ++s)
        {
        MX->freqs[s] = 1 + (unsigned) (floatPModel->freqs[s] * MX_PMODEL); 
        MX->sum     += MX->freqs[s];
        }

      ArithEncodeSymbol(sym, (int *)(MX->freqs), (int) MX->sum, Writter);

      cModelTotalWeight = 0;
      for(n = 0 ; n < P->nModels ; ++n)
        {
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;

        cModelTotalWeight += cModelWeight[n];
        
        //UPDATE REMOVE COUNTER
//        if(cModels[n]->ref == REFERENCE)
//        UpdateCModelCounterRM(cModels[n], sym);

        if(cModels[n]->ref == TARGET)
          {
          UpdateCModelCounter(cModels[n], sym);
          if(cModels[n]->ir == 1)                           // Inverted repeats
            {
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounterIr(cModels[n], irSym);
            }
          }
        }

      // Re-normalize the weights
      for(n = 0 ; n < P->nModels ; ++n)
        cModelWeight[n] /= cModelTotalWeight;

      if(++idx == BUFFER_SIZE)
        {
        memcpy(symbolBuffer - BGUARD, symbolBuffer + idx - BGUARD, BGUARD);
        idx = 0;
        }

      #ifdef PROGRESS
      CalcProgress(nSymbols, ++i);
      #endif
      }

  finish_encode(Writter);
  doneoutputtingbits(Writter);
  fclose(Writter);

  Free(MX);
  Free(name);
  Free(cModelWeight);
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      FreeCModel(cModels[n]);
  for(n = 0 ; n < P->nModels ; ++n)
    {
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
  Free(floatPModel->freqs);
  Free(floatPModel);
  Free(readerBuffer);
  Free(symbolBuffer-BGUARD);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID 

  return _bytes_output;
  }


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

CModel **LoadReference(Parameters *P)
  {
  FILE      *Reader = Fopen(P->ref, "r");
  uint32_t  n, k, idxPos;
  int32_t   idx = 0;
  uint8_t   *readerBuffer, *symbolBuffer, sym, irSym;
  CModel    **cModels;
  #ifdef PROGRESS
  uint64_t  i = 0, size = NBytesInFile(Reader);
  #endif

  if(P->verbose == 1)
    fprintf(stderr, "Building reference model ...\n");

  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  symbolBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + BGUARD+1, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModels       = (CModel **) Malloc(P->nModels * sizeof(CModel *)); 
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, 
      P->model[n].ir, REFERENCE);

  P->checksum   = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = DNASymToNum(readerBuffer[idxPos]);
      if(sym == 4)
        continue;
      symbolBuffer[idx] = sym;
      P->checksum = (P->checksum + (uint8_t) sym);

      for(n = 0 ; n < P->nModels ; ++n)
        if(P->model[n].type == REFERENCE)
          {
          GetPModelIdx(symbolBuffer+idx-1, cModels[n]);
          UpdateCModelCounter(cModels[n], sym);
          if(cModels[n]->ir == 1)                          // Inverted repeats
            {
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounterIr(cModels[n], irSym);
            }
          }

      if(++idx == BUFFER_SIZE)
        {
        memcpy(symbolBuffer - BGUARD, symbolBuffer + idx - BGUARD, BGUARD);
        idx = 0;
        }
      #ifdef PROGRESS
      CalcProgress(size, ++i);
      #endif
      }
 
  P->checksum %= CHECKSUMGF; 
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
  Free(readerBuffer);
  Free(symbolBuffer-BGUARD);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID  

  return cModels;
  }
  
//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[])
  {
  char        **p = *&argv;
  CModel      **refModels;
  uint32_t    n, k, refNModels;
  uint64_t    totalBytes;
  Parameters  *P;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2)
    {
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "Usage: GeCo [OPTIONS]... -r [FILE]  [FILE]:[...]    \n");
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "  -v                       verbose mode             \n");
    fprintf(stderr, "  -f                       force (be sure!)         \n");
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "  -rm <ctx>:<den>:<ir>     reference context model  \n");
    fprintf(stderr, "  -rm <ctx>:<den>:<ir>     reference context model  \n");
    fprintf(stderr, "  ...                                               \n");
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "  -tm <ctx>:<den>:<ir>     target context model     \n");
    fprintf(stderr, "  -tm <ctx>:<den>:<ir>     target context model     \n");
    fprintf(stderr, "  ...                                               \n");
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "  -g  <gamma>              gamma factor             \n");
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "  -r  <rFile>              reference file           \n");
    fprintf(stderr, "                                                    \n");
    fprintf(stderr, "  <tFile1>:<tFile2>:<...>  target file(s)         \n\n");
    return EXIT_SUCCESS;
    }

  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v" );
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-f" );

  P->nModels  = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0 || strcmp(argv[n], "-tm") == 0)
      P->nModels += 1;

  if(P->nModels == 0)
    {
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return 1;
    }

  P->model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));

  k = 0;
  refNModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0)
      {
      P->model[k++] = ArgsUniqModel(argv[n+1], 1);
      ++refNModels;
      }
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-tm") == 0)
      P->model[k++] = ArgsUniqModel(argv[n+1], 0);

  P->gamma = ArgsDouble (DEFAULT_GAMMA, p, argc, "-g");
  P->gamma = ((int)(P->gamma * 65536)) / 65536.0;
  P->ref   = ArgsString (NULL, p, argc, "-r");
  P->nTar  = ReadFNames (P, argv[argc-1]);
  if(P->verbose) 
    PrintArgs(P);

  if(refNModels == 0)
    refModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  else
    {
    if(P->ref == NULL)
      {
      fprintf(stderr, "Error: using reference model(s) in nonexistent "
      "reference sequence!\n");
      exit(1);
      }
    refModels = LoadReference(P);
    if(P->verbose)
      fprintf(stderr, "Checksum: %"PRIu64"\n", P->checksum);
    }

  uint64_t bytes[P->nTar];

  totalBytes = 0;
  for(n = 0 ; n < P->nTar ; ++n)
    {
    bytes[n]    = Compress(P, refModels, n, refNModels);
    totalBytes += bytes[n];
    }

  if(P->nTar > 1)
    for(n = 0 ; n < P->nTar ; ++n)
      fprintf(stderr, "File %d compressed bytes: %"PRIu64"\n", n+1, (uint64_t) 
      bytes[n]);

  fprintf(stderr, "Total bytes: %"PRIu64"\n", (uint64_t) totalBytes);  

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
