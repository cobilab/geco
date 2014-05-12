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
// - - - - - - - - - - - - - - D E C O M P R E S S O R - - - - - - - - - - - -

void Decompress(Parameters *P, CModel **cModels, uint8_t id)
  {
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".de");
  FILE        *Writter = Fopen(name, "w");
  uint64_t    bases = 0;
  uint32_t    n, s, k;
  double      *cModelWeight, cModelTotalWeight = 0;
  int32_t     idx = 0, idxOut = 0;
  uint8_t     *outBuffer, *symbolBuffer, sym = 0, irSym = 0;
  PModel      **pModel, *MX;
  FloatPModel *floatPModel;
  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stderr, "Decompressing %"PRIu64" bases of target %d ...\n", 
    P[id].size, id + 1);

  startinputtingbits();
  start_decode(Reader);

  P[id].watermark        = ReadNBits(32, Reader);
  garbage                = ReadNBits(32, Reader);
  P[id].size             = ReadNBits(32, Reader);
  P[id].gamma            = ReadNBits(32, Reader) / 65536.0;
  P[id].nModels          = ReadNBits(16, Reader);
  for(k = 0 ; k < P[id].nModels ; ++k)
    {
    P[id].model[k].ctx   = ReadNBits(16, Reader);
    P[id].model[k].den   = ReadNBits(16, Reader);
    P[id].model[k].ir    = ReadNBits( 1, Reader);
    P[id].model[k].type  = ReadNBits( 1, Reader);
    }

  bases         = P[id].size;
  pModel        = (PModel  **) Calloc(P[id].nModels, sizeof(PModel *));
  for(n = 0 ; n < P[id].nModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  floatPModel   = CreateFloatPModel(ALPHABET_SIZE);
  outBuffer     = (uint8_t  *) Calloc(BUFFER_SIZE,          sizeof(uint8_t));
  symbolBuffer  = (uint8_t  *) Calloc(BUFFER_SIZE + BGUARD, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModelWeight  = (double   *) Calloc(P[id].nModels,        sizeof(double ));

  for(n = 0 ; n < P[id].nModels ; ++n)
    {
    cModelWeight[n] = 1.0 / P[id].nModels;
    if(P[id].model[n].type == TARGET)
      cModels[n] = CreateCModel(P[id].model[n].ctx , P[id].model[n].den, 
      P[id].model[n].ir, TARGET);
    }

  while(bases--)
    {
    memset((void *) floatPModel->freqs, 0, ALPHABET_SIZE * sizeof(double));

    for(n = 0 ; n < P[id].nModels ; ++n)
      {
      GetPModelIdx(symbolBuffer + idx - 1, cModels[n]);
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

    symbolBuffer[idx] = sym = ArithDecodeSymbol(ALPHABET_SIZE, (int *) 
    MX->freqs, (int) MX->sum, Reader);
    outBuffer[idxOut] = NumToDNASym(sym);

    cModelTotalWeight = 0;
    for(n = 0 ; n < P[id].nModels ; ++n)
      {
      cModelWeight[n] = Power(cModelWeight[n], P[id].gamma) * (double)
      pModel[n]->freqs[sym] / pModel[n]->sum;

      cModelTotalWeight += cModelWeight[n];

      if(P[id].model[n].type == TARGET)
        {
        UpdateCModelCounter(cModels[n], sym);
        if(cModels[n]->ir == 1)                          // Inverted repeats
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

    if(++idxOut == BUFFER_SIZE)
      {
      fwrite(outBuffer, 1, idxOut, Writter);
      idxOut = 0;
      }

    #ifdef PROGRESS
    CalcProgress(P[id].size, ++i);
    #endif
    }
  if(idxOut != 0) 
    fwrite(outBuffer, 1, idxOut, Writter);

  finish_decode();
  doneinputtingbits();
  fclose(Writter);

  Free(MX);
  Free(name);
  Free(cModelWeight);
  for(n = 0 ; n < P[id].nModels ; ++n)
    if(P[id].model[n].type == REFERENCE)
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
  Free(outBuffer);
  Free(symbolBuffer-BGUARD);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID 
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

  cModels    = (CModel **) Malloc(P->nModels * sizeof(CModel *)); 
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
  uint32_t    n, k, *checksum, refNModels = 0;
  Parameters  *P;
  FILE        *Reader = NULL;
  uint8_t     help, verbose, force, nTar = 1;
  
  if((help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2)
    {
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage: GeDe [OPTIONS]... -r [FILE]  [FILE]:[...]   \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, " -v                       verbose mode             \n");
    fprintf(stderr, " -f                       force (be sure!)         \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, " -r  <rFile>              reference file           \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, " <tFile1>:<tFile2>:<...>  target file(s)         \n\n");
    return EXIT_SUCCESS;
    }

  verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v");
  force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-f");

  for(n = 0 ; n != strlen(argv[argc-1]) ; ++n)
    if(argv[argc-1][n] == ':')
      ++nTar;

  P        = (Parameters *) Malloc(nTar * sizeof(Parameters));
  checksum = (uint32_t   *) Calloc(nTar , sizeof(uint32_t));

  P[0].force   = force;
  P[0].verbose = verbose;
  P[0].nTar    = ReadFNames (P, argv[argc-1]);
  P[0].ref     = ArgsString (NULL, p, argc, "-r");
  for(n = 0 ; n < nTar ; ++n)
    {
    Reader = Fopen(P[0].tar[n], "r");
    startinputtingbits();
    start_decode(Reader);

    refNModels = 0;
    P[n].watermark = ReadNBits(32, Reader);
    if(P[n].watermark != WATERMARK)
      {
      fprintf(stderr, "Error: Invalid compressed file to decompress!\n");
      return 1;
      }
    checksum[n]    = ReadNBits(32, Reader);
    P[n].size      = ReadNBits(32, Reader);
    P[n].gamma     = ReadNBits(32, Reader) / 65536.0;
    P[n].nModels   = ReadNBits(16, Reader);
    P[n].model     = (ModelPar *) Calloc(P[n].nModels, sizeof(ModelPar));
    for(k = 0 ; k < P[n].nModels ; ++k)
      {
      P[n].model[k].ctx  = ReadNBits(16, Reader); 
      P[n].model[k].den  = ReadNBits(16, Reader); 
      P[n].model[k].ir   = ReadNBits( 1, Reader); 
      P[n].model[k].type = ReadNBits( 1, Reader);
      if(P[n].model[k].type == 1)
        ++refNModels;
      }

    finish_decode();
    doneinputtingbits();
    fclose(Reader);
    }

  if(P->verbose)
    PrintArgs(P);
 
  if(refNModels > 0 && P[0].ref == NULL)
    {
    fprintf(stderr, "Error: using reference model(s) in nonexistent "
    "reference sequence!\n");
    exit(1);
    }

  if(refNModels != 0)
    refModels = LoadReference(P);
  else
    refModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));

  if(P->verbose && refNModels != 0)
    fprintf(stderr, "Checksum: %"PRIu64"\n", P->checksum); 

  for(n = 0 ; n < nTar ; ++n)
    {
    if(refNModels != 0)
      {
      if(CmpCheckSum(checksum[n], P[0].checksum) == 0)
        Decompress(P, refModels, n);
      }
    else
      Decompress(P, refModels, n);
    }

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
