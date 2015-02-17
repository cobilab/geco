#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "mem.h"
#include "defs.h"
#include "buffer.h"
#include "alpha.h"
#include "common.h"
#include "context.h"
#include "gfcm.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - C O M P R E S S   R A W   S T R E A M - - - - - - - - - -

void CompressStream(FILE *F, GFCM *M, CBUF *B, uint8_t sym, uint8_t nSym){
  B->buf[B->idx] = sym;
  GetIdx(B->buf+B->idx-1, M);
  ComputeGFCM(M);
  AESym(sym, (int *) M->freqs, (int) M->freqs[nSym], F);
  UpdateGFCM(M, sym);
  UpdateCBuffer(B);
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void Compress(Parameters *P, CModel **cModels, uint8_t id, uint32_t 
refNModels, INF *I){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".co");
  FILE        *Writter = Fopen(name, "w");
  uint32_t    n, k, idxPos;
  int32_t     idx = 0;
  double       *cModelWeight, cModelTotalWeight = 0, mA, mC, mG, mT;
  uint8_t     *readerBuffer, *symbolBuffer, sym, irSym, *pos, extra = 0;
  PModel      **pModel, *MX;
  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stderr, "Analyzing data, creating alphabet and models ...\n");
 
  _bytes_output = 0;

  #ifdef EXTRA
  uint8_t ss;
  extra = 1;
  Alpha *A = CreateAlphabet();
  LoadAlphabet(Reader, A);
  uint64_t nSymbols = A->length;
  GFCM *EMod = CreateGFCM(EXTRA_MOD_CTX, EXTRA_MOD_DEN, A->nSym);
  GFCM *BMod = CreateGFCM(EXTRA_BIN_CTX, EXTRA_BIN_DEN, 2);
  GFCM *NMod = CreateGFCM(EXTRA_N_CTX,   EXTRA_N_DEN,   2);
  GFCM *LMod = CreateGFCM(EXTRA_L_CTX,   EXTRA_L_DEN,   2);
  CBUF *EBuf = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);
  CBUF *BBuf = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);
  CBUF *NBuf = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);
  CBUF *LBuf = CreateCBuffer(DEF_BUF_SIZE, DEF_BUF_GUARD);  
  #else
  extra = 0;
  uint64_t nSymbols = NDNASyminFile(Reader);
  #endif

  pModel        = (PModel  **) Calloc(P->nModels, sizeof(PModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  symbolBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + BGUARD, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModelWeight  = (double   *) Calloc(P->nModels, sizeof(double));

  for(n = 0 ; n < P->nModels ; ++n){
    cModelWeight[n] = 1.0 / P->nModels;
    if(P->model[n].type == TARGET)
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, 
      P->model[n].ir, TARGET, P->col);
    }

  if(P->verbose){
    fprintf(stderr, "Done!\n");
    fprintf(stderr, "Compressing target sequence %d [length: %"PRIu64"] ...\n", 
    id + 1, nSymbols);
    }

  startoutputtingbits();
  start_encode();

  WriteNBits(WATERMARK,                32, Writter);
  WriteNBits(P->checksum,              46, Writter);
  WriteNBits(nSymbols,                 46, Writter);
  WriteNBits((int) (P->gamma * 65536), 32, Writter);
  WriteNBits(P->col,                   32, Writter);
  WriteNBits(extra,                     1, Writter);
  #ifdef EXTRA
  WriteNBits(A->Ns,                     1, Writter);
  WriteNBits(A->NL,                     1, Writter);
  WriteNBits((int) A->lowBase,          8, Writter);
  for(n = 0 ; n < MAX_ALPHA ; ++n)
    WriteNBits(A->bin[n],               8, Writter);
  #endif
  WriteNBits(P->nModels,               16, Writter);
  for(n = 0 ; n < P->nModels ; ++n){
    WriteNBits(cModels[n]->ctx,        16, Writter);
    WriteNBits(cModels[n]->alphaDen,   16, Writter);
    WriteNBits(cModels[n]->ir,          1, Writter);
    WriteNBits(P->model[n].type,        1, Writter);
    }

  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      #ifdef PROGRESS
      CalcProgress(nSymbols, ++i);
      #endif

      #ifndef EXTRA
      if((sym = DNASymToNum(readerBuffer[idxPos])) == 4){
        continue;
        }
      symbolBuffer[idx] = sym;
      #else
      symbolBuffer[idx] = sym = S2NAlpha(ss = readerBuffer[idxPos], A);
      #endif
      mA = mC = mG = mT = 0;

      pos = &symbolBuffer[idx-1];
      for(n = 0 ; n < P->nModels ; ++n){
        GetPModelIdx(pos, cModels[n]);
        ComputePModel(cModels[n], pModel[n]);

        double factor = cModelWeight[n] / pModel[n]->sum;
        mA += (double) pModel[n]->freqs[0] * factor;
        mC += (double) pModel[n]->freqs[1] * factor;
        mG += (double) pModel[n]->freqs[2] * factor;
        mT += (double) pModel[n]->freqs[3] * factor;
        }

      MX->sum  = MX->freqs[0] = 1 + (unsigned) (mA * MX_PMODEL);
      MX->sum += MX->freqs[1] = 1 + (unsigned) (mC * MX_PMODEL);
      MX->sum += MX->freqs[2] = 1 + (unsigned) (mG * MX_PMODEL);
      MX->sum += MX->freqs[3] = 1 + (unsigned) (mT * MX_PMODEL);

      AESym(sym, (int *)(MX->freqs), (int) MX->sum, Writter);
      
      #ifdef ESTIMATE
      if(P->estim == 1)
        fprintf(stdout, "%.3g\n", PModelSymbolNats(MX, sym));
      #endif

      cModelTotalWeight = 0;
      for(n = 0 ; n < P->nModels ; ++n){

        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        
        if(cModels[n]->ref == TARGET){
          UpdateCModelCounter(cModels[n], sym);
          if(cModels[n]->ir == 1){                // REVERSE COMPLEMENTS
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounterIr(cModels[n], irSym);
            }
          }
        }

      for(n = 0 ; n < P->nModels ; ++n) // RENORMALIZE THE WEIGHTS
        cModelWeight[n] /= cModelTotalWeight;

      #ifdef EXTRA
      if(sym == A->lowBase){
        if(ss != 'A' && ss != 'C' && ss != 'G' && ss != 'T'){
          // A POSSIBILITY IS TO USE GOLOMB CODES HERE!
          // https://github.com/anirudhvr/golomb-coding/blob/master/encode.c
          CompressStream(Writter, BMod, BBuf, 1, 2);   //EXTRA EXISTS!
          switch(ss){                  // EVALUATE EXTRA SYMBOLIC SYMBOL
            case 'N':
            CompressStream(Writter, NMod, NBuf, 1, 2);
            break;
            case '\n':
            CompressStream(Writter, NMod, NBuf, 0, 2); // NOT 'N'
            CompressStream(Writter, LMod, LBuf, 1, 2); // FOUND '\n'
            break;
            default:
            CompressStream(Writter, NMod, NBuf, 0, 2); // NOT 'N'
            CompressStream(Writter, LMod, LBuf, 0, 2); // NOT '\n'
            CompressStream(Writter, EMod, EBuf, A->numeric[ss], EMod->nSym);
            // COMPRESS THE REMAINING EXTRA SYMBOLS
            }
          }
        else
          CompressStream(Writter, BMod, BBuf, 0, 2); // NOT FOUND EXTRA
        }
      #endif

      if(++idx == BUFFER_SIZE){
        memcpy(symbolBuffer-BGUARD, symbolBuffer+idx-BGUARD, BGUARD);
        idx = 0;
        }
      }

  finish_encode(Writter);
  doneoutputtingbits(Writter);
  fclose(Writter);

  #ifdef EXTRA
  FreeAlphabet(A);
  RemoveCBuffer(EBuf);
  RemoveCBuffer(BBuf);
  FreeGFCM(EMod);
  FreeGFCM(BMod);
  RemoveCBuffer(NBuf);
  FreeGFCM(NMod);
  RemoveCBuffer(LBuf);
  FreeGFCM(LMod);
  #endif

  Free(MX);
  Free(name);
  Free(cModelWeight);
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      FreeCModel(cModels[n]);
  for(n = 0 ; n < P->nModels ; ++n){
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
  Free(readerBuffer);
  Free(symbolBuffer-BGUARD);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stderr, "Done!                          \n");  // SPACES ARE VALID 

  I[id].bytes = _bytes_output;
  I[id].size  = nSymbols;
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
      P->model[n].ir, REFERENCE, P->col);

  P->checksum   = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = DNASymToNum(readerBuffer[idxPos]);
      if(sym == 4) continue;
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
  uint64_t    totalBytes, totalSize;
  Parameters  *P;
  INF         *I;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    fprintf(stderr, "                                                     \n");
    fprintf(stderr, "Usage: GeCo [OPTIONS]... -r [FILE]  [FILE]:[...]     \n");
    fprintf(stderr, "                                                     \n");
    fprintf(stderr, "  -v                       verbose mode,             \n");
    fprintf(stderr, "  -f                       force (be sure!),         \n");
    fprintf(stderr, "  -rm <ctx>:<den>:<ir>     reference context model,  \n");
    fprintf(stderr, "  -rm <ctx>:<den>:<ir>     reference context model,  \n");
    fprintf(stderr, "  ...                                                \n");
    fprintf(stderr, "  -tm <ctx>:<den>:<ir>     target context model,     \n");
    fprintf(stderr, "  -tm <ctx>:<den>:<ir>     target context model,     \n");
    fprintf(stderr, "  ...                                                \n");
    fprintf(stderr, "  -g  <gamma>              gamma factor,             \n");
    fprintf(stderr, "  -c  <collisions>         maximum hash collisions,  \n");
    #ifdef ESTIMATE
    fprintf(stderr, "  -e                       estimate only,            \n");
    #endif
    fprintf(stderr, "  -r  <rFile>              reference file,           \n");
    fprintf(stderr, "                                                     \n");
    fprintf(stderr, "  <tFile1>:<tFile2>:<...>  target file(s).         \n\n");
    return EXIT_SUCCESS;
    }

  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v" );
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-f" );
  P->estim    = ArgsState  (0,               p, argc, "-e" );
  P->col      = ArgsNum    (MAX_COLLISIONS,  p, argc, "-c", 1, 50000);

  P->nModels  = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0 || strcmp(argv[n], "-tm") == 0)
      P->nModels += 1;

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return 1;
    }

  P->model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));

  k = 0;
  refNModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0){
      P->model[k++] = ArgsUniqModel(argv[n+1], 1);
      ++refNModels;
      }
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-tm") == 0)
      P->model[k++] = ArgsUniqModel(argv[n+1], 0);

  P->gamma    = ArgsDouble (DEFAULT_GAMMA, p, argc, "-g");
  P->gamma    = ((int)(P->gamma * 65536)) / 65536.0;
  P->ref      = ArgsString (NULL, p, argc, "-r");
  P->nTar     = ReadFNames (P, argv[argc-1]);
  P->checksum = 0;
  if(P->verbose) 
    PrintArgs(P);

  if(refNModels == 0)
    refModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  else{
    if(P->ref == NULL){
      fprintf(stderr, "Error: using reference model(s) in nonexistent "
      "reference sequence!\n");
      exit(1);
      }
    refModels = LoadReference(P);
    if(P->verbose)
      fprintf(stderr, "Checksum: %"PRIu64"\n", P->checksum);
    }

  I = (INF *) Calloc(P->nTar, sizeof(INF));

  totalSize  = 0;
  totalBytes = 0;
  for(n = 0 ; n < P->nTar ; ++n){
    Compress(P, refModels, n, refNModels, I);
    totalSize  += I[n].size;
    totalBytes += I[n].bytes;
    }

  if(P->nTar > 1)
    for(n = 0 ; n < P->nTar ; ++n){
      fprintf(stderr, "File %d compressed bytes: %"PRIu64" (", n+1, (uint64_t) 
      I[n].bytes);
      PrintHRBytes(I[n].bytes);
      fprintf(stderr, ") , Distance: %.6g\n", (8.0*I[n].bytes)/(2*I[n].size));
      }

  fprintf(stderr, "Total bytes: %"PRIu64" (", totalBytes);
  PrintHRBytes(totalBytes);
  fprintf(stderr, ") , Distance: %.6g\n", (8.0*totalBytes)/(2*totalSize));  

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
