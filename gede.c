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
#include "common.h"
#include "context.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - D E C O M P R E S S O R - - - - - - - - - - - -

void Decompress(Parameters *P, CModel **cModels, uint8_t id){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = ReplaceSubStr(P->tar[id], ".co", ".de"); 
  FILE        *Writter = Fopen(name, "w");
  uint64_t    nSymbols = 0;
  uint32_t    n, k;
  double      *cModelWeight, cModelTotalWeight = 0, mA, mC, mG, mT;
  int32_t     idx = 0, idxOut = 0;
  uint8_t     *outBuffer, *symbolBuffer, sym = 0, irSym = 0, *pos;
  PModel      **pModel, *MX;
  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stderr, "Decompressing %"PRIu64" symbols of target %d ...\n", 
    P[id].size, id + 1);

  startinputtingbits();
  start_decode(Reader);

  P[id].watermark        = ReadNBits(32, Reader);
  garbage                = ReadNBits(46, Reader);
  P[id].size             = ReadNBits(46, Reader);
  P[id].gamma            = ReadNBits(32, Reader) / 65536.0;
  P[id].col              = ReadNBits(32, Reader);
  P[id].nModels          = ReadNBits(16, Reader);
  for(k = 0 ; k < P[id].nModels ; ++k){
    P[id].model[k].ctx   = ReadNBits(16, Reader);
    P[id].model[k].den   = ReadNBits(16, Reader);
    P[id].model[k].ir    = ReadNBits( 1, Reader);
    P[id].model[k].am    = ReadNBits( 8, Reader);
    P[id].model[k].type  = ReadNBits( 1, Reader);
    }

  nSymbols      = P[id].size;
  pModel        = (PModel  **) Calloc(P[id].nModels, sizeof(PModel *));
  for(n = 0 ; n < P[id].nModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  outBuffer     = (uint8_t  *) Calloc(BUFFER_SIZE,          sizeof(uint8_t));
  symbolBuffer  = (uint8_t  *) Calloc(BUFFER_SIZE + BGUARD, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModelWeight  = (double   *) Calloc(P[id].nModels,        sizeof(double ));

  for(n = 0 ; n < P[id].nModels ; ++n){
    cModelWeight[n] = 1.0 / P[id].nModels;
    if(P[id].model[n].type == TARGET)
      cModels[n] = CreateCModel(P[id].model[n].ctx , P[id].model[n].den, 
      P[id].model[n].ir, TARGET, P[id].col, P[id].model[n].am);
    }

  while(nSymbols--){
    #ifdef PROGRESS
    CalcProgress(P[id].size, ++i);
    #endif

    mA = mC = mG = mT = 0;
    pos = &symbolBuffer[idx-1];
    for(n = 0 ; n < P[id].nModels ; ++n){
      if(cModels[n]->am == 0)
        GetPModelIdx(pos, cModels[n]);
      else{
        GetPModelIdx(cModels[n]->correct.seq->buf+
        cModels[n]->correct.seq->idx-1, cModels[n]);
        GetPModelIdxCorr(pos, cModels[n]);
        }
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

    symbolBuffer[idx] = sym = ArithDecodeSymbol(ALPHABET_SIZE, (int *) 
    MX->freqs, (int) MX->sum, Reader);
    outBuffer[idxOut] = NumToDNASym(sym);

    for(n = 0 ; n < P[id].nModels ; ++n)
      if(cModels[n]->am != 0)
        cModels[n]->correct.seq->buf[cModels[n]->correct.seq->idx] = sym;

    cModelTotalWeight = 0;
    for(n = 0 ; n < P[id].nModels ; ++n){
      cModelWeight[n] = Power(cModelWeight[n], P[id].gamma) * (double)
      pModel[n]->freqs[sym] / pModel[n]->sum;
      cModelTotalWeight += cModelWeight[n];
      if(P[id].model[n].type == TARGET){
        if(cModels[n]->am != 0)
          UpdateCModelCounter(cModels[n], sym, cModels[n]->correct.idx);
        else
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
        if(cModels[n]->ir == 1){                      // Inverted repeats
          irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
          UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
          }
        }
      }

    for(n = 0 ; n < P->nModels ; ++n) // RENORMALIZE THE WEIGHTS
      cModelWeight[n] /= cModelTotalWeight;

      #ifdef CORRECT
      for(n = 0 ; n < P->nModels ; ++n){
        if(cModels[n]->am != 0){
          int32_t best = BestId(pModel[n]->freqs, pModel[n]->sum);
          switch(best){
            case -2:  // IT IS A ZERO COUNTER [NOT SEEN BEFORE]
              if(cModels[n]->correct.in != 0)
                Fail(cModels[n]);
            break;
            case -1:  // IT HAS AT LEAST TWO MAXIMUM FREQS [CONFUSION FREQS]
              if(cModels[n]->correct.in != 0)
                Fail(cModels[n]);
            break;
            default:  // IT HAS ONE MAXIMUM FREQ
              if(cModels[n]->correct.in == 0){ // IF IS OUT
                cModels[n]->correct.in = 1;
                memset(cModels[n]->correct.mask, 0, cModels[n]->ctx);
                }
              else{ // IF IS IN
                if(best == sym) Hit(cModels[n]);
                else{
                  Fail(cModels[n]);
                  cModels[n]->correct.seq->buf[cModels[n]->correct.seq->idx]
                  = best; // UPDATE BUFFER WITH NEW SYMBOL
                  }
                }
            }
          UpdateCBuffer(cModels[n]->correct.seq);
          }
        }
      #endif

    if(++idxOut == BUFFER_SIZE){
      fwrite(outBuffer, 1, idxOut, Writter);
      idxOut = 0;
      }

    if(++idx == BUFFER_SIZE){
      memcpy(symbolBuffer-BGUARD, symbolBuffer+idx-BGUARD, BGUARD);
      idx = 0;
      }
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
  for(n = 0 ; n < P->nModels ; ++n){
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
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
  uint64_t  nBases = 0;
  int32_t   idx = 0;
  uint8_t   *readerBuffer, *symbolBuffer, sym, irSym, type = 0, header = 1,
            line = 0, dna = 0;
  CModel    **cModels;
  #ifdef PROGRESS
  uint64_t  i = 0;
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
      P->model[n].ir, REFERENCE, P->col, P->model[n].am);

  sym = fgetc(Reader);
  switch(sym){
    case '>': type = 1; break;
    case '@': type = 2; break;
    default : type = 0;
    }
  rewind(Reader);

  switch(type){
    case 1:  nBases = NDNASymInFasta(Reader); break;
    case 2:  nBases = NDNASymInFastq(Reader); break;
    default: nBases = NDNASyminFile (Reader); break;
    }

  P->checksum   = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {
      sym = readerBuffer[idxPos];
      if(type == 1){  // IS A FAST[A] FILE
        if(sym == '>'){ header = 1; continue; }
        if(sym == '\n' && header == 1){ header = 0; continue; }
        if(sym == '\n') continue;
        if(sym == 'N' ) continue;
        if(header == 1) continue;
        }
      else if(type == 2){ // IS A FAST[Q] FILE
        switch(line){
          case 0: if(sym == '\n'){ line = 1; dna = 1; } break;
          case 1: if(sym == '\n'){ line = 2; dna = 0; } break;
          case 2: if(sym == '\n'){ line = 3; dna = 0; } break;
          case 3: if(sym == '\n'){ line = 0; dna = 0; } break;
          }
        if(dna == 0 || sym == '\n') continue;
        if(dna == 1 && sym == 'N' ) continue;
        }

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
        continue;

      symbolBuffer[idx] = sym = DNASymToNum(sym);      
      P->checksum = (P->checksum + (uint8_t) sym);

      for(n = 0 ; n < P->nModels ; ++n)
        if(P->model[n].type == REFERENCE){
          GetPModelIdx(symbolBuffer+idx-1, cModels[n]);
          UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
          if(cModels[n]->ir == 1){                        // Inverted repeats
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            }
          }

      if(++idx == BUFFER_SIZE){
        memcpy(symbolBuffer - BGUARD, symbolBuffer + idx - BGUARD, BGUARD);
        idx = 0;
        }
      #ifdef PROGRESS
      CalcProgress(nBases, ++i);
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

int32_t main(int argc, char *argv[]){
  char        **p = *&argv;
  CModel      **refModels; 
  uint32_t    n, k, *checksum, refNModels = 0;
  Parameters  *P;
  FILE        *Reader = NULL;
  uint8_t     help, verbose, force, nTar = 1;
  
  if((help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    fprintf(stderr,
    "Usage: GeDe [OPTION]... -r [FILE]  [FILE]:[...]                    \n"
    "Decompress a genomic sequence compressed by GeCo.                  \n"
    "                                                                   \n"
    "Non-mandatory arguments:                                           \n"
    "                                                                   \n"
    "  -h                    give this help,                            \n"
    "  -v                    verbose mode (more information),           \n"
    "                                                                   \n"
    "  -r <FILE>             reference file,                            \n"
    "                                                                   \n"
    "Mandatory arguments:                                               \n"
    "                                                                   \n"
    "  <FILE>                file to uncompress (last argument). For    \n"
    "                        more files use splitting \":\" characters. \n"
    "                                                                   \n"
    "Report bugs to <{pratas,ap,pjf}@ua.pt>.                            \n");
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    fprintf(stderr, "GeDe from GeCo %u.%u\n"
    "Copyright (C) 2015 University of Aveiro.\nThis is Free software. \nYou "
    "may redistribute copies of it under the terms of the GNU General \n"
    "Public License v2 <http://www.gnu.org/licenses/gpl.html>.\nThere is NO "
    "WARRANTY, to the extent permitted by law.\nWritten by Diogo Pratas, "
    "Armando J. Pinho and Paulo J. S. G. Ferreira.\n", VERSION, RELEASE);
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
  for(n = 0 ; n < nTar ; ++n){
    Reader = Fopen(P[0].tar[n], "r");
    startinputtingbits();
    start_decode(Reader);

    refNModels = 0;
    P[n].watermark = ReadNBits(32, Reader);
    if(P[n].watermark != WATERMARK){
      fprintf(stderr, "Error: Invalid compressed file to decompress!\n");
      return 1;
      }
    checksum[n]    = ReadNBits(46, Reader);
    P[n].size      = ReadNBits(46, Reader);
    P[n].gamma     = ReadNBits(32, Reader) / 65536.0;
    P[n].col       = ReadNBits(32, Reader);
    P[n].nModels   = ReadNBits(16, Reader);
    P[n].model     = (ModelPar *) Calloc(P[n].nModels, sizeof(ModelPar));
    for(k = 0 ; k < P[n].nModels ; ++k){
      P[n].model[k].ctx  = ReadNBits(16, Reader); 
      P[n].model[k].den  = ReadNBits(16, Reader); 
      P[n].model[k].ir   = ReadNBits( 1, Reader); 
      P[n].model[k].am   = ReadNBits( 8, Reader); 
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
 
  if(refNModels > 0 && P[0].ref == NULL){
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

  for(n = 0 ; n < nTar ; ++n){
    if(refNModels != 0){
      if(CmpCheckSum(checksum[n], P[0].checksum) == 0)
        Decompress(P, refModels, n);
      }
    else
      Decompress(P, refModels, n);
    }

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
