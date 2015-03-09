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
#include "levels.h"
#include "common.h"
#include "context.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void Compress(Parameters *P, CModel **cModels, uint8_t id, uint32_t 
refNModels, INF *I){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".co");
  FILE        *Writter = Fopen(name, "w");
  uint32_t    n, k, idxPos;
  int32_t     idx = 0;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0;
  double       *cModelWeight, cModelTotalWeight = 0, mA, mC, mG, mT;
  uint8_t     *readerBuffer, *symbolBuffer, sym, irSym, *pos, type = 0, 
              header = 1, line = 0, dna = 0;
  PModel      **pModel, *MX;
  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stderr, "Analyzing data and creating models ...\n");

  #ifdef ESTIMATE
  FILE *IAE = NULL;
  char *IAEName = NULL;
  if(P->estim == 1){
    IAEName = concatenate(P->tar[id], ".iae");
    IAE = Fopen(IAEName, "w");
    }
  #endif
  
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
 
  _bytes_output = 0;
  nSymbols      = NBytesInFile(Reader);
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
    if(P->model[n].type == TARGET){
      cModels[n] = CreateCModel(P->model[n].ctx, P->model[n].den, 
      P->model[n].ir, TARGET, P->col, P->model[n].am);
      }
    }

  if(P->verbose){
    fprintf(stderr, "Done!\n");
    fprintf(stderr, "Compressing target sequence %d [bases: %"PRIu64"] ...\n", 
    id + 1, nBases);
    }

  startoutputtingbits();
  start_encode();

  WriteNBits(WATERMARK,                32, Writter);
  WriteNBits(P->checksum,              46, Writter);
  WriteNBits(nBases,                   46, Writter);
  WriteNBits((int) (P->gamma * 65536), 32, Writter);
  WriteNBits(P->col,                   32, Writter);
  WriteNBits(P->nModels,               16, Writter);
  for(n = 0 ; n < P->nModels ; ++n){
    WriteNBits(cModels[n]->ctx,        16, Writter);
    WriteNBits(cModels[n]->alphaDen,   16, Writter);
    WriteNBits(cModels[n]->ir,          1, Writter);
    WriteNBits(cModels[n]->am,          8, Writter);
    WriteNBits(P->model[n].type,        1, Writter);
    }

  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){
      #ifdef PROGRESS
      CalcProgress(nSymbols, ++i);
      #endif

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
      mA = mC = mG = mT = 0;

      pos = &symbolBuffer[idx-1];
      for(n = 0 ; n < P->nModels ; ++n){
        if(cModels[n]->am == 0){
          GetPModelIdx(pos, cModels[n]);
          }
        else{
          cModels[n]->correct.seq->buf[cModels[n]->correct.seq->idx] = sym;
          GetPModelIdx(cModels[n]->correct.seq->buf+cModels[n]->correct.seq->idx
          -1, cModels[n]);
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

      AESym(sym, (int *)(MX->freqs), (int) MX->sum, Writter);
      #ifdef ESTIMATE
      if(P->estim != 0)
        fprintf(IAE, "%.3g\n", PModelSymbolNats(MX, sym));
      #endif

      cModelTotalWeight = 0;
      for(n = 0 ; n < P->nModels ; ++n){
        cModelWeight[n] = Power(cModelWeight[n], P->gamma) * (double) 
        pModel[n]->freqs[sym] / pModel[n]->sum;
        cModelTotalWeight += cModelWeight[n];
        if(cModels[n]->ref == TARGET){
          if(cModels[n]->am != 0)
            UpdateCModelCounter(cModels[n], sym, cModels[n]->correct.idx);
          else
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);

          if(cModels[n]->ir != 0){                // REVERSE COMPLEMENTS
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            }
          }
        }

      for(n = 0 ; n < P->nModels ; ++n){
        cModelWeight[n] /= cModelTotalWeight; // RENORMALIZE THE WEIGHTS
        if(cModels[n]->am != 0)
          CorrectCModel(cModels[n], pModel[n], sym);
        }

      if(++idx == BUFFER_SIZE){
        memcpy(symbolBuffer-BGUARD, symbolBuffer+idx-BGUARD, BGUARD);
        idx = 0;
        }

      ++compressed;
      }

  finish_encode(Writter);
  doneoutputtingbits(Writter);
  fclose(Writter);

  #ifdef ESTIMATE
  if(P->estim == 1){
    fclose(IAE);
    Free(IAEName);
    }
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
  I[id].size  = compressed;
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
  cModels       = (CModel **) Malloc(P->nModels * sizeof(CModel *)); 
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

  P->checksum = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

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
          if(cModels[n]->ir == 1){                         // Inverted repeats
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
  char        **p = *&argv, **xargv, *xpl = NULL;
  CModel      **refModels;
  int32_t     xargc = 0;
  uint32_t    n, k, refNModels, col;
  uint64_t    totalBytes, totalSize;
  double      gamma;
  
  Parameters  *P;
  INF         *I;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h")) == 1 || argc < 2){
    fprintf(stderr,
    "Usage: GeCo [OPTION]... -r [FILE]  [FILE]:[...]                        \n"
    "Compress and analyze a genomic sequence (by default, compress).        \n"
    "                                                                       \n"
    "Non-mandatory arguments:                                               \n"
    "                                                                       \n"
    "  -h                    give this help,                                \n"
    "  -x                    show several running examples,                 \n"
    "  -s                    show GeCo compression levels,                  \n"
    "  -v                    verbose mode (more information),               \n"
    "  -V                    display version number,                        \n"
    "  -f                    force overwrite of output,                     \n"
    "  -l <level>            level of compression [1;9] (lazy -tm setup),   \n"
    "  -g <gamma>            mixture decayment forgetting factor. It is     \n"
    "                        a real value in the interval [0;1),            \n"
    "  -c <cache>            maximum collisions for hash cache. Memory      \n"
    "                        values are higly dependent of the parameter    \n"
    "                        specification,                                 \n"
    #ifdef ESTIMATE
    "  -e                    it creates a file with the extension \".iae\"  \n" 
    "                        with the respective information content. If    \n" 
    "                        the file is FASTA or FASTQ it will only use    \n"
    "                        the \"ACGT\" (genomic) data,                   \n"
    #endif
    "  -r <FILE>             reference file (\"-rm\" are loaded here),      \n"
    "                                                                       \n"
    "Mandatory arguments:                                                   \n"
    "                                                                       \n"
    "  -rm <c>:<d>:<i>:<m>   reference context model (ex:-rm 13:100:0:0),   \n"
    "  -rm <c>:<d>:<i>:<m>   reference context model (ex:-rm 18:1000:0:1),  \n"
    "  ...                                                                  \n"
    "  -tm <c>:<d>:<i>:<m>   target context model (ex:-tm 4:1:0:0),         \n"
    "  -tm <c>:<d>:<i>:<m>   target context model (ex:-tm 18:20:1:1),       \n"
    "  ...                                                                  \n"
    "                        target and reference templates use <c> for     \n"
    "                        context-order size, <d> for alpha (1/<d>),     \n"
    "                        <i> (0 or 1) to set the usage of inverted      \n"
    "                        repeats (1 to use) and <m> to the maximum      \n"
    "                        allowed mutation on the context without        \n"
    "                        being discarded (usefull in deep contexts),    \n"
    "                                                                       \n"
    "  <FILE>                file to compress (last argument). For more     \n"
    "                        files use splitting \":\" characters.          \n"
    "                                                                       \n"
    "Report bugs to <{pratas,ap,pjf}@ua.pt>.                              \n");
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V")){
    fprintf(stderr, 
    "                                                                       \n"
    "                            ============                               \n"
    "                            | GeCo %u.%u |                             \n"
    "                            ============                               \n"
    "                                                                       \n"
    "[ A compressor and analyzer for multiple genomic [A,C,G,T] sequences ] \n"
    "Copyright (C) 2014-2015 University of Aveiro. This is a Free software. \n"
    "You may redistribute copies of it under the terms of the GNU - General \n"
    "Public License v2 <http://www.gnu.org/licenses/gpl.html>. There is NOT \n"
    "ANY WARRANTY, to the extent permitted by law. Developed and Written by \n"
    "Diogo Pratas, Armando J. Pinho and Paulo J. S. G. Ferreira.\n\n", VERSION, 
    RELEASE);
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s")){
    PrintLevels(); 
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_EXAMPLE, p, argc, "-x")){
    fprintf(stderr, 
    "                                                                       \n"
    "GeCo running examples:                                                 \n"
    "                                                                       \n"
    "Considerations: the decompression is symmetric, therefore the same     \n"
    "resources, namely time and memory will be used as in the compression.  \n"
    "The memory used, after creating the models, will be constant, even in  \n"
    "deeper context models (cache-hash context model).                      \n"
    "                                                                       \n"
    "[A]=> Compressing sequences C(X) or C(X,Y):                            \n"
    "                                                                       \n"
    "1) Compression of a human genome (using 5.8 GB RAM memory):            \n"
    "   ./GeCo -tm 6:1:0:0 -tm 13:20:1:0 -tm 19:50:1:2 -c 35 -g 0.8 HS      \n"          "                                                                       \n"
    "2) Compression of a human genome (using 3.8 GB of RAM memory):         \n"
    "   ./GeCo -tm 6:1:0:0 -tm 13:10:1:0 -tm 18:50:1:0 -c 20 -g 0.8 HS      \n" 
    "                                                                       \n"
    "3) Compression of a human genome (using 0.5 GB of RAM memory):         \n"
    "   ./GeCo -tm 6:1:0:0 -tm 13:10:1:0 -g 0.8 HS                          \n" 
    "                                                                       \n" 
    "   Decompression for A1, A2 and A3: ./GeDe HS.co                       \n"
    "   The decompressed file will be HS.de                                 \n"
    "                                                                       \n"
    "4) Compression of a human chromosome Y (repetitive nature):            \n"
    "   ./GeCo -tm 1:1:0:0 -tm 4:1:0:0 -tm 6:1:1:0 -tm 8:1:0:0 -tm 11:10:1:0\n"
    "   -tm 14:10:0:1 -tm 14:50:1:0 -tm 18:30:1:6 -c 10 -g 0.88 CY.fasta    \n"
    "   Decompression for A4: ./GeDe CY.fasta.co                            \n"
    "   The decompressed file will be CY.fasta.de                           \n"
    "                                                                       \n"
    "5) Highly-redundant genomic sequence (full ACGT from fastq)            \n"
    "   ./GeCo -tm 4:1:0:0 -tm 11:1:0:0 -tm 14:20:0:0 -tm 20:100:0:0 -c 50  \n"
    "   -g 0.8 SRR957627.fastq                                              \n"
    "   Decompression for A5: ./GeDe SRR957627.fastq.co                     \n"
    "   The decompressed file will be SRR957627.fastq.de                    \n"
    "                                                                       \n"
    "                                                                       \n"
    "[B]=> Conditional (referential) exclusive compression C(X||Y):         \n"
    "                                                                       \n"
    "1) Compression of the gorilla (GG8) chromosome 8 given exclusively     \n"
    "   information from chimpanzee (PT8):                                  \n"
    "   ./GeCo -rm 4:1:0:0 -rm 18:100:1:0 -c 20 -r PT8 GG8                  \n"
    "   Decompression for B1: ./GeDe -r PT8 GG8.co                          \n"
    "   The decompressed file will be GG8.de                                \n"
    "                                                                       \n"
    "2) Compression of the same file (for idempotency studies):             \n"
    "   ./GeCo -rm 20:1000:0:0 -c 30 -r File1.txt File1.txt                 \n"
    "   Decompression for B2: ./GeDe -r File1.txt File1.txt.co              \n"
    "   The decompressed file will be File1.txt.de                          \n"
    "                                                                       \n"
    "3) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  \n"
    "   chromsomes 5 given exclusively the gorilla (GG17) chromosome 17 as  \n"
    "   reference:                                                          \n"
    "   ./GeCo -rm 20:100:1:0 -c 20 -r GG17 HS5:PT5:PA5                     \n"
    "   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           \n"
    "   The decompressed files will be HS5.de, PT5.de and PA5.de            \n"
    "                                                                       \n"
    "                                                                       \n"
    "[C]=> Conditional compression C(X|Y) [use reference and target]:       \n"
    "                                                                       \n"
    "1) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  \n"
    "   chromsomes 5 given the gorilla (GG17) chromosome 17 as reference:   \n"
    "   ./GeCo -tm 4:1:0:0 -tm 8:1:0:0 -tm 11:10:0:0 -tm 14:20:1:0          \n"
    "   -rm 20:100:0:0 -c 20 -r GG17 HS5:PT5:PA5                            \n"
    "   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           \n"
    "   The decompressed files will be HS5.de, PT5.de and PA5.de            \n"
    "                                                                     \n");
    return EXIT_SUCCESS;
    }

  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v" );
  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-f" );
  P->estim    = ArgsState  (0,               p, argc, "-e" );
  P->level    = ArgsNum    (0, p, argc, "-l", MIN_LEVEL, MAX_LEVEL);

  P->nModels  = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0 || strcmp(argv[n], "-tm") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;
  
  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-rm") == 0 || strcmp(xargv[n], "-tm") == 0)
        P->nModels += 1;
    }

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
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-rm") == 0){
        P->model[k++] = ArgsUniqModel(xargv[n+1], 1);
        ++refNModels;
        }
    }

  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-tm") == 0)
      P->model[k++] = ArgsUniqModel(argv[n+1], 0);
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-tm") == 0)
        P->model[k++] = ArgsUniqModel(xargv[n+1], 0);
    }

  gamma = DEFAULT_GAMMA;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-g") == 0) 
      gamma = atof(xargv[n+1]);

  col = MAX_COLLISIONS;
  for(n = 1 ; n < xargc ; ++n) 
    if(strcmp(xargv[n], "-c") == 0) 
      col = atoi(xargv[n+1]);

  P->col      = ArgsNum    (col,   p, argc, "-c", 1, 10000);
  P->gamma    = ArgsDouble (gamma, p, argc, "-g");
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
  fprintf(stderr, "), %.4g bpb, Distance: %.6g\n", ((8.0*totalBytes)/
  totalSize), (4.0*totalBytes)/totalSize);  

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
