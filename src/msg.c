#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "                                                                       \n"
  "  -rm <c>:<d>:<i>:<m/e>  reference context model (ex:-rm 13:100:0:0/0), \n"
  "  -rm <c>:<d>:<i>:<m/e>  reference context model (ex:-rm 18:1000:0:1/1000),\n"
  "  ...                                                                  \n"
  "  -tm <c>:<d>:<i>:<m/e>  target context model (ex:-tm 4:1:0:0/0),      \n"
  "  -tm <c>:<d>:<i>:<m/e>  target context model (ex:-tm 18:20:1:2/10),   \n"
  "  ...                                                                  \n"
  "                         target and reference templates use <c> for    \n"
  "                         context-order size, <d> for alpha (1/<d>),    \n"
  "                         <i> (0 or 1) to set the usage of inverted     \n"
  "                         repeats (1 to use) and <m> to the maximum     \n"
  "                         allowed mutation on the context without       \n"
  "                         being discarded (usefull in deep contexts),   \n"
  "                         under the estimator <e>,                      \n");
  } 

void PrintMenu(void){
  fprintf(stderr,
  "Usage: GeCo [OPTION]... -r [FILE]  [FILE]:[...]                        \n"
  "Compress and analyze a genomic sequence (by default, compress).        \n"
  "                                                                       \n"
  "Non-mandatory arguments:                                               \n"
  "                                                                       \n"
  "  -h                     give this help,                               \n"
  "  -x                     show several running examples,                \n"
  "  -s                     show GeCo compression levels,                 \n"
  "  -v                     verbose mode (more information),              \n"
  "  -V                     display version number,                       \n"
  "  -f                     force overwrite of output,                    \n"
  "  -l <level>             level of compression [1;9] (lazy -tm setup),  \n"
  "  -g <gamma>             mixture decayment forgetting factor. It is    \n"
  "                         a real value in the interval [0;1),           \n"
  "  -c <cache>             maximum collisions for hash cache. Memory     \n"
  "                         values are higly dependent of the parameter   \n"
  "                         specification,                                \n");
  #ifdef ESTIMATE
  fprintf(stderr,
  "  -e                     it creates a file with the extension \".iae\" \n"
  "                         with the respective information content. If   \n"
  "                         the file is FASTA or FASTQ it will only use   \n"
  "                         the \"ACGT\" (genomic) data,                  \n");
  #endif
  ModelsExplanation();
  fprintf(stderr,
  "                                                                       \n"
  "  -r <FILE>              reference file (\"-rm\" are loaded here),     \n"
  "                                                                       \n"
  "Mandatory arguments:                                                   \n"
  "                                                                       \n"
  "  <FILE>                 file to compress (last argument). For more    \n"
  "                         files use splitting \":\" characters.         \n"
  "                                                                       \n"
  "Report bugs to <{pratas,ap,pjf}@ua.pt>.                              \n");
  }


void PrintVersion(void){
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
  }


void PrintExamples(void){
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
  "   ./GeCo -tm 6:1:0:0/0 -tm 13:20:1:0/0 -tm 19:50:1:2/10 -c 35 -g 0.8 HS\n"
  "                                                                       \n"
  "2) Compression of a human genome (using 3.8 GB of RAM memory):         \n"
  "   ./GeCo -tm 6:1:0:0/0 -tm 13:10:1:0/0 -tm 18:50:1:3/10 -c 20 -g 0.8 HS\n"
  "                                                                       \n"
  "3) Compression of a human genome (using 0.5 GB of RAM memory):         \n"
  "   ./GeCo -tm 6:1:0:0/0 -tm 13:10:1:0/0 -g 0.8 HS                      \n"
  "                                                                       \n"
  "   Decompression for A1, A2 and A3: ./GeDe HS.co                       \n"
  "   The decompressed file will be HS.de                                 \n"
  "                                                                       \n"
  "4) Compression of a human chromosome Y (repetitive nature):            \n"
  "   ./GeCo -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:1:0/0 -tm 8:1:0:0/0      \n"
  "   -tm 11:10:1:0/0 -tm 14:10:0:1/10 -tm 14:50:1:0/0 -tm 18:30:1:6/10   \n"
  "   -c 10 -g 0.88 CY.fasta                                              \n"
  "   Decompression for A4: ./GeDe CY.fasta.co                            \n"
  "   The decompressed file will be CY.fasta.de                           \n"
  "                                                                       \n"
  "5) Highly-redundant genomic sequence (full ACGT from fastq)            \n"
  "   ./GeCo -tm 4:1:0:0/0 -tm 11:1:0:0/0 -tm 14:20:0:0/0 -tm 20:100:0:1/10\n"
  "   -c 40 -g 0.8 SRR957627.fastq                                        \n"
  "   Decompression for A5: ./GeDe SRR957627.fastq.co                     \n"
  "   The decompressed file will be SRR957627.fastq.de                    \n"
  "                                                                       \n"
  "                                                                       \n"
  "[B]=> Conditional (referential) exclusive compression C(X||Y):         \n"
  "                                                                       \n"
  "1) Compression of the gorilla (GG8) chromosome 8 given exclusively     \n"
  "   information from chimpanzee (PT8):                                  \n"
  "   ./GeCo -rm 4:1:0:0/0 -rm 20:1000:1:1/100 -c 20 -r PT8 GG8           \n"
  "   Decompression for B1: ./GeDe -r PT8 GG8.co                          \n"
  "   The decompressed file will be GG8.de                                \n"
  "                                                                       \n"
  "2) Compression of the same file (for identity studies):                \n"
  "   ./GeCo -rm 20:1000:0:0/0 -c 30 -r File1.txt File1.txt               \n"
  "   Decompression for B2: ./GeDe -r File1.txt File1.txt.co              \n"
  "   The decompressed file will be File1.txt.de                          \n"
  "                                                                       \n"
  "3) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  \n"
  "   chromsomes 5 given exclusively the gorilla (GG17) chromosome 17 as  \n"
  "   reference:                                                          \n"
  "   ./GeCo -rm 20:1000:1:1/100 -c 20 -r GG17 HS5:PT5:PA5                \n"
  "   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           \n"
  "   The decompressed files will be HS5.de, PT5.de and PA5.de            \n"
  "                                                                       \n"
  "                                                                       \n"
  "[C]=> Conditional compression C(X|Y) [use reference and target]:       \n"
  "                                                                       \n"
  "1) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  \n"
  "   chromsomes 5 given the gorilla (GG17) chromosome 17 as reference:   \n"
  "   -rm 12:100:1:0/0 -rm 20:1000:1:1/100 -tm 4:1:0:0/0 -tm 14:20:1:1/10 \n"
  "   -c 20 -g 0.85 -r GG17 HS5:PT5:PA5                                   \n"
  "   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           \n"
  "   The decompressed files will be HS5.de, PT5.de and PA5.de            \n"
  "                                                                     \n");
  }

