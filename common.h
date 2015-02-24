#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include "context.h"
#include "defs.h"

typedef struct{
  uint64_t size;
  uint64_t bytes;
  }
INF;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

static const uint32_t    DEFAULT_CONTEXT[]  = { 2, 20, 2, 14 };
static const uint32_t    DEFAULT_ALPHADEN[] = { 1, 50, 1, 10 };       
static const uint32_t    DEFAULT_IR[]       = { 0,  1, 0,  1 };
static const uint32_t    DEFAULT_AM[]       = { 0,  4, 0,  3 };

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILE        *Fopen           (const char *, const char *);
void        ShiftBuffer      (uint8_t *, int, uint8_t);
void        FillLogTable     (uint32_t, uint32_t, uint32_t);
double      SearchLog        (uint32_t );
double      Power            (double, double);
uint32_t    FLog2            (uint64_t );
char        *ReplaceSubStr   (char *, char *, char *);
uint8_t     DNASymToNum      (uint8_t  );
uint8_t     NumToDNASym      (uint8_t  );
uint8_t     GetCompNum       (uint8_t  );
uint8_t     GetCompSym       (uint8_t  );
uint64_t    NDNASyminFile    (FILE *);
uint64_t    NDNASymInFastq   (FILE *);
uint64_t    NDNASymInFasta   (FILE *);
uint64_t    NBytesInFile     (FILE *);
uint64_t    FopenBytesInFile (const char *);
uint8_t     *ReverseStr      (uint8_t *, uint32_t);
char        *CloneString     (char *   );
char        *RepString       (const char *, const char *, const char *);
uint32_t    ArgsNum          (uint32_t , char *[], uint32_t, char *, uint32_t,
                              uint32_t);
ModelPar    ArgsUniqModel    (char *, uint8_t);
ModelPar    ArgsModel        (uint32_t , char *[], uint32_t, char *);
double      ArgsDouble       (double, char *[], uint32_t, char *);
uint8_t     ArgsState        (uint8_t  , char *[], uint32_t, char *);
char        *ArgsString      (char    *, char *[], uint32_t, char *);
char        *ArgsFiles       (char *[], uint32_t, char *);
void        TestReadFile     (char *);
uint8_t     CmpCheckSum      (uint32_t, uint32_t);
void        FAccessWPerm     (char    *);
uint32_t    ReadFNames       (Parameters *, char *);
inline void CalcProgress     (uint64_t , uint64_t);
void        PrintArgs        (Parameters *);
char        *concatenate     (char *   , char *);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
