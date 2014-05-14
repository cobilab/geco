#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

typedef unsigned long long ULL;

typedef struct
  {
  uint32_t ctx;
  uint32_t den;
  uint32_t ir;
  uint8_t  type;
  }
ModelPar;

typedef struct
  {
  uint8_t  help;
  uint8_t  verbose;
  uint8_t  force;
  ModelPar *model;
  char     *ref;
  char     **tar;
  uint8_t  nTar;
  uint64_t checksum;
  uint64_t size;
  uint32_t watermark;
  double   gamma;
  uint32_t nModels;
  }
Parameters;

uint32_t garbage;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define BUFFER_SIZE            65535      
#define PROGRESS_MIN           1000000
#define DEFAULT_HELP           0
#define DEFAULT_VERBOSE        0
#define DEFAULT_FORCE          0
#define MAX_CTX                31
#define MIN_CTX                1
#define MAX_DEN                1000000
#define MIN_DEN                1
#define BGUARD                 32
#define DEFAULT_MAX_COUNT      ((1 << (sizeof(ACCounter) * 8)) - 1)
#define MX_PMODEL              65535
#define ALPHABET_SIZE          4
#define CHECKSUMGF             1073741824
#define WATERMARK              16042014
#define DEFAULT_GAMMA          0.95
#define MAX_HISTORYSIZE        1000000
#define REFERENCE              1
#define TARGET                 0

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

