#ifndef LEVELS_H_INCLUDED
#define LEVELS_H_INCLUDED

#include "defs.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// COMPRESSION LEVELS FOR DNA 
//
#define LEVEL_1 " 1: -rm 20:1000:1:3/100 -g 0.95 "
#define LEVEL_2 " 2: -tm 4:1:0:0/0 -tm 11:10:0:0/0 -tm 13:20:1:0/0 -tm 19:20:1:2/10 -c 30 -g 0.85 "
// HUMAN ASSEMBLED CY -> 1.24 using 2.5 GB of memory
#define LEVEL_3 " 3: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 8:1:0:0/0 -tm 12:10:0:0/0 -tm 13:20:1:1/10 -tm 18:100:1:6/10 -c 9 -g 0.90 "
//#define LEVEL_3 " 3: -tm 1:1:0:0 -tm 4:1:0:0 -tm 6:1:0:0 -tm 8:1:0:0 -tm 13:10:1:0 -tm 14:10:1:1 -tm 18:10:1:6 -tm 18:100:1:0 -c 9 -g 0.9 "
#define LEVEL_4 " 4: -tm 4:1:0:0/0 -tm 6:1:1:0/0 -tm 13:20:1:0/0 -tm 18:20:1:3/10 -c 20 -g 0.9 "

// HUMAN ASSEMBLED GENOME (GRC) -> 4.8 G : 52m : 1.585 bpb
#define LEVEL_5 " 5: -tm 4:1:0:0/0 -tm 11:1:0:0/0 -tm 13:20:1:0/0 -tm 19:20:1:3/10 -c 30 -g 0.85 "
#define LEVEL_6 " 6: -tm 1:1:0:0 -tm 6:1:0:0 -tm 11:10:0:0 -tm 13:20:1:0 -tm 18:50:1:3 -g 0.9 -c 30 "
#define LEVEL_7 " 7: -tm 1:1:0:0 -tm 4:1:0:0 -tm 6:1:0:0 -tm 9:1:0:0 -tm 11:10:0:0 -tm 14:20:1:0 -tm 19:50:1:3 -g 0.9 -c 30"
#define LEVEL_8 " 8: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 12:10:1:0/0 -tm 14:50:1:0/0 -tm 19:100:1:3/10 -c 30 -g 0.9 "
#define LEVEL_9 " 9: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 12:10:1:0/0 -tm 14:50:1:0/0 -tm 19:100:1:3/10 -c 35 -g 0.9 "
#define LEVEL_10 " 10: -tm 1:1:0:0/0 -tm 4:1:0:0/0 -tm 6:1:0:0/0 -tm 9:1:0:0/0 -tm 11:1:0:0/0 -tm 13:10:0:0/0 -tm 14:20:1:0/0 -tm 19:20:1:3/10 -c 50 -g 0.8 "

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char    *GetLevels  (uint8_t);
void    PrintLevels (void);

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

