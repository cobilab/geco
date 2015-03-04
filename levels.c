#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "levels.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET PARAMETERS FROM LEVELS
//
char *GetLevels(uint8_t l){
  switch(l){
    case 1: return LEVEL_1;
    case 2: return LEVEL_2;
    case 3: return LEVEL_3;
    case 4: return LEVEL_4;
    case 5: return LEVEL_5;
    case 6: return LEVEL_6;
    case 7: return LEVEL_7;
    case 8: return LEVEL_8;
    case 9: return LEVEL_9;
    default: fprintf(stderr, "[x] Unknown level!\n");
    exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// PRINT LEVELS FROM LEVELS
//
void PrintLevels(void){
  fprintf(stderr, "Level %u: %s\n", 1, LEVEL_1);
  fprintf(stderr, "Level %u: %s\n", 2, LEVEL_2);
  fprintf(stderr, "Level %u: %s\n", 3, LEVEL_3);
  fprintf(stderr, "Level %u: %s\n", 4, LEVEL_4);
  fprintf(stderr, "Level %u: %s\n", 5, LEVEL_5);
  fprintf(stderr, "Level %u: %s\n", 6, LEVEL_6);
  fprintf(stderr, "Level %u: %s\n", 7, LEVEL_7);
  fprintf(stderr, "Level %u: %s\n", 8, LEVEL_8);
  fprintf(stderr, "Level %u: %s\n", 9, LEVEL_9);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

