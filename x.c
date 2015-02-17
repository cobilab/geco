#include <stdio.h>

int main()
  {
  unsigned long long a, b, ref = ((unsigned long long) 2<<40)+5;

  printf("ref=%llu\n-------------\n", ref);

  a = ref & 0xffffffffff;
  b = ref & 0xf0000000000;

  printf("a=%llu\n", a);
  printf("b=%llu\n", b);
  printf("a|b=%llu\n", a|b);

  return 0;
  }
