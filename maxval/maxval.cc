#include <stdio.h>
#include <stdlib.h>

int main (int argc, char * argv[])
{
  int nbits = atoi (argv[1]);
  const long minval = -(1L << (nbits-1)) + 1;
  const long maxval = (1L << (nbits-1)) - 1;

  long v = 1;

  printf (" nbits = %d, minval = %20ld, maxval = %20ld, v = %ld, v < minval = %d, v > maxval = %d\n", nbits, minval, maxval, v, v < minval, v > maxval);

  return 0;
}
