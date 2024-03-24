#include <algorithm>
#include <stdio.h>
#include <vector>


int main (int argc, char * argv[])
{

#pragma omp parallel
{
  struct tt
  {
    int a;
    int b;
  };

  std::vector<tt> vec = 
  {
     { 0, 52366 }, { 0, 52366 }, { 0, 52367 }, { 0, 52367 }, { 0, 52368 },
     { 0, 52368 }, { 0, 52368 }, { 0, 52369 }, { 0, 52369 }, { 0, 52369 },
     { 0, 52370 }, { 0, 52370 }, { 0, 52371 }, { 0, 52371 }, { 0, 52371 },
     { 0, 52372 }, { 0, 52372 }, { 0, 52372 }, { 0, 52373 }, { 0, 52373 },
     { 0, 52373 }, { 0, 52374 }, { 0, 52374 }, { 0, 52375 }, { 0, 52375 },
     { 0, 52375 }, { 0, 52376 }, { 0, 52376 }, { 0, 52377 }, { 0, 52377 },
     { 0, 52377 }, { 0, 52378 }, { 0, 52378 }, { 0, 52378 }, { 0, 52379 },
     { 0, 52379 }, { 0, 52380 }, { 0, 52380 }, { 0, 52380 }, { 0, 52381 },
     { 0, 52381 }, { 0, 52381 }, { 0, 52382 }, { 0, 52382 }, { 0, 52383 },
     { 0, 52383 }, { 0, 52383 }, { 0, 52384 }, { 0, 52384 }, { 0, 52384 },
  };

  std::vector<int> ord (vec.size ());

  std::iota (std::begin (ord), std::end (ord), 0);

  std::stable_sort (std::begin (ord), std::end (ord), 
            [&vec] (int i, int j) 
            { 
              printf("compare %8d and %8d\n", i, j);
              if (vec[i].a == vec[j].a)
                return vec[i].b < vec[j].b;
              return vec[i].a < vec[j].a;
            });

  for (int k = 0; k < vec.size (); k++)
    printf (" %8d %8d\n", k, vec[ord[k]].b);

}

  return 0;
}


