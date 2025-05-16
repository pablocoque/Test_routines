#include <stdio.h>
#include <math.h>
#include <assert.h>

void sum_with_pointer(int *var){
  *var += 1;
}

int main(){
  int var = 0;
  if(!var) printf("SIRVE ESTA VAINA \n");
  for (int i=0;i<3;i++){
    sum_with_pointer(&var);
  }
  
  printf("called %d times \n", var);

  double try = 0.0;
  if(!(try > 0)) printf("Try var is zero: %f \n", try);

  try += 100.0;
  if(!(try > 0)) printf("Should not happen: %f \n", try);

  try = 0.0/0.0;
  printf("Last try %f \n", try);

  assert(try <0.0);


}
