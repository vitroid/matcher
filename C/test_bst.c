#include <stdio.h>
#include <stdlib.h>
#include "bst.h"

int main()
{
  bnode* root = NULL;
  root = insert(root,5);
  root = insert(root,4);
  root = insert(root,8);
  root = insert(root,2);
  root = insert(root,7);
  root = insert(root,10);
  root = insert(root,9);
  root = insert(root,11);

  view(root);
  int* array = get_array(root);
  for(int i=0;i<size(root); i++){
    printf("%d %d\n", i, array[i]);
  }
  for(int i=0;i<20;i++){
    printf("%d: %d\n", i, lookup(root,i));
  }
  dispose(root);
}
  
  
