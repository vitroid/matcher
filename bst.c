#include <stdio.h>
#include <stdlib.h>
#include "bst.h"
//binary search tree as an implementation of Set type.

bnode* new_node(int value){
  bnode* node = (bnode*) malloc(sizeof(bnode));
  node->right = node->left = NULL;
  node->value = value;
  node->size  = 1;
  return node;
}


int size(bnode* root)
{
  if (root==NULL){
    return 0;
  }
  return root->size;
}


bnode*
insert(bnode* root, int value)
{
  if (root == NULL){
    return new_node(value);
  }
  if (root->value == value){
    //no update
    return root;
  }
  if (root->value < value){
    int orig = size(root->right);
    root->right = insert(root->right, value);
    int diff = root->right->size - orig;
    root->size += diff;
    return root;
  }
  int orig = size(root->left);
  root->left = insert(root->left, value);
  int diff = root->left->size - orig;
  root->size += diff;
  return root;
}



int
lookup(bnode* root, int value)
{
  if (root == NULL){
    return FALSE;
  }
  if (root->value == value){
    return TRUE;
  }
  if (root->value < value){
    return lookup(root->right, value);
  }
  return lookup(root->left, value);
}



int*
put_in_array(bnode* root, int* array)
{
  if (root == NULL){
    return array;
  }
  array = put_in_array(root->left, array);
  array[0] =root->value;
  array++;
  return put_in_array(root->right, array);
}


int*
get_array(bnode* root){
  int* array = (int*) malloc(sizeof(int)*(size(root)+1)); //one extra space.
  put_in_array(root, array);
  array[size(root)] = -1; //terminator
  return array;
}


void
dispose(bnode* root)
{
  if (root==NULL){
    return;
  }
  dispose(root->left);
  dispose(root->right);
  free(root);
}



void view(bnode* root){
  if ( root == NULL )
    return;
  printf("%d (size %d) [", root->value, size(root));
  if (root->left == NULL){
    printf("- , ");
  }
  else{
    printf("%d, ", root->left->value);
  }
  if (root->right == NULL){
    printf("- ]\n");
  }
  else{
    printf("%d]\n", root->right->value);
  }
  view(root->left);
  view(root->right);
}

