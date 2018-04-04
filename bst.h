#define TRUE 1
#define FALSE 0

typedef struct bnode_ {
  struct bnode_* left;
  struct bnode_* right;
  int value;
  int size;
} bnode;


bnode* new_node(int value);
bnode* insert(bnode* root, int value);
int    lookup(bnode* root, int value);
int    size(bnode* root);
int*   put_in_array(bnode* root, int* array);
int*   get_array(bnode* root);
void   dispose(bnode* root);
void   view(bnode* root);
