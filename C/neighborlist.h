#ifndef _PAIRLIST_NEIGHBORLIST_H
#define _PAIRLIST_NEIGHBORLIST_H

#define ADDRESS(x,y,z) (((z)*GY + (y))*GX + (x))
#define True 1
#define False 0


typedef struct
{
  int grid[3];
  int* nresidents;
  int** residents;//allocated list for residents in each district
}
  sAddressBook;


sAddressBook* AddressBook(int grid[3], int npos, double* rpos);
int
near_origin(sAddressBook* abook, int** residents);
void
dispose_addressbook(sAddressBook* abook);
int find_nearest(double rloc[3], sAddressBook* abook, double cell[3], double rpos[]);


#endif
