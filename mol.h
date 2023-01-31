#ifndef _mol_h // only used once!
#define M_PI 3.141592653589793238462643383279502884
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define _mol_h
// atom stucture
typedef struct atom
{
    char element[3];
    double x, y, z;
} atom;
// bond structure
typedef struct bond
{
    atom *a1, *a2;
    unsigned char epairs;
} bond;
// molecule structure
typedef struct molecule
{
    unsigned short atom_max, atom_no; // atom max stores max array size and atom no current array size
    atom *atoms, **atom_ptrs;         // atom_ptrs is a pointer to a pointer (we can later change order of atoms)
    unsigned short bond_max, bond_no; // bond max stores max array size and bond no current array size
    bond *bonds, **bond_ptrs;
} molecule;

typedef double xform_matrix[3][3];
// funtions begin (explanations later)
void atomset(atom *atom, char element[3], double *x, double *y, double *z); // use strcpy to coppy char element

void atomget(atom *atom, char element[3], double *x, double *y, double *z);

void bondset(bond *bond, atom *a1, atom *a2, unsigned char epairs);

void bondget(bond *bond, atom **a1, atom **a2, unsigned char *epairs);

molecule *molmalloc(unsigned short atom_max, unsigned short bond_max); // set values in molecure and allocate memory

molecule *molcopy(molecule *src); // copy molecule

void molfree(molecule *ptr); // free memory stored by moleule

void molappend_atom(molecule *molecule1, atom *atom1);

void molappend_bond(molecule *molecule1, bond *bond1);

void molsort(molecule *molecule);

void xrotation(xform_matrix xform_matrix, unsigned short deg);

void yrotation(xform_matrix xform_matrix, unsigned short deg);

void zrotation(xform_matrix xform_matrix, unsigned short deg);

void mol_xform(molecule *molecule, xform_matrix matrix);

#endif