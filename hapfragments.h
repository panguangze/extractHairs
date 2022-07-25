
#ifndef _HAPFRAGMENT_H
#define _HAPFRAGMENT_H
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include "readvariant.h"

extern int SINGLEREADS;
extern int DATA_TYPE;
extern int NEW_FORMAT;
extern int PRINT_COMPACT; // default = 1

typedef struct {
    char allele;
    char qv;
    int varid; // allele is 0/1    varid is index to varlist[varid] gives all information about the variant
} allele;

typedef struct {
    char* id;
    int variants;
    allele* alist;
	int read_qual;
    int blocks;
    int paired;
    int matepos;
    int absIS;
    char* barcode;
    char strand; // added 01/26/2018
    int rescued;
    float dm;
} FRAGMENT;

// when --vcf-phased specified, filter fragment by phasing info
int filter_by_phasing_info(FRAGMENT* fragment, VARIANT* varlist, std::map<int, int>* homo_recom);

int compare_fragments(const void *a, const void *b);

int compare_alleles(const void *a, const void *b);

int print_fragment(FRAGMENT* fragment, VARIANT* varlist, FILE* outfile, std::map<int, int>* homo_recom);

// make sure they are in the correct order, i+1 could be < i
int print_matepair(FRAGMENT* f1, FRAGMENT* f2, VARIANT* varlist, FILE* outfile);

// add mate bnd connection
int print_mate_bnd_fragment(std::unordered_map<std::string , std::pair<int, int>>& BNDs, FILE* outfile);
void clean_fragmentlist(FRAGMENT* flist, int* fragments, VARIANT* varlist, int currchrom, int currpos, int prevchrom);

#endif
