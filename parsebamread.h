#ifndef __PARSEBAMREAD_H__
#define __PARSEBAMREAD_H__

#include "hapfragments.h"
#include "bamread.h"
#include "readvariant.h"

extern int VARIANTS;
extern int PARSEINDELS;
extern int IFLAG;
extern int MAX_IS;

int compare_read_SNP(struct alignedread* read, VARIANT* varlist, int ss, int start, int l1, int l2, FRAGMENT* fragment);

int compare_read_INDEL(struct alignedread* read, VARIANT* varlist, int ss, int start, int l1, int l2, int clength, FRAGMENT* fragment, int CL, REFLIST* reflist);

int compare_read_BND(struct alignedread* read, VARIANT* varlist, int ss, int start, int l1, int l2, FRAGMENT* fragment);

int extract_variants_read(struct alignedread* read, HASHTABLE* ht, CHROMVARS* chromvars, VARIANT* varlist, int paired, FRAGMENT* fragment, int chrom, REFLIST* reflist);

int add_fragment(FRAGMENT* flist, FRAGMENT* fragment, struct alignedread* read, int fragments);
#endif // __PARSEBAMREAD_H__