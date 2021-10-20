
#include "readvariant.h"
#include "readfasta.h"
#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <unordered_map>


int count_variants_hts(char* vcffile, char* sampleid, int* samplecol){
    vcfFile *fp;
    char mode[5] = "r";
    int len = strlen(vcffile);
    if (vcffile[len-1] == 'z' && vcffile[len-2] == 'g')
    {
        mode[1] = 'z'; 
        mode[2] = '\0';
    }


    if ((fp = vcf_open(vcffile, mode)) == 0)
    {
        fprintf(stderr, "\nFail to open VCF file %s\n", vcffile);
        exit(-1);
    }

    bcf_hdr_t *header = bcf_hdr_read(fp);
    bcf1_t *record = bcf_init();
    int variants = 0;
    int ret = 0;
    
    while (true)
    {
        ret = bcf_read1(fp, header, record);
        if (ret < 0)
        break;
        variants++;
    }

    bcf_destroy(record);
    bcf_hdr_destroy(header);
    vcf_close(fp);
    return variants;
}

// count the # of variants in VCF file to allocate space for VCF variant array

int count_variants(char* vcffile, char* sampleid, int* samplecol) {
    FILE* fp = fopen(vcffile, "r");
    if (fp == NULL) {
        fprintf(stderr, "could not open file %s\n\n", vcffile);
        return -1;
    }
    int variants = 0;
    char buffer[500000];
    int i = 0, j = 0, cols = 0;

    while (fgets(buffer, 500000, fp)) {
        if (buffer[0] != '#') variants++; // this should work for non-VCF files as well.
        else if (buffer[0] == '#' && buffer[1] == '#') continue;
        else if (buffer[0] == '#' && buffer[1] == 'C' && buffer[2] == 'H' && buffer[3] == 'R' && buffer[4] == 'O' && buffer[5] == 'M') {
            // find the column of the sample we want to phase
            j = 0;
            while (buffer[i++] != '\n') {
                if ((buffer[i] == ' ' || buffer[i] == '\t') && j == 1) j = 0;
                else if (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n' && j == 0) {
                    j = 1;
                    cols++;
                }
            }
        }
    }
    fclose(fp);
    fprintf(stderr, "VCF file %s has %d variants \n", vcffile, variants);
    return variants;
}

// in VCF format all data lines are tab-delimited but we still allow spaces (why ??)
// we do not check the VCF file for consistency with format, assume that it is in correct format

int parse_variant(VARIANT* variant, char* buffer, int samplecol) {
    int i = 0, j = 0, k = 0, s = 0, e = 0;
    int col = 10;
    int flag = 0;
    char* tempstring;

    // additional variables added so that we can calculate genotype likelihoods and allele counts for each variant using all reads not just haplotype-informative reads
    variant->depth = 0;
    variant->A1 = 0;
    variant->A2 = 0;
    variant->H1 = 0;
    variant->H2 = 0;

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    variant->chrom = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) variant->chrom[j - s] = buffer[j];
    variant->chrom[j - s] = '\0';
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    tempstring = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) tempstring[j - s] = buffer[j];
    tempstring[j - s] = '\0';
    variant->position = atoi(tempstring);
    free(tempstring);

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;

    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i; // varid
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    variant->RA = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) variant->RA[j - s] = buffer[j];
    variant->RA[j - s] = '\0';
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    variant->AA = (char*) malloc(e - s + 1);
    for (j = s; j < e; j++) variant->AA[j - s] = buffer[j];
    variant->AA[j - s] = '\0';

    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;
    while (buffer[i] == ' ' || buffer[i] == '\t') i++;
    s = i;
    while (buffer[i] != ' ' && buffer[i] != '\t') i++;
    e = i;

    // assert that GT field is the first field
    assert(buffer[s] == 'G' && buffer[s+1] == 'T' && (buffer[s+2] == ':' || buffer[s+2] == '\t'));

    while (buffer[i] != '\n') {
        while (buffer[i] == ' ' || buffer[i] == '\t') i++;
        s = i;
        while (buffer[i] != ' ' && buffer[i] != '\t' && buffer[i] != '\n') i++;
        e = i;
        if (col == samplecol) {
            variant->genotype = (char*) malloc(e - s + 1);
            for (j = s; j < e; j++) variant->genotype[j - s] = buffer[j];
            variant->genotype[j - s] = '\0';
        } else col++;
    }

    int gl = strlen(variant->genotype);
    // check that the genotype field is diploid
    if ((gl >= 3 && (variant->genotype[1] == '/' || variant->genotype[1] == '|')) &&
       (gl == 3 || variant->genotype[3] == ':') &&
       (variant->genotype[0] == '0' || variant->genotype[0] == '1' || variant->genotype[0] == '2') &&
       (variant->genotype[2] == '0' || variant->genotype[2] == '1' || variant->genotype[2] == '2')){

        if (variant->genotype[0] != '2' && variant->genotype[2] != '2') // both alleles are 0/1
        {
            variant->allele1 = (char*) malloc(strlen(variant->RA) + 1);
            strcpy(variant->allele1, variant->RA);
            j = 0;
            while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
            variant->allele2 = (char*) malloc(j + 1);
            for (i = 0; i < j; i++) variant->allele2[i] = variant->AA[i];
            variant->allele2[i] = '\0';
            variant->type = strlen(variant->allele2) - strlen(variant->allele1);
        } else if (variant->genotype[0] == '0' || variant->genotype[2] == '0') // at least one allele is reference
        {
            variant->allele1 = (char*) malloc(strlen(variant->RA) + 1);
            strcpy(variant->allele1, variant->RA);
            j = 0;
            while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
            k = j + 1;
            while (variant->AA[k] != ',' && k < strlen(variant->AA)) k++;
            variant->allele2 = (char*) malloc(k - j + 1);
            for (i = j + 1; i < k; i++) variant->allele2[i - j - 1] = variant->AA[i];
            variant->allele2[i - j - 1] = '\0';
            variant->type = strlen(variant->allele2) - strlen(variant->allele1);
        } else // reference allele is missing 1/2 case
        {
            j = 0;
            while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
            variant->allele1 = (char*) malloc(j + 1);
            for (i = 0; i < j; i++) variant->allele1[i] = variant->AA[i];
            variant->allele1[i] = '\0';
            k = j + 1;
            while (variant->AA[k] != ',' && k < strlen(variant->AA)) k++;
            variant->allele2 = (char*) malloc(k - j + 1);
            for (i = j + 1; i < k; i++) variant->allele2[i - j - 1] = variant->AA[i];
            variant->allele2[i - j - 1] = '\0';
            variant->type = strlen(variant->allele2) - strlen(variant->allele1);
            //fprintf(stderr,"non-reference het variant %s %s \n",variant->allele1,variant->allele2);
            // need to allow for multiple alternate alleles or discard ones where there is a comma in alternate allele list
            //j=0; while (variant->AA[j] != ',' && j < strlen(variant->AA)) j++;
            //variant->allele1 = (char*)malloc(j+1); for (i=0;i<j;i++) variant->allele1[i] = variant->AA[i]; variant->allele1[i] = '\0';
            //variant->type = strlen(variant->allele2) -strlen(variant->allele1);
        }

        // reduce the length of the two alleles for VCF format outputted by samtoools, april 17 2012
        // basically CAAAA -> CAA  can be reduced to CA -> C
        i = strlen(variant->allele1) - 1;
        j = strlen(variant->allele2) - 1;
        flag = 0;
        while (i > 0 && j > 0) {
            if (variant->allele1[i] != variant->allele2[j]) break;
            i--;
            j--;
            flag++;
        }
        variant->allele1[i + 1] = '\0';
        variant->allele2[j + 1] = '\0';

        if (variant->type != 0) variant->position++; // add one to position for indels

        if (variant->type != 0)  //if indel, consider 0/1 only
        if ((variant->genotype[0] == '0' && variant->genotype[2] == '1') || (variant->genotype[0] == '1' && variant->genotype[2] == '0')) {
            //if (flag >0) fprintf(stderr,"%s %d %s %s \n",variant->chrom,variant->position,variant->allele1,variant->allele2);
            variant->heterozygous = '1'; // variant will be used for outputting hairs
            //fprintf(stdout,"variant %s %s %s %c\n",variant->allele1,variant->allele2,variant->genotype,variant->heterozygous);
            return 1;
        }
        if ((variant->genotype[0] == '0' && variant->genotype[2] == '2') || (variant->genotype[0] == '2' && variant->genotype[2] == '0')) {
            //if (flag >0) fprintf(stderr,"%s %d %s %s \n",variant->chrom,variant->position,variant->allele1,variant->allele2);
            variant->heterozygous = '1'; // variant will be used for outputting hairs
            //fprintf(stdout,"variant %s %s %s %c\n",variant->allele1,variant->allele2,variant->genotype,variant->heterozygous);
            return 1;
        } else if (variant->genotype[0] != variant->genotype[2] && variant->genotype[0] != '.' && variant->genotype[2] != '.') {
            if (TRI_ALLELIC == 1) fprintf(stderr, "non-ref het variant %d %s %s %s %s %s\n", variant->position, variant->RA, variant->AA, variant->allele1, variant->allele2, variant->genotype);
            // if both alleles are different from reference, we ignore it
            variant->heterozygous = '2';
            return 0;
        } else {
            variant->heterozygous = '0';
            return 0;
        }
    } else {
        fprintf(stdout, "\nERROR: Non-diploid VCF entry detected. Each VCF entry must have a diploid genotype (GT) field consisting of two alleles in the set {0,1,2} separated by either \'/\' or \'|\'. For example, \"1/1\", \"0/1\", and \"0|2\" are valid diploid genotypes for HapCUT2, but \"1\", \"0/3\", and \"0/0/1\" are not.\nThe invalid entry is: \n\n%s\n", buffer);
        exit(1);
    }
        //free(variant->genotype); free(variant->AA); free(variant->RA); free(variant->chrom);
}

// change this to VCF file now

//int read_variantfile(char* vcffile, VARIANT* varlist, HASHTABLE* ht, int* hetvariants, int samplecol) {
//    FILE* fp = fopen(vcffile, "r");
//    char buffer[500000];
//    int i = 0;
//    //	char allele1[256]; char allele2[256]; char genotype[256]; int quality;
//    char prevchrom[256];
//    strcpy(prevchrom, "----");
//    int chromosomes = 0; //int blocks=0;
//    *hetvariants = 0;
//    int het = 0;
//
//    while (fgets(buffer, 500000, fp)) {
//        if (buffer[0] == '#') continue;
//        else {
//            het = parse_variant(&varlist[i], buffer, samplecol);
//            (*hetvariants) += het;
//            //if (het ==0) continue; else (*hetvariants)++;
//            //	fprintf(stdout,"%s %d %s %s %s %s\n",varlist[i].chrom,varlist[i].position,varlist[i].RA,varlist[i].AA,varlist[i].genotype,prevchrom);
//            if (strcmp(varlist[i].chrom, prevchrom) != 0) {
//                //	fprintf(stderr,"chromosomes %d %d\n",chromosomes,i);
//                // insert chromname into hashtable
//                insert_keyvalue(ht, varlist[i].chrom, strlen(varlist[i].chrom), chromosomes);
//                strcpy(prevchrom, varlist[i].chrom);
//                chromosomes++;
//            }
//            i++;
//        }
//    }
//    fclose(fp); //chromosomes--;
//    fprintf(stderr, "vcffile %s chromosomes %d hetvariants %d %d\n", vcffile, chromosomes, *hetvariants, i);
//    return chromosomes;
//
//}

//TODO: clarify between 1/2 genotyped breakend and breakend with multiple mate, they both have two fields in alt alleles. might cause a bug
// now we don't consider multiple mate yet
int parse_bnd(VARIANT *variant)
{
    if (variant->bnd == 0)
    {
        fprintf(stderr, "bug here at position %i", variant->position);
        exit(1);
    }
    if (!STDBND) {
        variant->bnd_mate_chrom = variant->chrom;
        variant->bnd_type = BNDTYPE_PAIRED;
    } else {
        char *ori_bnd_str = variant->allele2;
        char *bnd_str = (char*) malloc(strlen(variant->allele2) + 1);
        strcpy(bnd_str, variant->allele2);
        char *pch;
        pch = strtok(bnd_str, "][:");

        if (strchr(ori_bnd_str, '[') != NULL)
        {
            if (ori_bnd_str[0] == '[')
            {
                variant->bnd_direction = BNDDIRECT_PRT; // [p[t
                variant->bnd_mate_chrom = pch;
                variant->bnd_mate_pos = atoi(strtok(NULL, "][:"));
            }
            else {
                variant->bnd_direction = BNDDIRECT_TRP; // t[p[
                variant->bnd_mate_chrom = strtok(NULL, "][:");
                variant->bnd_mate_pos = atoi(strtok(NULL, "][:"));
            }
        }
        else if (strchr(ori_bnd_str, ']') != NULL)
        {
            if (ori_bnd_str[0] == ']')
            {
                variant->bnd_direction = BNDDIRECT_PLT; // ]p]t
                variant->bnd_mate_chrom = pch;
                variant->bnd_mate_pos = atoi(strtok(NULL, "][:"));
            }
            else {
                variant->bnd_direction = BNDDIRECT_TLP; // t]p]
                variant->bnd_mate_chrom = strtok(NULL, "][:");
                variant->bnd_mate_pos = atoi(strtok(NULL, "][:"));
            }
        }

        else {
            //single breakend
            variant->bnd_type = BNDTYPE_SINGLE_END;
            variant->bnd_pair_distance = -1;
        }
        free(bnd_str);
    }

    if (variant->bnd_type != BNDTYPE_SINGLE_END)  
    {
        variant->bnd_type = BNDTYPE_PAIRED;
        variant->bnd_pair_distance = variant->bnd_mate_pos - variant->position;
        if (variant->bnd_pair_distance < 0) 
            variant->bnd_pair_distance *= -1;

        if (strcmp(variant->bnd_mate_chrom, variant->chrom) != 0)
        {    
            variant->bnd_type = BNDTYPE_INTRA_CHROMOSOME;
            variant->bnd_pair_distance = -1;
        }
    }
    fprintf(stderr, "BND\t%d\t%d\n", variant->position, variant->bnd_mate_pos);
    return 0;
}

int  parse_variant_hts(VARIANT *variant, bcf1_t *record, const bcf_hdr_t *header)
{
    variant->depth = 0;
    variant->A1 = 0;
    variant->A2 = 0;
    variant->H1 = 0;
    variant->H2 = 0;
    variant->phase_set = 0;

    const char *chrom = bcf_seqname(header, record);
    variant->chrom = (char*) malloc(strlen(chrom) + 3+1);
    strcpy(variant->chrom, "chr");
    strcat(variant->chrom, chrom);
    variant->id = (char*) malloc(strlen(record->d.id) + 1);
    strcpy(variant->id, record->d.id);
    variant->position = record->pos + 1;
    //fprintf(stderr, "%d\n", record->pos);
    variant->altalleles =  record->n_allele - 1;
    
    variant->RA = (char *) malloc(strlen(record->d.allele[0]) + 1);
    strcpy(variant->RA, record->d.allele[0]);

    variant->genotype = (char*) malloc(4);
    int ngt_arr = 0, ngt = 0, *gt = NULL;
    int nps_arr = 0, nps = 0, *ps = NULL;
    if (VCF_PHASED) {
        nps = bcf_get_format_int32(header, record, "PS", &ps, &nps_arr);
        if (nps == 1) variant->phase_set = *ps;
    }

    ngt = bcf_get_format_int32(header, record, "GT", &gt, &ngt_arr);
    if (ngt > 2)
    {
        fprintf(stdout, "\nERROR: Non-diploid VCF entry detected. Each VCF entry must have a diploid genotype (GT) field consisting of two alleles in the set {0,1,2} separated by either \'/\' or \'|\'. For example, \"1/1\", \"0/1\", and \"0|2\" are valid diploid genotypes for HapCUT2, but \"1\", \"0/3\", and \"0/0/1\" are not.\n");
        exit(1);
    }
    char t = bcf_gt_allele(gt[0]);
//    fixme for gvcf genotype ./. t = -1
    if (t == -1) return 0;
    variant->genotype[0] = bcf_gt_allele(gt[0]) + '0';
    variant->genotype[1] = '/';
    variant->genotype[2] = bcf_gt_allele(gt[1]) + '0';
    variant->genotype[3] = '\0';
    
    free(gt); 
    int gt_len = strlen(variant->genotype);

    int len = 0;
    int i = 1;
    for (; i <= variant->altalleles; i++)
        len += strlen(record->d.allele[i]) + 1; 
    variant->AA = (char *) malloc(len);
    
    int idx = 0;
    for (i = 1; i <= variant->altalleles; i++)
    {
        len = strlen(record->d.allele[i]);
        int j = 0;
        for (j = 0; j < len; j++)
            variant->AA[idx++] = record->d.allele[i][j];
        if (i < variant->altalleles)
            variant->AA[idx++] = ',';
        else variant->AA[idx++] = '\0';
    }

    int ninfo_arr = 0, ninfo = 0;
    char * info = NULL;

    ninfo = bcf_get_info_string(header, record, "SVTYPE", &info, &ninfo_arr);

    if (ninfo < 0)
        variant->bnd = 0;
    else {
        variant->bnd = 1;
//        if (strcmp(info, "BND") == 0)
//        {
//            variant->bnd = 1;
//        }
//        else    variant->bnd = 0;
    }

    free(info);
    


    int j, k;
    if (((variant->genotype[1] == '/' || variant->genotype[1] == '|')) && gt_len == 3){
        //for mhc phasing, you should not care , just find the two allele 
        if (variant->genotype[0] == '0' || variant->genotype[2] == '0')     // 0/x case
        {
            variant->allele1 = (char*) malloc(strlen(variant->RA) + 1);
            strcpy(variant->allele1, variant->RA);
            int index = 0;
            if (variant->genotype[0] == '0')
                index = variant->genotype[2] - '0';
            else 
                index = variant->genotype[0] - '0';
            
            variant->allele2 = (char*) malloc(strlen(record->d.allele[index]) + 1);
            strcpy(variant->allele2, record->d.allele[index]);
            variant->type = strlen(variant->allele2) - strlen(variant->allele1);
        }

        else {                                                              // x/y case
            int index1 = 0, index2 = 0;
            if (variant->genotype[0] > variant->genotype[2])
            {
                index1 = variant->genotype[2] - '0';
                index2 = variant->genotype[0] - '0';
            }
            else 
            {
                index1 = variant->genotype[0] - '0';
                index2 = variant->genotype[2] - '0';
            }

            variant->allele1 = (char *) malloc(strlen(record->d.allele[index1]) + 1);
            variant->allele2 = (char *) malloc(strlen(record->d.allele[index2]) + 1);
            strcpy(variant->allele1, record->d.allele[index1]);
            strcpy(variant->allele2, record->d.allele[index2]);
            //TODO: potential bug here, make sure type is only referenced when bnd == 0
            variant->type = strlen(variant->allele2) - strlen(variant->allele1);
        }
        // reduce the length of the two alleles for VCF format outputted by samtoools, april 17 2012
        // basically CAAAA -> CAA  can be reduced to CA -> C
        
        if (variant->bnd == 0)
        {
            i = strlen(variant->allele1) - 1;
            j = strlen(variant->allele2) - 1;
            int flag = 0;
            while (i > 0 && j > 0) {
                if (variant->allele1[i] != variant->allele2[j])
                    break;
                i--;
                j--;
                flag++;
            }
            variant->allele1[i + 1] = '\0';
            variant->allele2[j + 1] = '\0';

            if (variant->type != 0)
                variant->position++; // add one to position for indels
        }
        else {
            int end_info;
            int end_info_arr = 0;
            ninfo = bcf_get_info_string(header, record, "AC", &end_info, &end_info_arr);
            if (ninfo != 0)
                variant->bnd_mate_pos = end_info;
//            free(end_info);
            parse_bnd(variant);
        }
        
        if (variant->type != 0 || variant->bnd == 1)      // if indel or bnd, consider 0/1 only
        {
            if ((variant->genotype[0] == '0' && variant->genotype[2] == '1') || (variant->genotype[0] == '1' && variant->genotype[2] == '0')) {
                //if (flag >0) fprintf(stderr,"%s %d %s %s \n",variant->chrom,variant->position,variant->allele1,variant->allele2);
                variant->heterozygous = '1'; // variant will be used for outputting hairs
                //fprintf(stdout,"variant %s %s %s %c\n",variant->allele1,variant->allele2,variant->genotype,variant->heterozygous);
                return 1;
            }
            else {
                variant->heterozygous = '0';
                return 0;
            }
        }
        
        // snp only here 
        if ((variant->genotype[0] == '0' && variant->genotype[2] == '1') || (variant->genotype[0] == '1' && variant->genotype[2] == '0')) {
            //if (flag >0) fprintf(stderr,"%s %d %s %s \n",variant->chrom,variant->position,variant->allele1,variant->allele2);
            variant->heterozygous = '1'; // variant will be used for outputting hairs
            //fprintf(stdout,"variant %s %s %s %c\n",variant->allele1,variant->allele2,variant->genotype,variant->heterozygous);
            return 1;
        }
        if ((variant->genotype[0] == '0' && variant->genotype[2] == '2') || (variant->genotype[0] == '2' && variant->genotype[2] == '0')) {
            //if (flag >0) fprintf(stderr,"%s %d %s %s \n",variant->chrom,variant->position,variant->allele1,variant->allele2);
            variant->heterozygous = '1'; // variant will be used for outputting hairs
            //fprintf(stdout,"variant %s %s %s %c\n",variant->allele1,variant->allele2,variant->genotype,variant->heterozygous);
            return 1;
        } else if (variant->genotype[0] != variant->genotype[2] && variant->genotype[0] != '.' && variant->genotype[2] != '.') {
            if (TRI_ALLELIC == 1) fprintf(stderr, "non-ref het variant %d %s %s %s %s %s\n", variant->position, variant->RA, variant->AA, variant->allele1, variant->allele2, variant->genotype);
            // if both alleles are different from reference, we ignore it
            variant->heterozygous = '2';
            return 0;
        } else {
            variant->heterozygous = '0';
            return 0;
        }
        /*if ((variant->genotype[0] == '0' && variant->genotype[2] >= '2') || (variant->genotype[0] >= '2' && variant->genotype[2] == '0')) {
            //if (flag >0) fprintf(stderr,"%s %d %s %s \n",variant->chrom,variant->position,variant->allele1,variant->allele2);
            variant->heterozygous = '1'; // variant will be used for outputting hairs
            //fprintf(stdout,"variant %s %s %s %c\n",variant->allele1,variant->allele2,variant->genotype,variant->heterozygous);
            return 1;
        } else if (variant->genotype[0] != variant->genotype[2] && variant->genotype[0] != '.' && variant->genotype[2] != '.') {
            if (TRI_ALLELIC == 1) fprintf(stderr, "non-ref het variant %d %s %s %s %s %s\n", variant->position, variant->RA, variant->AA, variant->allele1, variant->allele2, variant->genotype);
            // if both alleles are different from reference, we ignore it
            variant->heterozygous = '2';
            return 0;
        } else {
            variant->heterozygous = '0';
            return 0;
        }
    } else {
        fprintf(stdout, "\nERROR: Non-diploid VCF entry detected. Each VCF entry must have a diploid genotype (GT) field consisting of two alleles in the set {0,1,2} separated by either \'/\' or \'|\'. For example, \"1/1\", \"0/1\", and \"0|2\" are valid diploid genotypes for HapCUT2, but \"1\", \"0/3\", and \"0/0/1\" are not.\n%d\n", strlen(variant->genotype));
        exit(1);
    }*/
    }
}

int read_variantfile_hts(char *vcffile, VARIANT *varlist, HASHTABLE *ht, int *hetvariants, std::unordered_map<std::string,std::pair<int, int>>& BNDs)
{
    vcfFile *fp; 
    if ((fp = vcf_open(vcffile, "r")) == 0)
    {
        fprintf(stderr, "\nFail to open VCF file %s\n", vcffile);
        exit(-1);
    }

    bcf_hdr_t *header = bcf_hdr_read(fp);
    bcf1_t *record = bcf_init();
    int i = 0;

    char prevchrom[256];
    strcpy(prevchrom, "----");
    int chromosomes = 0;
    *hetvariants = 0;
    int het = 0;

    while (bcf_read1(fp, header, record) >= 0)
    {
        bcf_unpack(record, BCF_UN_ALL);
        het = parse_variant_hts(&varlist[i], record, header);
        if(varlist[i].bnd == 1) {
//            TODO, here delimiter only work for svaba
            char* token = strtok(varlist[i].id, ":");
            if (BNDs.find(token) == BNDs.end()) {
                auto tmp = std::make_pair<int, int>(i+1, 0);
                BNDs[token] = tmp;
            } else {
                BNDs[token].second = i + 1;
            }
        }
        (*hetvariants) += het;

        //if (het ==0) continue; else (*hetvariants)++;
            //	fprintf(stdout,"%s %d %s %s %s %s\n",varlist[i].chrom,varlist[i].position,varlist[i].RA,varlist[i].AA,varlist[i].genotype,prevchrom);
        if (strcmp(varlist[i].chrom, prevchrom) != 0) {
                //	fprintf(stderr,"chromosomes %d %d\n",chromosomes,i);
                // insert chromname into hashtable
            insert_keyvalue(ht, varlist[i].chrom, strlen(varlist[i].chrom), chromosomes);
            strcpy(prevchrom, varlist[i].chrom);
            chromosomes++;
        }

        i++;
    }

    bcf_destroy(record);
    bcf_hdr_destroy(header);
    vcf_close(fp);

    fprintf(stderr, "vcffile %s chromosomes %d hetvariants %d %d\n", vcffile, chromosomes, *hetvariants, i);
    return chromosomes;
}

// build a physical map that maps  intervals on chromosomes to the first variant that precedes the start of that interval

void build_intervalmap(CHROMVARS* chromvars, int chromosomes, VARIANT* varlist, int snps) {
    int i = 0, j = 0, k = 0, blocks = 0;
    chromvars[j].first = 0;
    j = 0;
    for (i = 0; i < snps - 1; i++) {
        if (strcmp(varlist[i].chrom, varlist[i + 1].chrom) != 0) {
            chromvars[j].last = i;
            chromvars[j].variants = chromvars[j].last - chromvars[j].first + 1;
            //fprintf(stderr,"new chrom %d %d %s %s\n",j,chromvars[j].variants,varlist[i].chrom,varlist[i+1].chrom);
            j++;
            chromvars[j].first = i + 1;
        }
    }
    chromvars[j].last = i;
    //	int** intervalmap; // map 1000bp of chromosome to first snp in that region indexed by snp_array
    // first SNP to the right of the given base position including that position
    for (j = 0; j < chromosomes; j++) {
        blocks = (int) (varlist[chromvars[j].last].position / BSIZE) + 2;
        chromvars[j].blocks = blocks;
        //	fprintf(stderr,"chromosomes %d blocks %d \n",j,blocks);
        chromvars[j].intervalmap = (int*) malloc(sizeof (int)*blocks);
        for (i = 0; i < blocks; i++) chromvars[j].intervalmap[i] = -1;
        //fprintf(stderr,"blocks for chrom %d: %d \n",j,blocks);
        k = chromvars[j].first;
        for (i = 0; i < blocks; i++) {
            if (k == chromvars[j].last){
                chromvars[j].intervalmap[i] = k;
                continue;
            }
            while (varlist[k].position < BSIZE * i && k < chromvars[j].last) k++;
            if (varlist[k].position >= BSIZE * i && chromvars[j].intervalmap[i] == -1) chromvars[j].intervalmap[i] = k;
            //if (k == chromvars[j].last) break;

            //		if (chromvars[j].intervalmap[i] != -1) printf("FSNPtoright chrom %d block %d: %d %d \n",j,BSIZE*i,chromvars[j].intervalmap[i],varlist[chromvars[j].intervalmap[i]].position);
            //			else printf("FSNPtoright chrom %d block %d: %d \n",j,BSIZE*i,intervalmap[j][i]);
        }
    }
}

// this will only work for pure insertions and deletions, not for block substitutions
int calculate_rightshift(VARIANT* varlist, int ss, REFLIST* reflist) {
    int i = 0, j = 0;
    int a1 = 0, a2 = 0;
    int shift = 0;
    a1 = strlen(varlist[ss].allele1);
    a2 = strlen(varlist[ss].allele2);
    if (a1 > a2 && a2 == 1)
    {
        i = varlist[ss].position; // first base of deletion assuming position is +1 and not previous base
        j = varlist[ss].position + a1 - a2;
        while (i - 1 < reflist->lengths[reflist->current] && j - 1 < reflist->lengths[reflist->current] && reflist->sequences[reflist->current][i - 1] == reflist->sequences[reflist->current][j - 1]) {
            i++;
            j++;
            shift++;
        }
        varlist[ss].shift = shift;
        return shift;
    } 
    else if (a1 == 1 && a2 > a1) 
    {
        i = 1;
        j = varlist[ss].position;
        while (j - 1 < reflist->lengths[reflist->current] && varlist[ss].allele2[i] == reflist->sequences[reflist->current][j - 1] && i < a2) {
            i++;
            j++;
            shift++;
        }
        if (i == a2) // covered the full length of the inserted bases
        {
            i = varlist[ss].position;
            while (i - 1 < reflist->lengths[reflist->current] && j - 1 < reflist->lengths[reflist->current] && reflist->sequences[reflist->current][i - 1] == reflist->sequences[reflist->current][j - 1]) {
                i++;
                j++;
                shift++;
            }
        }
        varlist[ss].shift = shift;
        return shift;
    } 
    else 
	{
		varlist[ss].shift=0;
		return 0;
	}
}



