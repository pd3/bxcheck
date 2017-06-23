/*  bxcheck.c -- Collect chromium BX stats

    Copyright (C) 2017 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.  
*/

#include <unistd.h>     // for isatty()
#include <inttypes.h>   // for PRIu64
#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <htslib/sam.h>
#include "dist.h"
#include "cov.h"

#define FRAG_LEN_NORM   1
#define FRAG_LEN_EXP    2
#define FRAG_LEN_FIXED  3

typedef struct 
{
    uint32_t bx;        // binary representation of the BX string, e.g. CTACCGTTCCGCCTAT
    uint32_t suffix:2,  // the BX suffix, e.g. "1" from "CTACCGTTCCGCCTAT-1"
             pass:1,    // does the readpair pass local filters? (MQ, softclips, flags)
             tid:29, pos;
}
elem_t;

typedef struct
{
    uint32_t iread, len;
}
gap_t;

typedef struct
{
    int mbreaks, *breaks;   // indices of the gaps
    int mgaps, melem;
    gap_t *gaps;            // list of read indexes and the distance to the next read
    elem_t *elem;
    dist_t *d_nfrags, *d_frag_read_dist, *d_frag_size, *d_frag_nreads, *d_frag_min_dist, *d_bx_nreads;
    dist_t *d_sclip1, *d_sclip2, *d_sclip1_all, *d_sclip2_all;
    cov_t *cov_reads, *cov_frags;
    uint16_t exclude_flag, min_mq;
    uint64_t nreads, nprocessed, nexcluded;    // the number reads total and the number of unfiltered/filtered reads
    uint64_t nbarcodes, nflag_exclude, nflag[16], ntid, nmin_mq, nno_bx, nsoft_clips, nmin_barcode_reads;
    uint64_t nmin_read_spacing, nmin_reads_per_fragment, nbarcodeseq_has_n;
    uint64_t checksum_wr, checksum_ro;
    double frag_len;                // expected fragment size, normally distributed
    int frag_len_prior;
    htsFile *bam_fh;
    bam_hdr_t *bam_hdr;
    FILE **fh, *resume_fp;
    char **files, *tmp_dir, *resume_fname;
    uint32_t *cnt;
    int argc, nfiles, stop_after_sort, resume_after_sort, debug, max_soft_clips, min_barcode_reads;
    int min_reads_per_fragment, min_read_spacing;
    char **argv;
}
args_t;

static inline int cmp_elem(const void *a, const void *b)
{
    elem_t *x = (elem_t *) a;
    elem_t *y = (elem_t *) b;

    if ( x->bx < y->bx ) return -1;
    if ( x->bx > y->bx ) return 1;

    if ( x->suffix < y->suffix ) return -1;
    if ( x->suffix > y->suffix ) return 1;

    if ( x->tid < y->tid ) return -1;
    if ( x->tid > y->tid ) return 1;

    if ( x->pos < y->pos ) return -1;
    if ( x->pos > y->pos ) return 1;

    return 0;
}
static inline int cmp_gap_len_desc(const void *_a, const void *_b)
{
    gap_t *a = (gap_t *) _a;
    gap_t *b = (gap_t *) _b;
    if ( a->len < b->len ) return 1;
    else if ( a->len > b->len ) return -1;
    return 0;
}
static inline int cmp_int_asc(const void *_a, const void *_b)
{
    int *a = (int *) _a;
    int *b = (int *) _b;
    if ( *a < *b ) return -1;
    else if ( *a > *b ) return 1;
    return 0;
}

volatile sig_atomic_t fatal_error_in_progress = 0;
static args_t *cleanup_data = NULL;

void fatal_error_signal(int sig)
{
    if ( fatal_error_in_progress ) raise(sig);
    fatal_error_in_progress = 1;

    if ( cleanup_data )
    {
        fprintf(stderr,"Signal caught, cleaning up %s...\n", cleanup_data->tmp_dir);
        int i;
        for (i=0; i<cleanup_data->nfiles; i++)
            unlink(cleanup_data->files[i]);
        rmdir(cleanup_data->tmp_dir);
    }

    signal(sig, SIG_DFL);
    raise(sig);
}

static void error(const char *format, ...)
{
    if ( !format )
    {
        printf("About: The program collects BX statistics from BAM files\n");
        printf("Usage: bxcheck [OPTIONS] file.bam\n");
        printf("Options:\n");
        printf("        --debug                     Print debugging information\n");
        printf("        --prior-exp                 Fragment lengths distributed exponentially [this is the default]\n");
        printf("        --prior-fixed               Fragment lengths determined using a fixed threshold\n");
        printf("        --prior-norm                Fragment lengths distributed normally\n");
        printf("        --resume-after-sort DIR     For debugging: resume after the barcode sorting step\n");
        printf("        --stop-after-sort           For debugging: stop after the barcode sorting step\n");
        printf("    -b, --reads-per-barcode INT     Minimum number of good reads per barcode [2]\n");
        printf("    -e, --exclude-flag INT|STR      Reads to exclude [UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY]\n");
        printf("    -f, --fragment-size NUMBER      Expected fragment length, assuming normal distribution [5e4]\n");
        printf("    -m, --mapping-qual INT          Minimum mapping quality for the fragment size distribution [20]\n");
        printf("    -n, --n-temp-files INT          Create up to INT temporary files while sorting by barcode [100]\n");
        printf("    -r, --reads-per-fragment INT    Minimum number of good reads per fragment [2]\n");
        printf("    -R, --read-spacing INT          Minimum read spacing [10]\n");
        printf("    -s, --max-soft-clips INT        Exclude reads with more than INT soft-clipped bases [inf]\n");
        printf("    -t, --temp-dir STR              Temporary directory for sorting reads [/tmp/bxcheck.XXXXXX]\n");
        printf("\n");
    }
    else
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    exit(-1);
}

inline uint8_t char2nt(char nt)
{
    if ( nt=='A' || nt=='a' ) return 0;
    if ( nt=='C' || nt=='c' ) return 1;
    if ( nt=='G' || nt=='g' ) return 2;
    if ( nt=='T' || nt=='t' ) return 3;
    return 4;
}

void bx_int2str(uint32_t num, char *bx)
{
    int i,j;
    for (i=0; i<16; i++)
    {
        j = (num >> (i*2)) & 0x3;
        bx[i] = "ACGT"[j];
    }
    bx[16] = 0;
}

void process(args_t *args, const bam1_t *rec)
{
    uint32_t i;
    args->nreads++;

    // flag stats
    for (i=0; i<16; i++) if ( rec->core.flag&(1<<i) ) args->nflag[i]++;

    // exclude unmapped reads
    if ( rec->core.flag & args->exclude_flag ) { args->nflag_exclude++; args->nexcluded++; return; }
    if ( rec->core.tid < 0 ) { args->ntid++; args->nexcluded++; return; }

    // exclude reads with low mapping quality
    if ( rec->core.qual < args->min_mq ) { args->nmin_mq++; args->nexcluded++; return; }

    // exclude reads without the BX tag
    uint8_t *aux = bam_aux_get(rec,"BX");   
    if ( !aux ) { args->nno_bx++; args->nexcluded++; return; }
    const char *tag = bam_aux2Z(aux);
    if ( !tag ) { args->nno_bx++; args->nexcluded++; return; }

    // count soft clips
    int nsclip = 0;
    for (i=0; i<rec->core.n_cigar; i++)
    {
        int cig = bam_cigar_op(bam_get_cigar(rec)[i]);
        if ( cig!=BAM_CSOFT_CLIP ) continue;
        nsclip += bam_cigar_oplen(bam_get_cigar(rec)[i]);
    }
    if ( rec->core.flag & BAM_FREAD1 ) 
        dist_insert(args->d_sclip1_all, rec->core.l_qseq - nsclip);
    else if ( rec->core.flag & BAM_FREAD2 )
        dist_insert(args->d_sclip2_all, rec->core.l_qseq - nsclip);

    if ( nsclip > args->max_soft_clips ) { args->nsoft_clips++; args->nexcluded++; return; }

    // Add only first read from a pair and use their midpoint
    elem_t elem;
    elem.bx = 0;
    for (i=0; i<16; i++)
    {
        int base = char2nt(tag[i]);
        if ( base>3 ) { args->nbarcodeseq_has_n++; args->nexcluded++; return; }
        elem.bx |= base << (i*2);
    }
    if ( tag[16] )
    {
        if ( tag[17] != '1' ) error("TODO: expected \"-1\" suffix, found \"%s\"\n", tag);
        elem.suffix = tag[17] - '1';
    }
    else elem.suffix = 0;
    elem.tid = rec->core.tid;

    args->nprocessed++;

    cov_insert(args->cov_reads, rec->core.tid, rec->core.pos, rec->core.pos);

    if ( rec->core.flag & BAM_FREAD1 ) 
        dist_insert(args->d_sclip1, rec->core.l_qseq - nsclip);
    else if ( rec->core.flag & BAM_FREAD2 )
        dist_insert(args->d_sclip2, rec->core.l_qseq - nsclip);


    int use_midpoint = 1;
    if ( !(rec->core.flag & BAM_FPROPER_PAIR) ) use_midpoint = 0;
    if ( rec->core.flag & BAM_FMUNMAP ) use_midpoint = 0;
    if ( rec->core.tid != rec->core.mtid ) use_midpoint = 0;
    if ( abs(rec->core.pos - rec->core.mpos) > args->frag_len ) use_midpoint = 0;
    if ( use_midpoint )
    {
        if ( rec->core.flag & BAM_FREAD2 ) return;
        elem.pos = (rec->core.pos + rec->core.mpos)/2;
    }
    else
        elem.pos = rec->core.pos;

    // sanity check
    {
        char str[17];
        bx_int2str(elem.bx, str);
        if ( strncmp(str,tag,16) ) error("Wrong BX hashing!! %s vs %s\n", tag,str);
    }
    args->checksum_wr += elem.bx;

    i = elem.bx % args->nfiles;
    if ( fwrite(&elem, sizeof(elem_t), 1, args->fh[i]) != 1 ) error("Could not write %d bytes\n", sizeof(elem_t));
    args->cnt[i]++;
}

double calc_len_prob(args_t *args, int len)
{
    // normal or exponential distribution??
    // plot [0:1e5]  exp(-(x-5e4)**2/(0.25*5e4)**2), exp(-x/5e4)
    if ( args->frag_len_prior == FRAG_LEN_NORM )
    {
        // normal distribution
        double sigma2 = (0.25 * args->frag_len)*(0.25 * args->frag_len);
        return -0.5 * (len - args->frag_len)*(len - args->frag_len) / sigma2;
    }
    else if ( args->frag_len_prior == FRAG_LEN_EXP )
    {
        // exponential distribution
        return -(double)len/args->frag_len;
    }
    else if ( args->frag_len_prior == FRAG_LEN_FIXED )
    {
        // fixed threshold, twice the mean length
        return len >= 2*args->frag_len ? -HUGE_VAL : 0;
    }
    else error("Uh: frag_len_prior=%d\n",args->frag_len_prior);
    return 0;
}
double calc_bird_prob(args_t *ags, int pnt, int tot)
{
    if ( !tot ) return 0;

    double frac = (double)pnt / tot;

    // parameters fitted, see #1492521104 in ~/wtxt/logs/sandbox/10x-genomics/ChangeLog
    double mean  = 0.389507 + (303981.8 - 0.389507)/(1 + pow(frac/0.00001360044,1.085307));
    double sigma = 0.0002933294 + (79967.68 - 0.0002933294)/(1 + pow(frac/0.000208926,1.351487));
    if ( frac >= mean ) return 0;

    // Be more lenient, make the distribution wider to account for the fact that we 
    // do not know it exactly
    sigma *= 2;

    return -0.5*(frac-mean)*(frac-mean) / (sigma*sigma);
}

double calc_prob(args_t *args, elem_t *reads, int nreads, int *breaks, int nbreaks)
{
    // sort the breaks by the read index (which is monotonous with the read coordinate)
    if ( nbreaks ) qsort(breaks, nbreaks, sizeof(*breaks), cmp_int_asc);

    // add a top bound to form the last fragment
    if ( !nbreaks || breaks[nbreaks-1] != nreads - 1 ) breaks[nbreaks++] = nreads - 1;

    double prob = 0;
    int ibeg = 0, i,j;
    for (i=0; i<nbreaks; i++)
    {
        int pnt = 0, tot = 0, iend = breaks[i];
        if ( iend - ibeg == 1 )
        {
            pnt = tot = reads[iend].pos - reads[ibeg].pos;
            prob += calc_len_prob(args, tot) + calc_bird_prob(args, pnt, tot);
        }
        else
        {
            int must_paint = 1;
            int len_left = reads[ibeg+1].pos - reads[ibeg].pos;
            for (j=ibeg+1; j<iend; j++)
            {
                int len_right = reads[j+1].pos - reads[j].pos;

                // To paint a shorter segment only once, we paint only segments on the left
                // Shorter was the segment on the right in the previous step, paint in now
                if ( must_paint ) pnt += len_left;

                // The right segment is shorter, paint it next time
                if ( len_left > len_right ) must_paint = 1;

                // The left segment is shorter and it was not painted above
                else if ( !must_paint ) pnt += len_left; 

                // The left segment is shorter and it was painted
                else must_paint = 0;

                tot += len_left;
                len_left = len_right;
            }
            pnt += len_left;
            tot += len_left;
            prob += calc_len_prob(args, tot) + calc_bird_prob(args, pnt, tot);
        }
        ibeg = iend + 1;
    }
    return prob;
}

int analyze_fragment(args_t *args, elem_t *beg, int n)
{
    if ( n < args->min_reads_per_fragment )
    {
        args->nmin_reads_per_fragment++;
        return 0;
    }

    int i, dist = INT32_MAX;
    for (i=1; i<n; i++)
    {
        if ( dist > beg[i].pos - beg[i-1].pos ) dist = beg[i].pos - beg[i-1].pos;
    }
    dist_insert(args->d_frag_read_dist, dist);
    dist_insert(args->d_frag_nreads, n);
    dist_insert(args->d_frag_size, beg[n-1].pos - beg[0].pos);
    cov_insert(args->cov_frags, beg[0].tid, beg[0].pos,beg[n-1].pos);
    return 1;
}

/*
    beg: points to a list of reads sorted by coordinate
    n:   number of reads
*/
int analyze_fragments(args_t *args, elem_t *beg, int n)
{
    int nfrag = 0;
    if ( n==1 )
    {
        nfrag = analyze_fragment(args, beg, n);  // do not count singletons
        return nfrag;
    }

    int ngaps = n - 1;
    if ( ngaps > args->mgaps )
    {
        args->mgaps = ngaps;
        args->gaps  = (gap_t*) realloc(args->gaps, sizeof(*args->gaps)*args->mgaps);
    }

    int i,j;
    for (i=0; i<ngaps; i++)
    {
        args->gaps[i].iread = i;
        args->gaps[i].len   = beg[i+1].pos - beg[i].pos;
    }
    qsort(args->gaps, ngaps, sizeof(*args->gaps), cmp_gap_len_desc);

    if ( args->mbreaks < 10 )
    {
        args->mbreaks = 10;
        args->breaks = (int*) realloc(args->breaks, sizeof(*args->breaks)*args->mbreaks);
    }
    int nbreaks = 0;
    double prob_prev = calc_prob(args, beg, ngaps+1, args->breaks, nbreaks);

    if ( args->debug )
    {
        char str[17];
        bx_int2str(beg->bx, str);
        printf("%s\tn=%d\n", str,n);
        printf("\t%f\n", prob_prev);
    }

    for (i=0; i<ngaps; i++)
    {
        nbreaks++;
        hts_expand(int, (nbreaks+1), args->mbreaks, args->breaks);
        args->breaks[nbreaks-1] = args->gaps[i].iread;
        double prob = calc_prob(args, beg, ngaps+1, args->breaks, nbreaks);
        if ( args->debug )
        {
            int j;
            printf("\t%f", prob);
            for (j=0; j<nbreaks; j++) printf("\t%d:%d", args->gaps[j].iread,args->gaps[j].len);
            printf("\n");
        }
        if ( prob < prob_prev ) 
        {
            nbreaks--;
            break;
        }
        prob_prev = prob;
    }

    if ( nbreaks )
    {
        for (i=0; i<nbreaks; i++) args->breaks[i] = args->gaps[i].iread;
        qsort(args->breaks, nbreaks, sizeof(*args->breaks), cmp_int_asc);
    }

    i = 0;
    while (i <= nbreaks)
    {
        int ibeg = i==0 ? 0 : args->breaks[i-1] + 1;
        int iend = i==nbreaks ? ngaps : args->breaks[i];
        nfrag += analyze_fragment(args, beg+ibeg, iend-ibeg+1);
        i++;

        if ( args->debug )
        {
            printf("\t[dist=%.1fMb len=%.1fkb, %d, %d-%d]",ibeg==0?0:(beg[ibeg].pos-beg[ibeg-1].pos)*1e-6,(beg[iend].pos - beg[ibeg].pos)*1e-3,beg[ibeg].pos,ibeg,iend);
            for (j=ibeg; j<=iend; j++) 
            {
                if ( j==ibeg ) printf(" 0");
                else printf(" %d", beg[j].pos - beg[j-1].pos);
            }
            printf("\n");
        }
    }
    if ( args->debug ) printf("\n\n");

    if ( nbreaks>0 )
    {
        double min = args->gaps[0].len;
        for (i=1; i<nbreaks; i++)
            if ( min > args->gaps[i].len ) min = args->gaps[i].len;
        dist_insert(args->d_frag_min_dist, min);
    }

    return nfrag;
}

void analyze_barcode(args_t *args, elem_t *elem, uint32_t cnt)
{
    args->nbarcodes++;

    if ( cnt < args->min_barcode_reads )
    {
        args->nmin_barcode_reads++;
        return;
    }

    hts_expand(elem_t, cnt, args->melem, args->elem);
    int i, j = 1;
    args->elem[0] = elem[0];
    for (i=1; i<cnt; i++)
    {
        if ( elem[i].tid==args->elem[j-1].tid && elem[i].pos - args->elem[j-1].pos < args->min_read_spacing )
        {
        assert( i!=j-1 );
            args->nmin_read_spacing++;
            continue;
        }
        args->elem[j++] = elem[i];
    }
    cnt = j;

    if ( cnt < args->min_barcode_reads )
    {
        args->nmin_barcode_reads++;
        return;
    }

    int iend, ibeg = 0, nfrag = 0;
    while ( ibeg < cnt )
    {
        for (iend=ibeg+1; iend<cnt; iend++)
            if ( args->elem[ibeg].tid != args->elem[iend].tid ) break;

        nfrag += analyze_fragments(args, &args->elem[ibeg], iend-ibeg);

        ibeg = iend;
    }
    dist_insert(args->d_nfrags, nfrag);
    if ( nfrag > 0 )
        dist_insert(args->d_bx_nreads, cnt);
}

void stop_after_sort(args_t *args)
{
    uint32_t i, max = 0;
    for (i=0; i<args->nfiles; i++) if ( max < args->cnt[i] ) max = args->cnt[i];

    elem_t *elem = (elem_t*) malloc(sizeof(elem_t)*max);

    for (i=0; i<args->nfiles; i++)
    {
        if ( fclose(args->fh[i])!=0 ) error("Failed to close the file: %s\n", args->files[i]);
        fprintf(args->resume_fp,"%s\t%d\n", args->files[i], args->cnt[i]);

        FILE *fp = fopen(args->files[i],"r");
        if ( !fp ) error("Failed to open temporary file for reading %s: %s\n", args->files[i],strerror(errno));
        int n = fread(elem, sizeof(elem_t), args->cnt[i], fp);
        if ( n!=args->cnt[i] ) error("Read failed: read %d out of %d elements: %s\n", n,args->cnt[i],args->files[i]);
        if ( fclose(fp)!=0 ) error("Failed to close the file: %s\n", args->files[i]);

        qsort(elem, args->cnt[i], sizeof(elem_t), cmp_elem);

        char *tmp = (char*) malloc(strlen(args->files[i])+6);
        sprintf(tmp,"%s.part",args->files[i]);
        fp = fopen(tmp,"w");
        if ( !fp ) error("Failed to open temporary file for writing %s: %s\n", tmp,strerror(errno));
        n = fwrite(elem, sizeof(elem_t), args->cnt[i], fp);
        if ( n!=args->cnt[i] ) error("Write failed: wrote %d out of %d elements: %s\n", n,args->cnt[i],tmp);
        if ( fclose(fp)!=0 ) error("Failed to close the file: %s\n", tmp);

        if ( rename(tmp,args->files[i])!=0 ) error("Rename failed: %s -> %s: %s\n", tmp,args->files[i],strerror(errno));

        unlink(tmp);
        free(tmp);
    }
    if ( fclose(args->resume_fp)!=0 ) error("Failed to close the file: %s\n", args->resume_fname);
    fprintf(stderr,"To resume, run the program as: bxcheck --resume-after-sort %s\n", args->resume_fname);

    free(elem);
}

void analyze(args_t *args)
{
    uint32_t i, max = 0;
    for (i=0; i<args->nfiles; i++) if ( max < args->cnt[i] ) max = args->cnt[i];

    elem_t *elem = (elem_t*) malloc(sizeof(elem_t)*max);

    for (i=0; i<args->nfiles; i++)
    {
        if ( !args->resume_after_sort && fclose(args->fh[i])!=0 ) error("Failed to close the file: %s\n", args->files[i]);

        FILE *fp = fopen(args->files[i],"r");
        if ( !fp ) error("Failed to open temporary file for reading %s: %s\n", args->files[i],strerror(errno));
        int n = fread(elem, sizeof(elem_t), args->cnt[i], fp);
        if ( n!=args->cnt[i] ) error("Read failed, read %d out of %d elements: %s .. %s\n", n,args->cnt[i],args->files[i],strerror(errno));
        if ( fclose(fp)!=0 ) error("Failed to close the file: %s\n", args->files[i]);

        if ( !args->resume_after_sort )
        {
            unlink(args->files[i]);
            qsort(elem, args->cnt[i], sizeof(elem_t), cmp_elem);
        }

        uint32_t ibeg = 0, iend;
        while ( ibeg < args->cnt[i] )
        {
            for (iend = ibeg; iend < args->cnt[i]; iend++)
            {
                if ( elem[ibeg].bx != elem[iend].bx ) break;
                if ( elem[ibeg].suffix != elem[iend].suffix ) break;
                args->checksum_ro += elem[iend].bx;
            }
            analyze_barcode(args, elem+ibeg, iend-ibeg);
            ibeg = iend;
        }
    }

    free(elem);

    if ( args->checksum_wr != args->checksum_ro )
    {
        fprintf(stderr,"Error: checksum failed .. %"PRIu64" vs %"PRIu64"\n", args->checksum_wr,args->checksum_ro);
        printf("ERR\tchecksum_failed\t%"PRIu64"\t%"PRIu64"\n", args->checksum_wr,args->checksum_ro);
    }

    printf("CMD\t%s", args->argv[0]);
    for (i=1; i<args->argc; i++) printf(" %s",args->argv[i]);
    printf("\n");

    printf("SN\tn_reads\t%"PRIu64"\n", args->nreads);
    for (i=0; i<16; i++)
    {
        char *str = bam_flag2str(1<<i);
        if ( *str ) printf("SN\tn_flag_%s\t%"PRIu64"\n", str, args->nflag[i]);
        free(str);
    }
    printf("SN\tn_processed\t%"PRIu64"\n", args->nprocessed);
    printf("SN\tn_excluded\t%"PRIu64"\n", args->nexcluded);
    printf("# Note: n_flag_excluded can be smaller than the sum of n_flag_* counts because reads can have multiple bits set\n");
    printf("SN\tn_flag_excluded\t%"PRIu64"\n", args->nflag_exclude);
    printf("SN\tn_negative_tid\t%"PRIu64"\n", args->ntid);
    printf("SN\tn_min_mq\t%"PRIu64"\t%d\n", args->nmin_mq,args->min_mq);
    printf("SN\tn_max_soft_clips\t%"PRIu64"\t%d\n", args->nsoft_clips,args->max_soft_clips);
    printf("SN\tn_no_BX\t%"PRIu64"\n", args->nno_bx);
    printf("SN\tn_barcodes\t%"PRIu64"\n", args->nbarcodes);
    printf("SN\tn_min_reads_per_barcode\t%"PRIu64"\t%d\n", args->nmin_barcode_reads,args->min_barcode_reads);
    printf("SN\tn_min_read_spacing\t%"PRIu64"\t%d\n", args->nmin_read_spacing,args->min_read_spacing);
    printf("SN\tn_min_reads_per_fragment\t%"PRIu64"\t%d\n", args->nmin_reads_per_fragment,args->min_reads_per_fragment);
    printf("SN\tn_barcode_seq_has_n\t%"PRIu64"\n", args->nbarcodeseq_has_n);

    uint32_t beg, end;

    printf("# DIST_NFRAGS, Number of fragments per barcode\n");
    int n = dist_n(args->d_nfrags);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_nfrags, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_NFRAGS\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_FRAG_NREADS, Number of reads per fragment\n");
    n = dist_n(args->d_frag_nreads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_frag_nreads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_FRAG_NREADS\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_FRAG_SIZE, Distribution of fragment size\n");
    n = dist_n(args->d_frag_size);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_frag_size, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_FRAG_SIZE\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_FRAG_SPACING, Minimum distance between two fragments within a barcode (for debugging)\n");
    n = dist_n(args->d_frag_min_dist);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_frag_min_dist, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_FRAG_SPACING\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_BX_NREADS, Number of read pairs per barcode\n");
    n = dist_n(args->d_bx_nreads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_bx_nreads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_BX_NREADS\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_READ_SPACING, minimum spacing between reads in a fragment\n");
    n = dist_n(args->d_frag_read_dist);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_frag_read_dist, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_READ_SPACING\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_SCLIP1, Number of bases left after soft-clipping in first reads (in processed reads)\n");
    n = dist_n(args->d_sclip1);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_sclip1, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_SCLIP1\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_SCLIP2, Number of bases left after soft-clipping in second reads (in processed reads)\n");
    n = dist_n(args->d_sclip2);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_sclip2, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_SCLIP2\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_ALL_SCLIP1, Number of bases left after soft-clipping in first reads (all reads)\n");
    n = dist_n(args->d_sclip1_all);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_sclip1_all, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_ALL_SCLIP1\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# DIST_ALL_SCLIP2, Number of bases left after soft-clipping in second reads (all reads)\n");
    n = dist_n(args->d_sclip2_all);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->d_sclip2_all, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_ALL_SCLIP2\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# COV_FRAGS, Genome coverage (filtered fragments)\n");
    int rid, nrid = cov_nrid(args->cov_frags);
    for (rid=0; rid<nrid; rid++)
    {
        n = cov_nbin(args->cov_frags, rid);
        if ( !n ) continue;
        for (i=0; i<n; i++)
        {
            uint16_t cnt = cov_get(args->cov_frags, rid, i, &beg, &end);
            if ( !cnt ) continue;
            printf("COV_FRAGS\t%s\t%u\t%u\t%u\n", args->bam_hdr->target_name[rid], beg, end, cnt);
        }
    }

    printf("# COV_READS, Genome coverage (processed reads)\n");
    nrid = cov_nrid(args->cov_reads);
    for (rid=0; rid<nrid; rid++)
    {
        n = cov_nbin(args->cov_reads, rid);
        if ( !n ) continue;
        for (i=0; i<n; i++)
        {
            uint16_t cnt = cov_get(args->cov_reads, rid, i, &beg, &end);
            if ( !cnt ) continue;
            printf("COV_READS\t%s\t%u\t%u\t%u\n", args->bam_hdr->target_name[rid], beg, end, cnt);
        }
    }
}

void init(args_t *args)
{
    args->files = (char**) calloc(args->nfiles, sizeof(char*));
    args->fh    = (FILE**) calloc(args->nfiles, sizeof(FILE*));
    args->cnt   = (uint32_t*) calloc(args->nfiles, sizeof(uint32_t));

    int i,len = strlen(args->tmp_dir);
    if ( len<6 || strcmp("XXXXXX",args->tmp_dir+len-6) )
    {
        char *tmp = (char*) malloc(len+7);
        sprintf(tmp,"%sXXXXXX",args->tmp_dir);
        args->tmp_dir = tmp;
    }
    else
        args->tmp_dir = strdup(args->tmp_dir);
    char *tmp_dir = mkdtemp(args->tmp_dir);
    if ( !tmp_dir ) error("mkdtemp(%s) failed: %s\n", args->tmp_dir,strerror(errno));
    fprintf(stderr,"Writing to tempdir %s\n", tmp_dir);

    if ( args->stop_after_sort )
    {
        args->resume_fname = (char*) malloc(strlen(tmp_dir)+13);
        sprintf(args->resume_fname,"%s/resume.txt",tmp_dir);
        args->resume_fp = fopen(args->resume_fname,"w");
    }

    cleanup_data = args;
    signal(SIGINT, fatal_error_signal);

    len = strlen(tmp_dir);
    for (i=0; i<args->nfiles; i++)
    {
        args->files[i] = (char*) malloc(len+11);
        snprintf(args->files[i],len+11,"%s/%d",tmp_dir,i);
        args->fh[i] = fopen(args->files[i],"w");
        if ( !args->fh[i] ) error("Failed to open temporary file %s: %s\n", args->files[i],strerror(errno));
    }

    args->d_nfrags         = dist_init(4);
    args->d_frag_read_dist = dist_init(4);
    args->d_frag_size      = dist_init(4);
    args->d_frag_nreads    = dist_init(4);
    args->d_frag_min_dist  = dist_init(4);
    args->d_bx_nreads      = dist_init(4);
    args->d_sclip1         = dist_init(3);        
    args->d_sclip2         = dist_init(3);
    args->d_sclip1_all     = dist_init(3);        
    args->d_sclip2_all     = dist_init(3);
    args->cov_reads        = cov_init(1000);
    args->cov_frags        = cov_init(10000);
}
void init_resume_after_sort(args_t *args)
{
    args->files = hts_readlist(args->tmp_dir, 1, &args->nfiles);
    if ( !args->files ) error("Failed to read the file %s: %s\n", args->tmp_dir,strerror(errno));

    args->fh  = (FILE**) calloc(args->nfiles, sizeof(FILE*));
    args->cnt = (uint32_t*) calloc(args->nfiles, sizeof(uint32_t));

    int i;
    for (i=0; i<args->nfiles; i++)
    {
        int len = strlen(args->files[i]);
        char *tmp = args->files[i] + len - 1;
        while ( tmp > args->files[i] && isdigit(*tmp) ) tmp--;
        args->cnt[i] = strtol(tmp+1,NULL,10);
        *tmp = 0;
    }

    args->d_nfrags         = dist_init(4);
    args->d_frag_read_dist = dist_init(4);
    args->d_frag_size      = dist_init(4);
    args->d_frag_nreads    = dist_init(4);
    args->d_frag_min_dist  = dist_init(4);
    args->d_bx_nreads      = dist_init(4);
    args->d_sclip1         = dist_init(3);        
    args->d_sclip2         = dist_init(3);
    args->d_sclip1_all     = dist_init(3);        
    args->d_sclip2_all     = dist_init(3);
}
void destroy_resume_after_sort(args_t *args)
{
    dist_destroy(args->d_nfrags);
    dist_destroy(args->d_frag_read_dist);
    dist_destroy(args->d_frag_size);
    dist_destroy(args->d_frag_nreads);
    dist_destroy(args->d_frag_min_dist);
    dist_destroy(args->d_bx_nreads);
    dist_destroy(args->d_sclip1);
    dist_destroy(args->d_sclip2);
    dist_destroy(args->d_sclip1_all);
    dist_destroy(args->d_sclip2_all);
    free(args->breaks);
    free(args->gaps);
    free(args->elem);
    int i;
    for (i=0; i<args->nfiles; i++) free(args->files[i]);
    free(args->cnt);
    free(args->fh);
    free(args->files);
    free(args->resume_fname);
    free(args);
}

void destroy(args_t *args)
{
    dist_destroy(args->d_nfrags);
    dist_destroy(args->d_frag_read_dist);
    dist_destroy(args->d_frag_size);
    dist_destroy(args->d_frag_nreads);
    dist_destroy(args->d_frag_min_dist);
    dist_destroy(args->d_bx_nreads);
    dist_destroy(args->d_sclip1);
    dist_destroy(args->d_sclip2);
    dist_destroy(args->d_sclip1_all);
    dist_destroy(args->d_sclip2_all);
    cov_destroy(args->cov_reads);
    cov_destroy(args->cov_frags);
    rmdir(args->tmp_dir);
    free(args->tmp_dir);
    free(args->breaks);
    free(args->gaps);
    free(args->elem);
    int i;
    for (i=0; i<args->nfiles; i++) free(args->files[i]);
    free(args->cnt);
    free(args->fh);
    free(args->files);
    free(args->resume_fname);
    free(args);
}

int main(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->min_mq = 20;
    args->exclude_flag = BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FSUPPLEMENTARY;
    args->tmp_dir = "/tmp/bxcheck.XXXXXX";
    args->nfiles = 100;
    args->frag_len = 5e4;
    args->frag_len_prior = FRAG_LEN_EXP;
    args->max_soft_clips = INT32_MAX;
    args->min_barcode_reads = 2;
    args->min_reads_per_fragment = 2;
    args->min_read_spacing = 10;
    args->argc = argc;
    args->argv = argv;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"stop-after-sort", no_argument, NULL, 1},
        {"resume-after-sort", required_argument, NULL, 2},
        {"debug", no_argument, NULL, 3},
        {"prior-norm", no_argument, NULL, 4},
        {"prior-exp", no_argument, NULL, 5},
        {"prior-fixed", no_argument, NULL, 6},
        {"exclude-flag", required_argument, NULL, 'e'},
        {"reads-per-barcode", required_argument, NULL, 'b'},
        {"reads-per-fragment", required_argument, NULL, 'r'},
        {"read-spacing", required_argument, NULL, 'R'},
        {"max-soft-clips", required_argument, NULL, 's'},
        {"fragment-size", required_argument, NULL, 'f'},
        {"mapping-qual", required_argument, NULL, 'm'},
        {"n-temp-files", required_argument, NULL, 'n'},
        {"temp-dir", required_argument, NULL, 't'},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    int opt;
    while ( (opt=getopt_long(argc,argv,"?hm:n:t:e:s:b:r:R:",loptions,NULL))>0 )
    {
        switch (opt)
        {
            case  1 :
                args->stop_after_sort = 1;
                break;
            case  2 :
                args->resume_after_sort = 1;
                args->tmp_dir = optarg;
                break;
            case  3 :
                args->debug = 1;
                break;
            case  4 :
                args->frag_len_prior = FRAG_LEN_NORM;
                break;
            case  5 :
                args->frag_len_prior = FRAG_LEN_EXP;
                break;
            case  6 :
                args->frag_len_prior = FRAG_LEN_FIXED;
                break;
            case 'e': 
                args->exclude_flag = bam_str2flag(optarg);
                if ( args->exclude_flag<0 ) { fprintf(stderr,"Could not parse --exclude-flag %s\n", optarg); return 1; }
                break;
            case 't': 
                args->tmp_dir = optarg;
                break;
            case 's': 
                if ( strcasecmp("inf",optarg) )
                {
                    args->max_soft_clips = strtod(optarg, &tmp); 
                    if ( tmp==optarg || *tmp ) error("Could not parse the argument: --max-soft-clips %s\n", optarg);
                }
                break;
            case 'R': 
                args->min_reads_per_fragment = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --reads-per-fragment %s\n", optarg);
                break;
            case 'r': 
                args->min_read_spacing = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --read-spacing %s\n", optarg);
                break;
            case 'b': 
                args->min_barcode_reads = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --reads-per-barcode %s\n", optarg);
                break;
            case 'f': 
                args->frag_len = strtod(optarg, &tmp); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --fragment-size %s\n", optarg);
                break;
            case 'm': 
                args->min_mq = strtol(optarg, &tmp, 10); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --mappinq-qual %s\n", optarg);
                break;
            case 'n': 
                args->nfiles = strtol(optarg, &tmp, 10); 
                if ( tmp==optarg || *tmp ) error("Could not parse the argument: --nfiles %s\n", optarg);
                break;
            case '?':
            case 'h': error(NULL);
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( args->resume_after_sort )
    {
        init_resume_after_sort(args);
        analyze(args);
        destroy_resume_after_sort(args);
        return 0;
    }

    char *bam_fname = NULL;
    if ( optind<argc )
        bam_fname = argv[optind++];
    if ( !bam_fname )
    {
        if ( isatty(STDIN_FILENO) ) error(NULL);
        bam_fname = "-";
    }
    if ((args->bam_fh = sam_open(bam_fname, "r")) == 0) error("Failed to open: %s\n", bam_fname);
    args->bam_hdr = sam_hdr_read(args->bam_fh);
    bam1_t *bam_line = bam_init1();

    init(args);

    while (sam_read1(args->bam_fh, args->bam_hdr, bam_line) >= 0) process(args, bam_line);

    if ( args->stop_after_sort )
        stop_after_sort(args);
    else
        analyze(args);

    bam_destroy1(bam_line);
    bam_hdr_destroy(args->bam_hdr);
    if ( hts_close(args->bam_fh)!=0 ) error("Close failed: %s\n",strcmp("-",bam_fname)?bam_fname:"STDIN");
    destroy(args);

    return 0;
}
