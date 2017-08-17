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
    uint32_t bx;            // binary representation of the BX string, e.g. CTACCGTTCCGCCTAT
    uint32_t
        read1:1,            // is this the first read in pair
        flag_exclude:1,     // exclude read because of .. BAM flag filters
        sclip_exclude:1,    // .. number soft clipped bases > max_soft_clips
        mq_exclude:1,       // .. MQ < min_mq
        tid:28,
        pos;
    uint8_t read_len;       // read length after soft clipping
}
elem_t;

typedef struct
{
    uint32_t iread, len;
}
gap_t;

typedef struct
{
    uint64_t 
        nmapped_bases,             // number of mapped bases (soft-clips removed)
        nmin_read_spacing,         // reads too close (< min_read_spacing), possible duplicates
        nmin_fragments,            // barcodes excluded because no good fragment was found (all isolated reads)
        nreads,                    // number of reads total
        nprocessed,                // total reads
        nbarcodes,                 // number of unique barcodes
        nflag_exclude,             // reads excluded because of BAM flag filters
        nmq_exclude,               // excluded because of low MQ (< min_mq)
        nsoft_clips,               // excluded because too many soft clips (> max_soft_clips)
        nmin_barcode_reads;        // too few reads per barcode (< nmin_barcode_reads), possibly an error in the barcode sequence
    cov_t 
        *cov_reads,           // for coverage distribution by reads
        *cov_frags;           // for coverage distribution by fragments
    dist_t
        *d_nfrags,           // distribution of fragment count per barcode
        *d_frag_nreads,      // distribution of read count per fragment
        *d_bx_nreads,        // distribution of read count per barcode
        *d_frag_size,        // fragment size distribution
        *d_sclip1,           // distribution of soft clips in first reads
        *d_sclip2;           // distribution of soft clips in second reads
}
stats_t;

typedef struct
{
    uint64_t
        nflag[16],              // BAM flag stats
        nbarcodeseq_has_n,      // the barcode was not corrected, contains N's
        nno_bx;                 // excluded because no BX tag
    stats_t all, good;
    gap_t *gaps;            // list of read indexes and the distance to the next read
    int mbreaks, *breaks;   // gap indices
    int mgaps, melem;
    elem_t *elem;
    uint16_t exclude_flag, min_mq;
    uint64_t checksum_wr, checksum_ro;
    double frag_len;                // expected fragment size, normally distributed
    int frag_len_prior;
    htsFile *bam_fh;
    bam_hdr_t *bam_hdr;
    FILE **fh;
    char **files, *tmp_dir;
    uint32_t *cnt;
    int argc, nfiles, debug, max_soft_clips, min_barcode_reads;
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
    int i, pass = 1;

    // flag stats
    for (i=0; i<16; i++) if ( rec->core.flag&(1<<i) ) args->nflag[i]++;

    elem_t elem;
    memset(&elem, 0, sizeof(elem));

    // exclude masked reads
    if ( rec->core.flag & args->exclude_flag )
    {
        pass = 0;
        elem.flag_exclude = 1;
        args->all.nflag_exclude++;
    }

    // exclude reads with low mapping quality
    else if ( rec->core.qual < args->min_mq )
    {
        pass = 0;
        elem.mq_exclude = 1;
        args->all.nmq_exclude++;
    }

    // exclude reads without the BX tag
    const char *tag = NULL;
    uint8_t *aux = bam_aux_get(rec,"BX");   
    if ( !aux || !(tag = bam_aux2Z(aux)) ) args->nno_bx++;

    // count soft clips
    int nsclip = 0;
    for (i=0; i<rec->core.n_cigar; i++)
    {
        int cig = bam_cigar_op(bam_get_cigar(rec)[i]);
        if ( cig!=BAM_CSOFT_CLIP ) continue;
        nsclip += bam_cigar_oplen(bam_get_cigar(rec)[i]);
    }
    if ( nsclip > args->max_soft_clips )
    {
        pass = 0;
        elem.sclip_exclude = 1;
        args->all.nsoft_clips++;
    }
    elem.read1 = rec->core.flag & BAM_FREAD1 ? 1 : 0;
    elem.read_len = rec->core.l_qseq - nsclip < 256 ? rec->core.l_qseq - nsclip : 255;

    if ( tag )
    {
        for (i=0; i<16; i++)
        {
            int base = char2nt(tag[i]);
            if ( base>3 ) { args->nbarcodeseq_has_n++; tag = NULL; break; }
            elem.bx |= base << (i*2);
        }
        if ( tag && tag[16] )
        {
            if ( tag[16] != '-' || tag[17] != '1' ) error("TODO: expected \"-1\" suffix, found \"%s\"\n", tag);
        }
    }
    if ( !tag ) return;
    if ( rec->core.flag&BAM_FUNMAP ) return;

    elem.pos = rec->core.pos;
    elem.tid = rec->core.tid;
    dist_insert(elem.read1 ? args->all.d_sclip1 : args->all.d_sclip2, elem.read_len);
    cov_insert(args->all.cov_reads, elem.tid, elem.pos, elem.pos + elem.read_len - 1);
    args->all.nmapped_bases += elem.read_len;
    if ( pass )
    {
        args->good.nmapped_bases += elem.read_len;
        cov_insert(args->good.cov_reads, elem.tid, elem.pos, elem.pos + elem.read_len - 1);
    }

    // sanity check: does the numeric representation of the barcode work?
    char str[17];
    bx_int2str(elem.bx, str);
    if ( strncmp(str,tag,16) ) error("Wrong BX hashing!! %s vs %s\n", tag,str);
    args->checksum_wr += elem.bx;

    i = elem.bx % args->nfiles;
    if ( fwrite(&elem, sizeof(elem_t), 1, args->fh[i]) != 1 ) error("Could not write %d bytes\n", sizeof(elem_t));
    args->cnt[i]++;
}

/*
    beg: points to a list of reads sorted by coordinate
    n:   number of reads
*/
int analyze_fragments(args_t *args, elem_t *elem, int n)
{
    if ( n==1 ) return 0;
    dist_insert(args->all.d_frag_size, elem[n-1].pos - elem[0].pos);
    dist_insert(args->all.d_frag_nreads, n);
    return 1;
}

void analyze_barcode(args_t *args, elem_t *elem, uint32_t cnt)
{
    args->all.nbarcodes++;
    args->all.nreads += cnt;
    dist_insert(args->all.d_bx_nreads, cnt);

    // use only good reads in a barcode
    int nelem = 0, i;
    hts_expand(elem_t, cnt, args->melem, args->elem);
    for (i=0; i<cnt; i++)
    {
        if ( elem[i].flag_exclude ) continue;
        if ( elem[i].sclip_exclude ) continue;
        if ( elem[i].mq_exclude ) continue;
        args->elem[nelem++] = elem[i];
    }

    int good_barcode = 1;
    if ( nelem < args->min_barcode_reads ) { args->good.nmin_barcode_reads++; good_barcode = 0; }

    int iend, ibeg = 0, nfrag = 0;
    while ( ibeg < nelem )
    {
        for (iend=ibeg+1; iend<nelem; iend++)
            if ( args->elem[ibeg].tid != args->elem[iend].tid ) break;
    
        nfrag += analyze_fragments(args, &args->elem[ibeg], iend-ibeg);
    
        ibeg = iend;
    }
    if ( nfrag )
        dist_insert(args->good.d_nfrags, nfrag);
    else
    {
        args->good.nmin_fragments++;
        good_barcode = 0;
    }

    if ( good_barcode )
    {
        args->good.nbarcodes++;
        dist_insert(args->good.d_bx_nreads, cnt);
    }
}

void analyze(args_t *args)
{
    uint32_t i, max = 0;
    for (i=0; i<args->nfiles; i++) if ( max < args->cnt[i] ) max = args->cnt[i];
    elem_t *elem = (elem_t*) malloc(sizeof(elem_t)*max);
    for (i=0; i<args->nfiles; i++)
    {
        if ( fclose(args->fh[i])!=0 ) error("Failed to close the file: %s\n", args->files[i]);

        FILE *fp = fopen(args->files[i],"r");
        if ( !fp ) error("Failed to open temporary file for reading %s: %s\n", args->files[i],strerror(errno));
        int n = fread(elem, sizeof(elem_t), args->cnt[i], fp);
        if ( n!=args->cnt[i] ) error("Read failed, read %d out of %d elements: %s .. %s\n", n,args->cnt[i],args->files[i],strerror(errno));
        if ( fclose(fp)!=0 ) error("Failed to close the file: %s\n", args->files[i]);

        unlink(args->files[i]);
        qsort(elem, args->cnt[i], sizeof(elem_t), cmp_elem);

        uint32_t ibeg = 0, iend;
        while ( ibeg < args->cnt[i] )
        {
            for (iend = ibeg; iend < args->cnt[i]; iend++)
            {
                if ( elem[ibeg].bx != elem[iend].bx ) break;
                args->checksum_ro += elem[iend].bx;
            }
            analyze_barcode(args, elem+ibeg, iend-ibeg);
            ibeg = iend;
        }
    }

    free(elem);
}

void report(args_t *args)
{
    if ( args->checksum_wr != args->checksum_ro )
    {
        fprintf(stderr,"Error: checksum failed .. %"PRIu64" vs %"PRIu64"\n", args->checksum_wr,args->checksum_ro);
        printf("ERR\tchecksum_failed\t%"PRIu64"\t%"PRIu64"\n", args->checksum_wr,args->checksum_ro);
    }

    uint32_t i, n, beg, end;
    printf("CMD\t%s", args->argv[0]);
    for (i=1; i<args->argc; i++) printf(" %s",args->argv[i]);
    printf("\n");

    printf("LM\tflags\t0x%x\n", args->exclude_flag);
    printf("LM\tminMQ\t%d\n", args->min_mq);
    printf("SN\tn_all_reads\t%"PRIu64"\n", args->all.nreads + args->nno_bx);
    printf("SN\tn_bx_reads\t%"PRIu64"\n", args->all.nreads);
    printf("SN\tn_good_reads\t%"PRIu64"\n", args->good.nreads);
    printf("SN\tn_mapped_bx_bases\t%"PRIu64"\n", args->all.nmapped_bases);
    printf("SN\tn_mapped_good_bases\t%"PRIu64"\n", args->good.nmapped_bases);
    printf("SN\tn_all_barcodes\t%"PRIu64"\n", args->all.nbarcodes);
    printf("SN\tn_excluded\t%"PRIu64"\n", args->all.nmq_exclude + args->all.nsoft_clips + args->all.nflag_exclude + args->nbarcodeseq_has_n);
    printf("SN\tn_excluded_bx_has_n\t%"PRIu64"\n", args->nbarcodeseq_has_n);
    printf("SN\tn_excluded_mq\t%"PRIu64"\t%d\n", args->all.nmq_exclude,args->min_mq);
    if ( args->max_soft_clips!=INT32_MAX ) printf("SN\tn_excluded_soft_clips\t%"PRIu64"\t%d\n", args->all.nsoft_clips,args->max_soft_clips);
    printf("# Note: n_excluded_flag can be smaller than the sum of n_flag_* counts because one read can have multiple bits set\n");
    printf("SN\tn_excluded_flag\t%"PRIu64"\n", args->all.nflag_exclude);
    for (i=0; i<16; i++)
    {
        char *str = bam_flag2str(1<<i);
        if ( *str )
        {
            if ( args->exclude_flag & (1<<i) )
            {
                printf("SN\tn_flag_%s\t%"PRIu64"\texcluded\n", str, args->nflag[i]);
            }
            else
                printf("SN\tn_flag_%s\t%"PRIu64"\n", str, args->nflag[i]);
        }
        free(str);
    }
    printf("SN\tn_good_barcodes\t%"PRIu64"\n", args->good.nbarcodes);
    printf("SN\tn_barcodes_excluded\t%"PRIu64"\n", args->all.nbarcodes - args->good.nbarcodes);
    printf("SN\tn_barcodes_excluded_min_reads\t%"PRIu64"\t%d\n", args->good.nmin_barcode_reads, args->min_barcode_reads);
    printf("SN\tn_barcodes_excluded_min_frags\t%"PRIu64"\n", args->good.nmin_fragments);

    printf("# DIST_ALL_SCLIP1, Number of bases left after soft-clipping in first reads (all reads)\n");
    n = dist_n(args->all.d_sclip1);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->all.d_sclip1, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_ALL_SCLIP1\t%u\t%"PRIu64"\n", beg, cnt);
    }

    printf("# DIST_ALL_SCLIP2, Number of bases left after soft-clipping in second reads (all reads)\n");
    n = dist_n(args->all.d_sclip2);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->all.d_sclip2, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_ALL_SCLIP2\t%u\t%"PRIu64"\n", beg, cnt);
    }

    printf("# DIST_COV_ALL_READS, Histogram of genome coverage, all reads\n");
    n = cov_n(args->all.cov_reads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = cov_get(args->all.cov_reads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_COV_ALL_READS\t%u\t%"PRIu64"\n", beg, cnt);
    }

    printf("# DIST_COV_GOOD_READS, Histogram of genome coverage, good reads (filtered by flags, MQ, ...)\n");
    n = cov_n(args->good.cov_reads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = cov_get(args->good.cov_reads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("DIST_COV_GOOD_READS\t%u\t%"PRIu64"\n", beg, cnt);
    }

    printf("# BX_NALL_READS, Number of reads per barcode (all reads)\n");
    n = dist_n(args->all.d_bx_nreads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->all.d_bx_nreads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("BX_NALL_READS\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# BX_NGOOD_READS, Number of good reads per barcode\n");
    n = dist_n(args->good.d_bx_nreads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->good.d_bx_nreads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("BX_NGOOD_READS\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# FRAG_SIZE, Fragment size\n");
    n = dist_n(args->all.d_frag_size);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->all.d_frag_size, i, &beg, &end);
        if ( !cnt ) continue;
        printf("FRAG_SIZE\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }

    printf("# FRAG_NREADS, Number of reads per fragment\n");
    n = dist_n(args->all.d_frag_nreads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = dist_get(args->all.d_frag_nreads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("FRAG_NREADS\t%u\t%u\t%"PRIu64"\n", beg, end, cnt);
    }


#if 0
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
#endif
}

void init_stats(stats_t *stats)
{
    memset(stats, 0, sizeof(*stats));
    stats->d_nfrags         = dist_init(4);
    stats->d_frag_size      = dist_init(4);
    stats->d_frag_nreads    = dist_init(4);
    stats->d_bx_nreads      = dist_init(4);
    stats->d_sclip1         = dist_init(3);        
    stats->d_sclip2         = dist_init(3);
    stats->cov_reads        = cov_init(1000);
    stats->cov_frags        = cov_init(10000);
}
void destroy_stats(stats_t *stats)
{
    dist_destroy(stats->d_nfrags);
    dist_destroy(stats->d_frag_size);
    dist_destroy(stats->d_frag_nreads);
    dist_destroy(stats->d_bx_nreads);
    dist_destroy(stats->d_sclip1);
    dist_destroy(stats->d_sclip2);
    cov_destroy(stats->cov_reads);
    cov_destroy(stats->cov_frags);
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

    init_stats(&args->all);
    init_stats(&args->good);

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
}

void destroy(args_t *args)
{
    destroy_stats(&args->good);
    destroy_stats(&args->all);

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

#if 0
    args->good.cov_reads = cov_init();
    cov_insert(args->good.cov_reads, 2, 1, 1);
    cov_insert(args->good.cov_reads, 2, 5, 6);
    cov_insert(args->good.cov_reads, 2, 6, 7);
    cov_insert(args->good.cov_reads, 2, 10, 19);
    cov_insert(args->good.cov_reads, 2, 10, 19);
    cov_insert(args->good.cov_reads, 2, 10, 19);
    cov_insert(args->good.cov_reads, 2, 10, 19);

    uint32_t beg,end;
    int i,n = cov_n(args->good.cov_reads);
    for (i=0; i<n; i++)
    {
        uint64_t cnt = cov_get(args->good.cov_reads, i, &beg, &end);
        if ( !cnt ) continue;
        printf("COV_HIST_GOOD_READS\t%d\t%"PRIu64"\n", beg, cnt);
    }
    cov_destroy(args->good.cov_reads);
    return 0;
#endif

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
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

    analyze(args);
    report(args);

    bam_destroy1(bam_line);
    bam_hdr_destroy(args->bam_hdr);
    if ( hts_close(args->bam_fh)!=0 ) error("Close failed: %s\n",strcmp("-",bam_fname)?bam_fname:"STDIN");
    destroy(args);

    return 0;
}
