/*
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "version.h"
#include "bxcheck.h"

int main_trim(int argc, char *argv[]);
int main_stats(int argc, char *argv[]);

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help;
}
cmd_t;

static cmd_t cmds[] =
{
    { .func = main_trim,
      .alias = "trim",
      .help = "trim and error-correct barcode sequences, create fastq files suitable for mapping"
    },
    { .func = main_stats,
      .alias = "stats",
      .help = "collect BAM file stats, the output can be plotted with plot-bxcheck.py"
    },
    { .func  = NULL,
      .alias = NULL,
      .help  = NULL
    }
};

char *bxcheck_version(void)
{
    return BXCHECK_VERSION;
}

static void usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Program: bxcheck (Tools for processing 10x Genomics Chromium data)\n");
    fprintf(fp, "License: MIT/Expat license\n");
    fprintf(fp, "Version: %s\n", bxcheck_version());
    fprintf(fp, "\n");
    fprintf(fp, "Usage:   bxcheck [--version|--version-only] [--help] <command> <argument>\n");
    fprintf(fp, "\n");
    fprintf(fp, "Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].alias)
    {
        if ( !cmds[i].func ) sep = cmds[i].alias;
        if ( sep )
        {
            fprintf(fp, "\n -- %s\n", sep);
            sep = NULL;
        }
        if ( cmds[i].func && cmds[i].help[0]!='-' ) fprintf(fp, "    %-12s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }
    fprintf(fp,"\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(stderr); return 1; }

    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0) 
    {
        printf("bxcheck %s\nCopyright (C) 2016 Genome Research Ltd.\n", bxcheck_version());
        printf("License: The MIT/Expat license\n");
        printf("This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n");
        return 0;
    }
    else if (strcmp(argv[1], "--version-only") == 0) 
    {
        printf("%s\n", bxcheck_version());
        return 0;
    }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        if (argc == 2) { usage(stdout); return 0; }
        // Otherwise change "bxcheck help COMMAND [...]" to "bxcheck COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }

    int i = 0;
    while (cmds[i].alias)
    {
        if (cmds[i].func && strcmp(argv[1],cmds[i].alias)==0)
        {
            return cmds[i].func(argc-1,argv+1);
        }
        i++;
    }
    fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
    return 1;
}

