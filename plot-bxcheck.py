#!/usr/bin/env python
#
#   Author: petr.danecek@sanger
#   About:  Script for plotting the output of bxcheck
#   Usage:  
#       bxcheck trim -l known-barcodes.txt file.fq.gz -o trimmed
#       bxcheck stats -l known-barcodes.txt mapped.bam > stats.txt
#       plot-bxcheck.py trimmed.txt stats.txt -d plots
#

import sys

def usage():
    print 'Usage: plot-bxcheck.py [OPTIONS] bx-output.txt'
    print 'Options:'
    print '   -d, --dir STR     output directory'
    print '   -n, --name STR    main title [bxcheck stats]'
    print 'Example:'
    print '   bxcheck trim -l known-barcodes.txt file.fq.gz -o trimmed'
    print '   bxcheck stats -l known-barcodes.txt mapped.bam > stats.txt'
    print '   plot-bxcheck.py trimmed.txt stats.txt -d plots'
    print ''
    sys.exit(1)


fnames = []
dir    = None
itype  = 'png'
main_title = None
replot = True

if len(sys.argv) < 2: usage()
args = sys.argv[1:]
while len(args):
    if args[0]=='-d' or args[0]=='--dir': 
        args = args[1:]
        dir  = args[0]
    elif args[0]=='-n' or args[0]=='--name': 
        args  = args[1:]
        main_title = args[0]
    elif args[0]=='--dont-replot':  # for debugging, do not call matplotlib if data file older than existing image
        replot = False
    else:
        fnames.append(args[0])
    args = args[1:]

if dir==None: usage()
if len(fnames)==0: usage()
if main_title==None: main_title = fnames[0]

color = [
    '#E24A33', # ggplot red
    '#777777', # ggplot gray
    '#348ABD', # ggplot blue
    '#FBC15E', # ggplot orange
    '#8EBA42', # ggplot green
]

graphs = []

dists = \
{
    'frag_nreads':
    {
        'key': 'FRAG_NREADS',
        'ylabel': 'Number of pairs',
        'xlabel': 'Read pairs per fragment',
        'hist':   1,
        'xlim':   0.99,
        'ysci':   1,
        'plot_density': 1,
        'handler': 'frag_nreads_frag2readhist',
        'alpha': 0.8,
    },
    'frag_cov':
    {
        'key': 'FRAG_COV', 
        'ylabel': 'Number of fragments',
        'xlabel': 'Fragment coverage',
        'xlim':   0.95,
        'xsci':   1,
        'ysci':   1,
        'alpha': 0.8,
    },
    'frag_size_bases':
    {
        'key': 'FRAG_SIZE',
        'ylabel': 'Number of bases',
        'xlabel': 'Fragment length',
        'plot_density': 1,
        'ysci':   1,
        'xsci':   1,
        'xlog':   1,
        #'cfrac':  1,
        'handler': 'frag_size_to_bases',
        'alpha': 0.8,
        'xlim': [100],
    },
    'frag_size_seq_bases':
    {
        'key': 'FRAG_SIZE_SEQ',
        'ylabel': 'Sequenced bases',
        'xlabel': 'Fragment length',
        'plot_density': 1,
        'ysci':   1,
        'xsci':   1,
        'xlog':   1,
        #'cfrac':  1,
        'handler': 'smooth_dist',
        'alpha': 0.8,
        'xlim': [100],
    },
    'frag_size':
    {
        'key': 'FRAG_SIZE',
        'ylabel': 'Number of fragments',
        'xlabel': 'Fragment length',
        'ylim':   [0.1,1],
        'xlim':   [100],
        'xsci':   1,
        'xlog':   1,
        'ylog':   1,
        'plot_density': 1,
        'handler': 'smooth_dist',
        #'cfrac':  1,
        'alpha': 0.8,
        'line_plot': 1,
    },
    'insert_size':
    {
        'key': 'ISIZE',
        'ylabel': 'Number of read pairs',
        'xlabel': 'Insert size',
        'xlim':   0.98,
        'hist':   1,
        'ysci':   1,
        'plot_density': 1,
        'alpha': 0.8,
    },
}

dists2 = \
{
    'dist_sclip':
    {
        'key1': 'DIST_ALL_SCLIP1',
        'key2': 'DIST_ALL_SCLIP2',
        'label1': 'First reads',
        'label2': 'Second reads',
        'xlabel': 'Number of bases after soft-clipping',
        'ylabel': 'Number of reads',
        'ylog':   1,
    },
    'cov_reads':
    {
        'key1': 'DIST_COV_ALL_READS',
        'key2': 'DIST_COV_GOOD_READS',
        'ylabel': 'Number of sites',
        'xlabel': 'Coverage',
        'label1': 'All mapped reads',
        'label2': 'Good reads',
        'hist':   1,
        'xlim':   0.99,
        'ysci':   1,
        'alpha': 0.8,
    },
    'bx_reads':
    {
        'key1': 'BX_NALL_READS',
        'key2': 'BX_NGOOD_READS',
        'ylabel': 'Number of pairs',
        'xlabel': 'Read pairs per barcode',
        'label1': 'All mapped reads',
        'label2': 'Good reads',
        'ysci':   1,
        'xlim':   0.93,
        'ylim2': 1,
        'hist':  1,
        'alpha': 0.8,
        'plot_density': 1,
        'handler': 'bx_reads_bx2readhist',
    },
}


import matplotlib as mpl
from matplotlib.image import BboxImage
from matplotlib.transforms import Bbox, TransformedBbox
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv, os, numpy, base64, re
csv.register_dialect('tab', delimiter='\t', quoting=csv.QUOTE_NONE)
plt.style.use('ggplot')

mod_time = 0
for fname in fnames:
    mtime = os.path.getctime(fname)
    if mod_time < mtime: mod_time = mtime

try:
    os.makedirs(dir)
except:
    pass

def get_xlim(xdat,ydat,frac):
    xlim = 0
    sum = 0
    tot = 0
    for cnt in ydat: tot += cnt
    for i in range(len(xdat)):
        sum += ydat[i]
        xlim = xdat[i]
        if sum/tot > frac: break
    return xlim

def crop_frag_size(dist):
    global dat
    xlim = get_xlim(dat['ISIZE']['xdat'],dat['ISIZE']['ydat'],dat['ISIZE']['xlim'])
    i = 0
    for i in range(len(dist['xdat'])):
        if dist['xdat'][i] > xlim: break
    dist['xdat'] = dist['xdat'][i:]
    dist['ydat'] = dist['ydat'][i:]
    dat['ISIZE']['xlim_bp'] = xlim

def frag_nreads_frag2readhist(dist):
    global dat
    for i in range(len(dist['xdat'])):
        dist['ydat'][i] = dist['ydat'][i] * dist['xdat'][i]

def frag_size_to_bases(dist):
    xdat = []
    ydat = []
    norm = 1
    for i in range(len(dist['xdat'])):
        x = dist['xdat'][i]
        y = dist['ydat'][i]
        y = y * x
        if len(xdat)==0 or x - xdat[-1] >= 10:
            if len(ydat)>0: ydat[-1] /= norm
            norm = 1
            xdat.append(x)
            ydat.append(y)
        else:
            ydat[-1] += y
            norm += 1
    dist['xdat'] = xdat
    dist['ydat'] = ydat

def smooth_dist(dist):
    xdat = []
    ydat = []
    norm = 1
    for i in range(len(dist['xdat'])):
        x = dist['xdat'][i]
        y = dist['ydat'][i]
        if len(xdat)==0 or x - xdat[-1] >= 10:
            if len(ydat)>0: ydat[-1] /= norm
            norm = 1
            xdat.append(x)
            ydat.append(y)
        else:
            ydat[-1] += y
            norm += 1
    dist['xdat'] = xdat
    dist['ydat'] = ydat

def plot_dist(dist):
    if 'handler' in dist:
        globals()[dist['handler']](dist)
    
    xdat = dist['xdat']
    ydat = dist['ydat']
    
    wh = (7,3)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('#aaaaaa')
    ax1.spines['left'].set_color('#aaaaaa')

    alpha = 0.8
    if 'alpha' in dist: alpha = dist['alpha']

    if 'cfrac' in dist:
        xlim = 1
        if 'xlim' in dist: xlim = dist['xlim']
        frac = 0
        tot  = sum(ydat)
        i = 0
        for i in range(len(xdat)):
            frac += ydat[i]/tot
            ydat[i] = frac
            if frac >= xlim: break
        xdat = xdat[0:i]
        ydat = ydat[0:i]
        if 'ylim' in dist:
            i = 0
            for i in range(len(xdat)):
                if ydat[i] >= dist['ylim'][0]: break
            xdat = xdat[i:]
            ydat = ydat[i:]
        ax1.plot(xdat,ydat,'o-',color=color[0],mec=color[0],alpha=alpha)
        ax1.set_ylim(top=1.05)
    elif 'line_plot' in dist:
        ax1.plot(xdat,ydat,'o-',color=color[0],mec=color[0],alpha=alpha)
    else:
        width = xdat[1] - xdat[0]
        if len(xdat) > 100:
            #ax1.fill_between([x-width*0.5 for x in xdat],[0 for y in ydat],ydat,color=color[0],alpha=alpha)
            ax1.fill_between(xdat,[0 for y in ydat],ydat,color=color[0],alpha=alpha)
        else:
            ax1.bar([x-width*0.5 for x in xdat],ydat,width=width,color=color[0],edgecolor=color[0],alpha=alpha)
        if 'xlim' in dist:
            xlim = get_xlim(xdat,ydat,dist['xlim'])
            ax1.set_xlim(-width,xlim)
        else:
            ax1.set_xlim(-width)

    if 'xlim' in dist and type(dist['xlim']).__name__ == 'list':
        if len(dist['xlim'])==1:
            ax1.set_xlim(dist['xlim'][0])
        else:
            ax1.set_xlim(dist['xlim'])
    
    if 'xsci' in dist and dist['xsci']:
        ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='x')
    
    if 'ysci' in dist and dist['ysci']:
        ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    
    if 'ylog' in dist:
        ax1.set_yscale('symlog')

    if 'xlog' in dist:
        ax1.set_xscale('symlog')

    if 'ylabel' in dist: ax1.set_ylabel(dist['ylabel'])
    if 'xlabel' in dist: ax1.set_xlabel(dist['xlabel'])

    plt.subplots_adjust(bottom=0.18)
    if not replot:
        img = dir+'/'+dist['img']+'.png'
        if os.path.exists(img) and os.path.getctime(img) > mod_time: return
    plt.savefig(dir+'/'+dist['img']+'.png')
    plt.savefig(dir+'/'+dist['img']+'.pdf')
    plt.close()

def dat_to_cfrac(xdat,ydat,xlim):
    frac = 0
    tot  = sum(ydat)
    i = 0
    for i in range(len(xdat)):
        if xdat[i] > xlim: break
        frac += ydat[i]/tot
        ydat[i] = frac
    xdat = xdat[0:i]
    ydat = ydat[0:i]
    return (xdat,ydat)

def bx_reads_bx2readhist(dist):
    global dat
    for i in range(len(dist['dat1']['xdat'])):
        if i==0: continue
        dist['dat1']['ydat'][i] = dist['dat1']['ydat'][i] * dist['dat1']['xdat'][i]
    for i in range(len(dist['dat2']['xdat'])):
        if i==0: continue
        dist['dat2']['ydat'][i] = dist['dat2']['ydat'][i] * dist['dat2']['xdat'][i]

def plot_dist2(dist):
    if 'handler' in dist:
        globals()[dist['handler']](dist)

    wh = (7,3)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('#aaaaaa')
    ax1.spines['left'].set_color('#aaaaaa')

    xdat1 = dist['dat1']['xdat']
    ydat1 = dist['dat1']['ydat']
    xdat2 = dist['dat2']['xdat']
    ydat2 = dist['dat2']['ydat']

    line1 = 'o-'
    line2 = 'o-'
    
    col2 = color[0]
    col1 = color[1]
    
    alpha = 0.8
    if 'alpha' in dist: alpha = dist['alpha']

    plots  = None
    labels = [dist['label1'],dist['label2']]
    loc = 'best'

    if 'hist' in dist and dist['hist']==1:
        width = xdat1[1] - xdat1[0]
        #   try:
        #       width = xdat1[1] - xdat1[0]
        #   except:
        #       print dist
        #       sys.exit(1)
        if len(xdat2)>1 and width > xdat2[1] - xdat2[0]: width = xdat2[1] - xdat2[0]
        if len(xdat1) > 100 or len(xdat2) > 100:
            ax1.fill_between([x-width*0.5 for x in xdat1],[0 for y in ydat1],ydat1,color=col1,alpha=alpha)
            ax1.fill_between([x-width*0.5 for x in xdat2],[0 for y in ydat2],ydat2,color=col2,alpha=alpha)
            plots = [plt.Rectangle((0,0),1,1,color=col1,alpha=alpha),plt.Rectangle((0,0),1,1,color=col2,alpha=alpha)]
        else:
            plots  = ax1.bar([x-width*0.5 for x in xdat1],ydat1,width=width,color=col1,edgecolor=col1,alpha=alpha)
            plots += ax1.bar([x-width*0.5 for x in xdat2],ydat2,width=width,color=col2,edgecolor=col2,alpha=alpha)
        if 'xlim' in dist:
            xlim1 = get_xlim(xdat1,ydat1,dist['xlim'])
            xlim2 = get_xlim(xdat2,ydat2,dist['xlim'])
            if xlim1 < xlim2: xlim1 = xlim2
            ax1.set_xlim(-width,xlim1)
        else:
            ax1.set_xlim(-width)
    elif 'cfrac' in dist and dist['cfrac']==1:
        loc = 'lower right'
        xlim = 1
        if 'xlim' in dist: 
            xlim  = get_xlim(xdat1,ydat1,dist['xlim'])
            xlim2 = get_xlim(xdat2,ydat2,dist['xlim'])
            if xlim < xlim2: xlim = xlim2
        (xdat1,ydat1) = dat_to_cfrac(xdat1,ydat1,xlim)
        (xdat2,ydat2) = dat_to_cfrac(xdat2,ydat2,xlim)
        plots  = ax1.plot(xdat1,ydat1,line1,color=col1,mec=col1,alpha=alpha)
        plots += ax1.plot(xdat2,ydat2,line2,color=col2,mec=col2,alpha=alpha)
    else:
        plots  = ax1.plot(dist['dat1']['xdat'],dist['dat1']['ydat'],line1,color=col1,mec=col1,alpha=alpha)
        plots += ax1.plot(dist['dat2']['xdat'],dist['dat2']['ydat'],line2,color=col2,mec=col2,alpha=alpha)
    
    plt.legend(plots,labels,numpoints=1,markerscale=1,loc=loc,prop={'size':10},frameon=False)
    if 'ylog' in dist and dist['ylog']==1: ax1.set_yscale('symlog')
    if 'ylim2' in dist and dist['ylim2']==1: ax1.set_ylim(0,max(ydat2))
    ax1.set_ylabel(dist['ylabel'])
    ax1.set_xlabel(dist['xlabel'])
    if 'title' in dist:
        ax1.set_title(dist['title'])
    if 'ysci' in dist and dist['ysci']:
        ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    plt.subplots_adjust(bottom=0.17)
    if not replot:
        img = dir+'/'+dist['img']+'.png'
        if os.path.exists(img) and os.path.getctime(img) > mod_time: return
    plt.savefig(dir+'/'+dist['img']+'.png')
    plt.savefig(dir+'/'+dist['img']+'.pdf')
    plt.close()

def plot_err_correct(dist):
    wh = (7,3.5)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('#aaaaaa')
    ax1.spines['left'].set_color('#aaaaaa')
    ax2 = ax1.twinx()

    line1 = 'o-'
    line2 = '--'
    col1 = color[1]
    col2 = color[2]
    col3 = color[0]
    col4 = color[3]
    col5 = color[1]
    alpha = 0.8
    plots  = None
    labels = ['Before, known barcodes','Before, unknown','After, known','After, unknown','Modified']
    loc = 'best'

    xdat1 = [float(x[1]) for x in dist['RAW_KNOWN']['dat']]
    ydat1 = [float(x[3]) for x in dist['RAW_KNOWN']['dat']]
    xdat2 = [float(x[1]) for x in dist['RAW_UNKNOWN']['dat']]
    ydat2 = [float(x[3]) for x in dist['RAW_UNKNOWN']['dat']]
    xdat3 = [float(x[1]) for x in dist['OUT_KNOWN']['dat']]
    ydat3 = [float(x[3]) for x in dist['OUT_KNOWN']['dat']]
    xdat4 = [float(x[1]) for x in dist['OUT_UNKNOWN']['dat']]
    ydat4 = [float(x[3]) for x in dist['OUT_UNKNOWN']['dat']]
    xdat5 = [float(x[1]) for x in dist['MOD_NNEW']['dat']]
    ydat5 = [float(x[3]) for x in dist['MOD_NNEW']['dat']]

    dist['xlim'] = 0.99;
    xlim = get_xlim(xdat1,ydat1,dist['xlim'])
    tmp  = get_xlim(xdat2,ydat2,dist['xlim'])
    if xlim < tmp: xlim = tmp
    ax1.set_xlim(1,xlim)
    ax2.set_xlim(1,xlim)
    ax1.set_xscale('symlog')
    ax2.set_xscale('symlog')
    ax2.grid(False)
    ax2.set_ylabel('Modified reads', color=col5)

    plots  = ax1.plot(xdat1,ydat1,line1,color=col1,mec=col1,alpha=alpha)
    plots += ax1.plot(xdat2,ydat2,line1,color=col2,mec=col2,alpha=alpha)
    plots += ax1.plot(xdat3,ydat3,line1,color=col3,mec=col3,alpha=alpha)
    plots += ax1.plot(xdat4,ydat4,line1,color=col4,mec=col4,alpha=alpha)
    plots += ax2.plot(xdat5,ydat5,line2,color=col5,mec=col5,alpha=alpha,lw=2)
    
    plt.legend(plots,labels,numpoints=1,markerscale=1,loc=loc,prop={'size':10},frameon=False)
    ax1.set_ylabel('Number of reads')
    ax1.set_xlabel('Read pairs per barcode')
    ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    ax2.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    plt.subplots_adjust(bottom=0.18)

    dist['img'] = 'err-correct'
    if not replot:
        img = dir+'/'+dist['img']+'.png'
        if os.path.exists(img) and os.path.getctime(img) > mod_time: return
    plt.savefig(dir+'/'+dist['img']+'.png')
    plt.savefig(dir+'/'+dist['img']+'.pdf')
    plt.close()

def parse_file(fname):
    if fname==None: f = sys.stdin
    else: f = open(fname, 'r')
    dat = {}
    dat['SN'] = {}
    dat['LM'] = {}
    dat['title'] = re.sub(r'^.*/','', re.sub(r'\.\s*$', '', fname))
    reader = csv.reader(f, 'tab')
    for row in reader:
        key = row[0]
        if key[0]=='#': continue
        done = False

        for img in dists:
            if dists[img]['key'] != key: continue
            if img not in dat: 
                dat[img] = {}
                dat[img].update(dists[img])
                dat[img]['img'] = img
            if 'xdat' not in dat[img]:
                dat[img]['xdat'] = []
                dat[img]['ydat'] = []
            dat[img]['xdat'].append(float(row[1]))
            if 'plot_density' in dists[img]:
                dat[img]['ydat'].append(float(row[3])/(float(row[2])-float(row[1])))
            else:
                dat[img]['ydat'].append(float(row[3]))
            done = True
        if done: continue

        for img in dists2:
            dataset = None
            if dists2[img]['key1'] == key: dataset = 'dat1'
            if dists2[img]['key2'] == key: dataset = 'dat2'
            if dataset==None: continue
            if img not in dat:
                dat[img] = {}
                dat[img].update(dists2[img])
                dat[img]['img'] = img
            if dataset not in dat[img]: dat[img][dataset] = {}
            if 'xdat' not in dat[img][dataset]:
                dat[img][dataset]['xdat'] = []
                dat[img][dataset]['ydat'] = []
            dat[img][dataset]['xdat'].append(float(row[1]))
            if 'plot_density' in dists2[img]:
                dat[img][dataset]['ydat'].append(float(row[3])/(float(row[2])-float(row[1])))
            else:
                dat[img][dataset]['ydat'].append(float(row[2]))
            done = True
        if done: continue

        if row[0]=="CMD":
            dat['CMD'] = row[1]
        elif row[0]=="SN":
            dat['SN'][row[1]] = row[2]
            if len(row)==4: dat['LM'][row[1]] = row[3]
        elif row[0]=="LM":
            dat['LM'][row[1]] = row[2]
        else:
            if key not in dat: dat[key] = {}
            if 'dat' not in dat[key]:
                dat[key]['dat'] = []
            dat[key]['dat'].append(row)

    #
    # Number of reads in good fragments
    #
    nreads = 0
    if 'frag_size' in dat:
        for i in range(len(dat['frag_nreads']['xdat'])): nreads += dat['frag_nreads']['xdat'][i]*dat['frag_nreads']['ydat'][i]
    dat['SN']['n_reads_in_good_fragments'] = 2*int(nreads)
    #
    # N50
    #
    tmp = 0
    dat['SN']['N50'] = 0
    if 'frag_size' in dat:
        for i in range(len(dat['frag_size']['xdat'])): tmp += dat['frag_size']['xdat'][i]*dat['frag_size']['ydat'][i]
        tmp *= 0.5
        for i in range(len(dat['frag_size']['xdat'])): 
            tmp -= dat['frag_size']['xdat'][i]*dat['frag_size']['ydat'][i]
            if tmp > 0: continue
            dat['SN']['N50'] = int(dat['frag_size']['xdat'][i])
            break
    #
    # N10x, N20x
    #
    tmp = 0
    dat['SN']['N10x'] = 0
    dat['SN']['N20x'] = 0
    if 'frag_size' in dat:
        genome_len = float(dat['LM']['genome_length'])
        for i in range(len(dat['frag_size']['xdat'])-1,-1,-1):
            tmp += dat['frag_size']['xdat'][i]*dat['frag_size']['ydat'][i]
            if tmp / genome_len < 10: continue
            dat['SN']['N10x'] = int(dat['frag_size']['xdat'][i])
            break
        tmp = 0
        for i in range(len(dat['frag_size']['xdat'])-1,-1,-1):
            tmp += dat['frag_size']['xdat'][i]*dat['frag_size']['ydat'][i]
            if tmp / genome_len < 20: continue
            dat['SN']['N20x'] = int(dat['frag_size']['xdat'][i])
            break
    return dat

def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

def percent(part, total):
    return "%.1f%%" % (float(part)*100./float(total));

def embed_image(name):
    if name==None: return ''
    fh = open(dir+'/'+name+'.'+itype, "rb")
    return "<img class='plot' src='data:image/png;base64," + base64.b64encode(fh.read()) + "'>"

help_id = 0
def help_text(text):
    global help_id
    help_id = help_id + 1
    return \
        "<div style='float:right;padding:0.5em;'><div onclick='toggle_help(event,"+str(help_id)+")' class='help'>?</div></div>" + \
        "<div onclick='toggle_help(event,"+str(help_id)+")' id='help"+str(help_id)+"' class='help_text'>"+text+"</div>"

def write_trim(dat):
    out = []
    if dat.get('SN',{}).get('nraw_reads')==None: return out    # not an output of `trim`
    plot_err_correct(dat)
    out.append(help_text("""
                        <dl>
                            <dt>Read pairs total</dt>
                                <dd>Total number of reads, each read pair counts as one.</dd>
                            <dt>Reads with known barcodes</dt>
                                <dd>Number of reads with a known barcode</dd>
                            <dt>Reads trusted</dt>
                                <dd>Number of reads trusted because they come from a known and well populated barcode
                                with at least """ + dat['LM'].get('trust_nreads','40') + """ reads.</dd>
                            <dt>Reads modified</dt>
                                <dd>Number of reads with a modified barcode. Note that barcodes were
                                changed only if the new barcode had at least """ + dat['LM'].get('good_barcode_nreads','5') + """ reads
                                </dd>
                            <dt>Reads written</dt>
                                <dd>Total number of reads written. (Reads with unknown barcodes are discarded by default.)
                                </dd>
                            <dt>Barcoded reads written</dt>
                                <dd>Total number of barcoded reads written. (Reads with unknown barcodes are discarded by default.)
                                </dd>
                            <dt>Read pairs per barcode<dt>
                                <dd>Distribution of read pairs per known and unknown barcode, before and after error correction.
                                Modified reads show the number of read pairs in the new barcode prior to error correction.
                                </dd>
                        </dl>
                    """) +"""
                    <h2 style='margin-bottom:0;'>Fastq processing</h2>
                    <div style='margin-bottom:1em;text-align:center'>""" + dat['title'] + """</div>

                    <table class='numbers'>
                        <tr><td class='laln'>Read pairs total <td>""" 
                            + bignum(dat['SN']['nraw_reads']) + """<td></tr>
                        <tr><td class='laln'>Reads with known barcode<td>""" 
                            + bignum(dat['SN']['nraw_reads_known_bx']) + """<td>""" \
                            + percent(dat['SN']['nraw_reads_known_bx'],dat['SN']['nraw_reads']) + """</tr>
                        <tr><td class='laln'>Reads trusted<td>""" 
                            + bignum(dat['SN']['ntrusted_reads']) + """<td>""" \
                            + percent(dat['SN']['ntrusted_reads'],dat['SN']['nraw_reads']) + """</tr>
                        <tr><td class='laln'>Reads modified<td>""" 
                            + bignum(dat['SN']['nchanged_reads']) + """<td>""" \
                            + percent(dat['SN']['nchanged_reads'],dat['SN']['nraw_reads']) + """</tr>
                        <tr><td class='laln'>Reads written<td>""" 
                            + bignum(dat['SN']['nwr_reads_total']) + """<td>""" \
                            + percent(dat['SN']['nwr_reads_total'],dat['SN']['nraw_reads']) + """</tr>
                        <tr><td class='laln'>Barcoded reads written<td>""" 
                            + bignum(dat['SN']['nwr_reads_bx']) + """<td>""" \
                            + percent(dat['SN']['nwr_reads_bx'],dat['SN']['nraw_reads']) + """</tr>
                    </table>
                    <div class='sep topsep'></div> """
                    + embed_image(dat.get('img'))  + """
            """);
    return out

def write_stats(dat):
    out = []
    if dat.get('SN',{}).get('n_all_reads')==None: return out   # not an output of `stats`
    excluded = ''
    excluded_help = ''
    if 'n_excluded_soft_clips' in dat['SN'] and dat['SN']['n_excluded_soft_clips']!="0":
        excluded += """<tr><td class='lalni'>.. soft-clips&gt;"""+dat['LM']['n_excluded_soft_clips'] + """<td>""" \
                            + bignum(dat['SN']['n_excluded_soft_clips']) + """<td>""" \
                            + percent(dat['SN']['n_excluded_soft_clips'],dat['SN']['n_all_reads']) + """</tr>""";
        excluded_help += """
            <dt>Soft-clips</dt>
            <dd>Number of reads excluded because of too many soft-clipped bases.</dd>"""
    if 'n_excluded_unlisted_bx' in dat['SN'] and dat['SN']['n_excluded_unlisted_bx']!="0":
        excluded += """<tr><td class='lalni'>.. unlisted barcode<td>""" \
                            + bignum(dat['SN']['n_excluded_unlisted_bx']) + """<td>""" \
                            + percent(dat['SN']['n_excluded_unlisted_bx'],dat['SN']['n_all_reads']) + """</tr>""";
        excluded_help += """
            <dt>Unlisted barcode</dt>
            <dd>Barcode sequence not present in the list of """ + bignum(dat['LM']['n_barcodes_excluded_unlisted_bx']) + """ known barcode sequences.</dd>"""
    if 'n_excluded_anomalous_pair' in dat['SN'] and dat['SN']['n_excluded_anomalous_pair']!="0":
        excluded += """<tr><td class='lalni'>.. anomalous pair<td>""" \
                            + bignum(dat['SN']['n_excluded_anomalous_pair']) + """<td>""" \
                            + percent(dat['SN']['n_excluded_anomalous_pair'],dat['SN']['n_all_reads']) + """</tr>""";
        excluded_help += """
            <dt>Anomalous pair</dt>
            <dd>Reads from a pair mapped to different chromosomes.</dd>"""
    if 'n_excluded_bx_has_n' in dat['SN'] and dat['SN']['n_excluded_bx_has_n']!="0":
        excluded += """<tr><td class='lalni'>.. barcode has N<td>""" \
                            + bignum(dat['SN']['n_excluded_bx_has_n']) + """<td>""" \
                            + percent(dat['SN']['n_excluded_bx_has_n'],dat['SN']['n_all_reads']) + """</tr>""";
        excluded_help += """
            <dt>barcode has N</dt>
            <dd>Unknown base in the barcode sequence.</dd>"""

    filters = ''
    filters_help = ''
    if 'n_flag_DUP' in dat['SN'] and dat['SN']['n_flag_DUP']!="0":
        filters += """<tr><td class='lalnii'>.. duplicates<td>""" \
                    + bignum(dat['SN']['n_flag_DUP']) + """<td>""" \
                    + percent(dat['SN']['n_flag_DUP'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Duplicates</dt>
            <dd>Number of reads marked as duplicate.</dd>"""
    if 'n_flag_UNMAP' in dat['SN'] and dat['SN']['n_flag_UNMAP']!="0":
        filters += """<tr><td class='lalnii'>.. unmapped<td>""" \
                    + bignum(dat['SN']['n_flag_UNMAP']) + """<td>""" \
                    + percent(dat['SN']['n_flag_UNMAP'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Unmapped</dt>
            <dd>Number of unmapped reads.</dd>"""
    if 'n_flag_MUNMAP' in dat['SN'] and dat['SN']['n_flag_MUNMAP']!="0":
        filters += """<tr><td class='lalnii'>.. mate unmapped<td>""" \
                    + bignum(dat['SN']['n_flag_MUNMAP']) + """<td>""" \
                    + percent(dat['SN']['n_flag_MUNMAP'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Mate unmapped</dt>
            <dd>Number of reads with unmapped mate from the read pair.</dd>"""
    if 'n_flag_SUPPLEMENTARY' in dat['SN'] and dat['SN']['n_flag_SUPPLEMENTARY']!="0":
        filters += """<tr><td class='lalnii'>.. supplementary<td>""" \
                    + bignum(dat['SN']['n_flag_SUPPLEMENTARY']) + """<td>""" \
                    + percent(dat['SN']['n_flag_SUPPLEMENTARY'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Supplementary</dt>
            <dd>Number of supplementary reads.</dd>"""
    if 'n_flag_SECONDARY' in dat['SN'] and dat['SN']['n_flag_SECONDARY']!='0':
        filters += """<tr><td class='lalnii'>.. secondary<td>""" \
                    + bignum(dat['SN']['n_flag_SECONDARY']) + """<td>""" \
                    + percent(dat['SN']['n_flag_SECONDARY'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>Secondary</dt>
            <dd>Number of secondary reads.</dd>"""
    if 'n_flag_QCFAIL' in dat['SN'] and dat['SN']['n_flag_QCFAIL']!='0':
        filters += """<tr><td class='lalnii'>.. QC failed<td>""" \
                    + bignum(dat['SN']['n_flag_QCFAIL']) + """<td>""" \
                    + percent(dat['SN']['n_flag_QCFAIL'],dat['SN']['n_all_reads']) + """</tr>""";
        filters_help += """
            <dt>QC failed</dt>
            <dd>Number of reads not passing quality controls.</dd>"""

    frag_size_xlim = 0
    if 'frag_size' in dat: frag_size_xlim = dat['frag_size']['xlim']

    out.append(help_text("""
                        <dl>
                            <dt>Reads total</dt>
                                <dd>Total number of reads, each read pair counts as two.</dd>
                            <dt>Barcoded reads</dt>
                                <dd>Reads with the BX tag</dd>
                            <dt>Reads excluded</dt>
                                <dd>Reads excluded because of flag, low mapping quality, etc.</dd>
                            <dt>MQ&lt;"""+dat['LM']['n_excluded_mq'] + """</dt>
                                <dd>Number of reads excluded because of low mapping quality.</dd>
                            """ + excluded_help + """
                            <dt>Filtered</dt>
                                <dd>Number of reads excluded because of BAM flags.</dd>
                            """ + filters_help + """
                            <dt>Coverage</dt>
                                <dd>Good reads were filtered by flag (""" + dat['LM']['flags'] + """) and mapping quality (&ge;""" + dat['LM']['minMQ'] + """).
                                    The x axis range was set to show at least """ + str(100*dat['cov_reads']['xlim']) + """% of the data.
                                </dd>
                            <dt>Number of bases after soft-clipping</dt>
                                <dd>Histogram of read-lengths after trimming soft-clipped bases. Includes all mapped reads. </dd>
                            <dt>Unclipped reads</dt>
                                <dd>Number of mapped reads without soft-clips</dd>
                        </dl>
                    """) +"""
                    <h2>Sequencing</h2>
                    <table class='numbers'>
                        <tr><td class='laln'>Reads total <td>""" 
                            + bignum(dat['SN']['n_all_reads']) + """<td></tr>
                        <tr><td class='laln'>Barcoded reads <td>"""
                            + bignum(dat['SN']['n_bx_reads']) + """<td>"""
                            + percent(dat['SN']['n_bx_reads'],dat['SN']['n_all_reads']) +"""</tr>
                        <tr><td class='laln'>Reads excluded <td>"""
                            + bignum(dat['SN']['n_excluded']) + """<td>""" 
                            + percent(dat['SN']['n_excluded'],dat['SN']['n_all_reads']) + """</tr>
                        <tr><td class='lalni'>.. MQ&lt;"""+dat['LM']['n_excluded_mq'] + """<td>"""
                            + bignum(dat['SN']['n_excluded_mq']) + """<td>""" 
                            + percent(dat['SN']['n_excluded_mq'],dat['SN']['n_all_reads']) + """</tr>
                        """ + excluded + """
                        <tr><td class='lalni'>.. filtered<td>"""
                            + bignum(dat['SN']['n_excluded_flag']) + """<td>""" 
                            + percent(dat['SN']['n_excluded_flag'],dat['SN']['n_all_reads']) + """</tr>
                        """ + filters + """
                        <tr><td class='laln'>Unclipped reads<td>""" 
                            + bignum(dat['SN']['n_unclipped_reads']) + """<td>"""
                            + percent(dat['SN']['n_unclipped_reads'],dat['SN']['n_mapped_reads']) + """</tr>
                    </table>
                    <div class='sep topsep'></div> """
                    + embed_image(dat['cov_reads']['img'])  + """ <div class='sep'></div> """
                    + embed_image(dat['dist_sclip']['img']) + """
                """);

    out.append(help_text("""
                        <dl>
                            <dt>Unique barcodes total</dt>
                                <dd>Number of barcode sequences (the BX tag) found in the data.</dd>
                            <dt>Good GEMs</dt>
                                <dd>"GEM" stands for Gel Bead in Emulsion. In the sequencing data GEM corresponds to a set of reads with the same BX tag.
                                    Good GEMs are barcodes which were not excluded from reasons given below.
                                </dd>
                            <dt>Excluded barcodes</dt>
                                <dd>Number of excluded barcodes. Barcodes can be excluded because their sequence is not present
                                    in the list of """ + bignum(dat['LM']['n_barcodes_excluded_unlisted_bx']) + """ known barcode sequences, because
                                    there are too few read pairs with
                                    the same barcode sequence (&lt;""" + dat['LM']['n_barcodes_excluded_min_pairs'] + """)
                                    or because there is no good fragment (all reads with the same barcode map
                                    to different chromosomes).
                                </dd>
                            <dt>Reads in good GEMs</dt>
                                <dd>Number of good reads in good GEMs (but not necessarily in good fragments).
                                </dd>
                            <dt>Reads in good fragments</dt>
                                <dd>Number of good reads in good fragments.
                                At least """ + dat['LM'].get('min_readpairs_per_fragment','2') + """ read pairs per fragment are required.
                                </dd>
                            <dt>Fragment size N50</dt>
                                <dd>Shortest fragment at 50% of the total length: sort fragments by their size
                                in ascending order, mark the center of the entire length and report the size of the fragment
                                which happens to be in the middle.
                                </dd>
                            <dt>Fragment size N10x, N20x</dt>
                                <dd>Longest fragment at 10x (20x) genome coverage: sort fragments by their size
                                in descending order, mark the point where the distance from the beginning divided by the
                                genome length is at least 10 (20) and report the size of the fragment
                                which happens to be at that point. Note that the genome length is determined from the
                                BAM header as the sum of all contigs (""" +bignum(dat['LM']['genome_length'])+ """ bp).
                                </dd>
                            <dt>Read pairs per barcode</dt>
                                <dd>Read pair frequency in barcodes with this many read pairs per barcode.
                                The x axis range was set to include """ + str(100*dat['bx_reads']['xlim']) + """% of the data.
                                </dd>
                            <dt>Read pairs per fragment</dt>
                                <dd>Read pair frequency in fragments with this many read pairs per fragment.
                                The x axis range was set to include """ + str(100*dat.get('FRAG_NREADS',{}).get('xlim',0)) + """% of the data.
                                </dd>
                            <dt>Fragment coverage</dt>
                                <dd>Average depth within fragments. 
                                The x axis range was set to include """ + str(100*dat.get('FRAG_COV',{}).get('xlim',0)) + """% of the data.
                                </dd>
                        </dl>
                    """) +"""
                    <h2>GEMs</h2>
                    <table class='numbers'>
                        <tr><td class='laln'>Unique barcodes total <td>""" 
                            + bignum(dat['SN']['n_all_barcodes']) + """<td></tr>
                        <tr><td class='laln'>Good GEMs<td>""" 
                            + bignum(dat['SN']['n_good_barcodes']) + """<td>"""
                            + percent(dat['SN']['n_good_barcodes'],dat['SN']['n_all_barcodes']) +"""</tr>
                        <tr><td class='laln'>Excluded barcodes<td>"""
                            + bignum(dat['SN']['n_barcodes_excluded']) + """<td>"""
                            + percent(dat['SN']['n_barcodes_excluded'],dat['SN']['n_all_barcodes']) +"""</tr>
                        <tr><td class='lalni'>.. fewer than """ + dat['LM']['n_barcodes_excluded_min_pairs'] + """ good read pairs<td>"""
                            + bignum(dat['SN']['n_barcodes_excluded_min_pairs']) + """<td>"""
                            + percent(dat['SN']['n_barcodes_excluded_min_pairs'],dat['SN']['n_all_barcodes']) +"""</tr>
                        <tr><td class='lalni'>.. no good fragments<td>"""
                            + bignum(dat['SN']['n_barcodes_excluded_min_frags']) + """<td>"""
                            + percent(dat['SN']['n_barcodes_excluded_min_frags'],dat['SN']['n_all_barcodes']) +"""</tr>
                        <tr><td class='lalni'>.. unlisted barcode sequence<td>"""
                            + bignum(dat['SN']['n_barcodes_excluded_unlisted_bx']) + """<td>"""
                            + percent(dat['SN']['n_barcodes_excluded_unlisted_bx'],dat['SN']['n_all_barcodes']) +"""</tr>
                        <tr><td class='laln'>Reads in good GEMs<td>""" 
                            + bignum(dat['SN']['n_reads_in_good_barcodes']) + """<td>"""
                            + percent(dat['SN']['n_reads_in_good_barcodes'],dat['SN']['n_all_reads']) +"""</tr>
                        <tr><td class='laln'>Reads in good fragments<td>""" 
                            + bignum(dat['SN']['n_reads_in_good_fragments']) + """<td>"""
                            + percent(dat['SN']['n_reads_in_good_fragments'],dat['SN']['n_all_reads']) +"""</tr>
                        <tr><td class='laln'>Fragment size N50<td>"""
                            + bignum(dat['SN']['N50']) + """<td></tr>
                        <tr><td class='lalni'>.. N10x<td>"""
                            + bignum(dat['SN']['N10x']) + """<td></tr>
                        <tr><td class='lalni'>.. N20x<td>"""
                            + bignum(dat['SN']['N20x']) + """<td></tr>
                    </table>
                    <div class='sep topsep'></div> """
                    + embed_image(dat['bx_reads']['img']) + """ <div class='sep'></div> """
                    + embed_image(dat.get('frag_nreads',{}).get('img')) + """ <div class='sep'></div> """
                    + embed_image(dat.get('frag_cov',{}).get('img')) + """
            """);
    out.append(help_text("""
                        <dl>
                            <dt>Insert size</dt>
                                <dd>Insert size distribution, includes all mapped read pairs. 
                                    The x axis range was set to show at least """ + str(100*dat['insert_size']['xlim']) + """% of 
                                    the data with insert sizes bigger than """ + ("%.0f" % dat['insert_size'].get('xlim_bp',0)) + """ bp.
                                </dd>
                            <dt>Number of bases against fragment length</dt>
                                <dd>Total number of bases in fragments of given length (includes bases with zero coverage).
                                The fragment
                                size is calculated as <i>D + D/(N-1)</i>, where <i>N</i> is the number of read pairs within
                                the fragment and <i>D</i> is the distance between the first and
                                the last pair. All read pairs must map to the same chromosome with
                                the maximum gap of """ + bignum(dat['LM']['max_frag_gap']) + """ bp.
                                At least """ + dat['LM'].get('min_readpairs_per_fragment','2') + """ read pairs per fragment are required.
                                </dd>
                            <dt>Sequences bases against fragment length</dt>
                                <dd>Total number of sequence in fragments of given length.
                                </dd>
                            <dt>Number of fragments against fragment length</dt>
                                <dd>Cumulative frequency of estimated DNA fragment size.  The fragment length is calculated as above.
                                </dd>
                        </dl>
                    """)
                    + """<h2>Input DNA</h2>""" 
                    + embed_image(dat['insert_size']['img']) + """ <div class='sep'></div> """
                    + embed_image(dat.get('frag_size_bases',{}).get('img')) + """ <div class='sep'></div> """
                    + embed_image(dat.get('frag_size_seq_bases',{}).get('img')) + """ <div class='sep'></div> """
                    + embed_image(dat.get('frag_size',{}).get('img')) + """
                """);
    return out

def write_html(fname, data):
    fh = open(fname,"w")
    fh.write("""<!DOCTYPE html><html>
        <style>
            body { 
                background: #efefef; 
                text-align: center; 
                padding: 0px;
                margin: 0px;
            }
            #container { 
                background: #fff;
                display: inline-block;
                height: 100%;
                min-height: 100vh;
                padding: 0em 2em 0em 2em;
                margin: 0px;
            }
            .title {
                display: block;
                font-weight: bold;
                font-size: large;
                margin-top: 1em;
                color: """ + color[1] + """;
            }
            .column {
                margin-top:1em;
                background: #fff;
                clear: left;
                width: 490px;
                display: inline-block;
                vertical-align: top;
            }
            .box {
                padding: 0.5em;
                padding-bottom: 1em;
                margin: 0.5em;
                margin-bottom: 1em;
                border: solid 1px #ddd;
                border-radius: 0.2em;
            }
            .plot {
                width:100%;
            }
            div.sep {
                margin-top:1em;
                margin-bottom:1em;
                margin-left: auto;
                margin-right: auto;
                width: 80%; 
                height: 1px; 
                background: #ddd;
                overflow: hidden;
            }
            div.topsep {
                margin-top:2em;
            }
            h2 {
                color: #E24A33;
            }
            dt {
                font-weight: bold;
            }
            dd {
                margin-bottom: 1em;
                margin-left: 1em;
            }
            .help {
                display: inline-block;
                margin-left:0.5em;
                width: 1.5em;
                height: 1.5em;
                line-height: 1.5em;
                font-size: 0.5em;
                font-weight: bold;
                text-align: center;
                color:white;
                background-color: #ddd;
                border-radius: 50%;
                cursor: pointer;
            }
            .help_text {
                font-size: small;
                font-weight: normal;
                text-align: left;
                width: 90%;
                display: none;
                background-color: white;
                border-radius: 0.2em;
                border: solid 1px #ddd;
                cursor: pointer;
                padding: 1em;
                box-shadow: 3px 3px 3px #ddd;
            }
            .laln { text-align: left; }
            .raln { text-align: right; }
            .lalni { text-align: left; padding-left: 1.5em; }
            .lalnii { text-align: left; padding-left: 3em; }
            .tsep { padding-top:0.5em; }
            table.numbers {
                width: 90%;
                margin: auto;
                text-align: right;
            }
        </style>
        <script type="text/javascript">
            function toggle_help(e, id)
            {
                div = document.getElementById('help'+id)
                if ( div.style.display=='block' ) div.style.display='none';
                else div.style.display='block';
            }
        </script>
        <body> 
            <div id="container">
            <span class="title">"""+ main_title +"""</span>""");

    boxes1 = []
    for trim in data['trim']:
        boxes1 += write_trim(trim);
    boxes2 = write_stats(data.get('stats',{}));

    # odd boxes in the first column
    fh.write('<span class="column">')
    for box in boxes1[0:][::2]:
        fh.write('<div class="box">'+box+'</div>')
    for box in boxes2[0:][::2]:
        fh.write('<div class="box">'+box+'</div>')
    fh.write('</span>')

    # even boxes in the second column
    fh.write('<span class="column">')
    for box in boxes1[1:][::2]:
        fh.write('<div class="box">'+box+'</div>')
    for box in boxes2[1:][::2]:
        fh.write('<div class="box">'+box+'</div>')
    fh.write('</span>')

    fh.write("""
            </div>
        </body></html>
        """)
    fh.close()

#--------------------

plots = []
data = { 'trim':[] }
for fname in fnames:
    fdat = parse_file(fname)
    if fdat.get('CMD','')[0:5] == "trim ": data['trim'].append(fdat)
    else: data['stats'] = fdat

if 'stats' in data:
    dat = data['stats']
    for key in dists2:
        if key not in dat: continue
        if 'img' in dists2[key]: continue
        if key not in graphs: graphs.append(key)
        plot_dist2(dat[key])
    
    for key in dists:
        if key not in dat: continue
        if key not in graphs: graphs.append(key)
        plot_dist(dat[key])

write_html(dir+'/'+"plots.html", data)

