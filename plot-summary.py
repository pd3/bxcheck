#!/usr/bin/env python
#
#   Author: petr.danecek@sanger
#   About:  Script for creating summary from multiple bxcheck outputs
#   Usage:  
#       bxcheck trim -l known-barcodes.txt file.fq.gz -o trimmed
#       bxcheck stats -l known-barcodes.txt mapped.bam > stats.txt
#       plot-bxcheck.py trimmed.txt stats.txt -d plots
#

import sys

def usage():
    print 'Usage: plot-summary.py [OPTIONS] output.html'
    print 'Options:'
    print '   -i, --in FILE     list of bxcheck stat files, the first column is the label'
    print 'Example:'
    print '   plot-summary.py -i inputs.txt output.html'
    print ''
    sys.exit(1)


inputs = None
replot = True
output = None

if len(sys.argv) < 2: usage()
args = sys.argv[1:]
while len(args):
    if args[0]=='-i' or args[0]=='--in': 
        args = args[1:]
        inputs = args[0]
        dir  = args[0]
    elif args[0]=='--dont-replot':  # for debugging, do not call matplotlib if data file older than existing image
        replot = False
    else:
        output = args[0]
    args = args[1:]

if inputs==None: usage()
if output==None: usage()

color = [
    '#E24A33', # ggplot red
    '#777777', # ggplot gray
    '#348ABD', # ggplot blue
    '#FBC15E', # ggplot orange
    '#8EBA42', # ggplot green
]

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os, numpy, base64, re
plt.style.use('ggplot')

#-------------------------

def embed_image(name):
    if name==None: return ''
    fh = open(name, "rb")
    img = "<img class='plot' src='data:image/png;base64," + base64.b64encode(fh.read()) + "'>"
    os.remove(name)
    return img

def write_html_header(fh):
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
                width: 490px;
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
            table {
                margin: auto;
            }
            td {
                padding: 0.1em;
                padding-right: 1em;
                text-align: right;
            }
            th {
                text-align: center;
                vertical-align: top;
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
            <div id="container">""");

def write_html_footer(fh):
    fh.write("""
            </div>
        </body></html>
        """)

help_id = 0
def help_text(text):
    global help_id
    help_id = help_id + 1
    return \
        "<div style='float:right;padding:0.5em;'><div onclick='toggle_help(event,"+str(help_id)+")' class='help'>?</div></div>" + \
        "<div onclick='toggle_help(event,"+str(help_id)+")' id='help"+str(help_id)+"' class='help_text'>"+text+"</div>"

def bignum(num):
    s = str(num); out = ''; slen = len(s)
    for i in range(slen):
        out += s[i]
        if i+1<slen and (slen-i-1)%3==0: out += ','
    return out

def percent(part, total):
    return "%.1f%%" % (float(part)*100./float(total));

def color_values(dat,desc={}):
    keys = dat.keys()
    for key in keys:
        min = max = dat[key][0]
        for val in dat[key]:
            if min > val: min = val
            if max < val: max = val
        td = 'td_'+key
        dat[td] = []
        for i in range(len(dat[key])):
            val = (dat[key][i]-float(min))/(max - float(min))
            if key in desc: val = 1 - val
            if val <= 0.5:
                col = 'rgb(226,74,51,%.1f)' % (2*0.8*(0.5-val) + 0.1)  # red
            else:
                col = 'rgb(142,186,66,%.1f)' % (2*0.9*(val-0.5) + 0.1) # green
            dat[td].append('<td style="background-color:'+col+'">')

def write_help_text(fh):
    fh.write(help_text("""
                        <dl>
                            <dt>Total reads</dt>
                                <dd>The total number of reads in all fastq files</dd>
                            <dt>Barcoded reads</dt>
                                <dd>Number of reads with barcodes written by the <tt>bxcheck trim</tt> command</dd>
                            <dt>Bad reads</dt>
                                <dd>Reads excluded for various reasons (low mapping quality, flag, mate mapped to 
                                a different chromosome)
                                </dd>
                            <dt>Reads in good fragments</dt>
                                <dd>Number of good reads in good fragments
                                </dd>
                            <dt>N50</dt>
                                <dd>Shortest fragment at 50% of the total length: sort fragments by their size
                                    in ascending order, mark the center of the entire length and report the size of the fragment
                                    which happens to be in the middle.
                                </dd>
                            <dt>N10x</dt>
                                <dd>Longest fragment at 10x genome coverage: sort fragments by their size in
                                descending order, mark the point where the distance from the beginning divided by the
                                genome length is at least 10 and report the size of the fragment which happens to be at
                                that point. Note that the genome length is determined from the BAM header as the sum of
                                all contigs.
                                </dd>
                            <dt>Fragment length against Nx coverage</dt>
                                <dd>N10x fragment length calculated for all dephts.
                                </dd>
                            <dt>Number of fragments (density) againsts fragment length</dt>
                                <dd>Fragment size distribution plotted as density. 
                                The fragment size is calculated as <i>D + D/(N-1)</i>, where <i>N</i> is the number of read pairs
                                within the fragment and <i>D</i> is the distance between the first and the last pair. All read pairs
                                must map to the same chromosome with the maximum gap of """ + bignum(dat[dat.keys()[0]].get('max_frag_gap','100000')) + """ bp.
                                Only fragments with more than """+ dat[dat.keys()[0]].get('min_readpairs_per_fragment','2') +""" read pairs are included.
                                </dd>
                        </dl>
                    """))
    fh.write("<div class='topsep'></div>")

def calc_N(dat,percent):
    tmp = 0
    for i in range(len(dat['xdat'])):
        tmp += dat['xdat'][i]*dat['ydat'][i]
    tmp *= 1 - percent/100.
    for i in range(len(dat['xdat'])):
        tmp -= dat['xdat'][i]*dat['ydat'][i]
        if tmp <= 0: return dat['xdat'][i]
    return 0

def calc_Nx(dat,genome_len):
    xdat = []
    ydat = []
    sum = 0
    for i in range(len(dat['xdat'])-1,-1,-1):
        sum += dat['xdat'][i]*dat['ydat'][i]
        xdat.append(float(sum)/genome_len)
        ydat.append(dat['xdat'][i])
    return {'xdat':xdat,'ydat':ydat}

def get_Nx(dat,coverage):
    for i in range(len(dat['xdat'])):
        if dat['xdat'][i] >= coverage: return dat['ydat'][i]
    return 0

def ymax(xdat,ydat,xlim):
    ymax = 0
    for i in range(len(xdat)):
        if xlim[0] > xdat[i]: continue
        if xlim[1]!=None and xlim[1] < xdat[i]: continue
        ymax = max(ymax,ydat[i])
    return ymax

def smooth_dist(dist):
    xdat = []
    ydat = []
    n = 20
    xbuf = []
    ybuf = []
    sum  = 0
    for i in range(len(dist['xdat'])):
        if len(xbuf)==n:
            xdat.append(xbuf[0])
            ydat.append(sum/n)
            sum -= ybuf[0]
            xbuf = xbuf[1:]
            ybuf = ybuf[1:]
        xbuf.append(dist['xdat'][i])
        ybuf.append(dist['ydat'][i])
        sum += dist['ydat'][i]
    dist['xdat'] = xdat
    dist['ydat'] = ydat

def plot_dist(dat, names, args):
    wh = (7,3)
    fig, ax1 = plt.subplots(1, 1, figsize=wh)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()
    ax1.spines['bottom'].set_color('#aaaaaa')
    ax1.spines['left'].set_color('#aaaaaa')

    ylim = 0
    alpha = 0.8
    for i in range(len(names)):
        if 'handler' in args: globals()[args['handler']](dat[i])
        ax1.plot(dat[i]['xdat'],dat[i]['ydat'],'.-',alpha=alpha,label=names[i])
        if 'xlim' in args:
            ylim = max(ylim, ymax(dat[i]['xdat'],dat[i]['ydat'],args['xlim']))
    if 'xlim' in args: ax1.set_ylim(0,ylim*1.05)
    ax1.ticklabel_format(style='sci', scilimits=(-2,2), axis='y')
    if 'ylog' in args: ax1.set_yscale('symlog')
    if 'xlog' in args: ax1.set_xscale('symlog')
    if 'xlim' in args: ax1.set_xlim(args['xlim'])
    if 'ylabel' in args: ax1.set_ylabel(args['ylabel'])
    if 'xlabel' in args: ax1.set_xlabel(args['xlabel'])
    plt.legend(numpoints=1,markerscale=1,loc='best',prop={'size':10},frameon=False)

    plt.subplots_adjust(bottom=0.18)
    plt.savefig(output+'.png')
    plt.close()
    return output+'.png'

def write_html(fname,dat):
    fh = open(fname,"w")
    write_html_header(fh)
    write_help_text(fh)

    names = sorted(dat.keys())
    out = {}
    out['raw'] = []
    out['bx']  = []
    out['bad'] = []
    out['nfrag'] = []
    out['n10x'] = []
    out['n50'] = []
    for name in names: 
        out['raw'].append(dat[name]['nraw_reads'])
        out['bx'].append(float(dat[name]['nwr_reads_bx'])/float(dat[name]['nraw_reads']))
        out['bad'].append(float(dat[name]['n_excluded'])/float(dat[name]['n_all_reads']))
        out['nfrag'].append(float(dat[name]['n_reads_in_good_fragments'])/float(dat[name]['n_all_reads']))
        dat[name]['Nx'] = calc_Nx(dat[name]['frag_size'],dat[name]['genome_length'])
        out['n10x'].append(get_Nx(dat[name]['Nx'],10))
        out['n50'].append(calc_N(dat[name]['frag_size'],50))
    color_values(out,{'bad':-1})
    fh.write('''<table><tr><td>
        <th>Total<br>reads
        <th>Barcoded<br>reads
        <th>Bad<br>reads
        <th>Reads in<br>good fragments
        <th>N50
        <th>N10x
        ''')
    for i in range(len(names)):
        name = names[i]
        fh.write('<tr><td>'+name
            +out['td_raw'][i]+bignum(dat[name]['nraw_reads'])
            +out['td_bx'][i]+percent(dat[name]['nwr_reads_bx'],dat[name]['nraw_reads'])
            +out['td_bad'][i]+percent(dat[name]['n_excluded'],dat[name]['n_all_reads'])
            +out['td_nfrag'][i]+percent(dat[name]['n_reads_in_good_fragments'],dat[name]['n_all_reads'])
            +out['td_n50'][i]+bignum(out['n50'][i])
            +out['td_n10x'][i]+bignum(out['n10x'][i])
            )
    fh.write('</table>')

    img = plot_dist([dat[name]['Nx'] for name in names],names,{'xlabel':'Nx coverage','ylabel':'Fragment length','xlim':[1,60]});
    fh.write(""" <div class='topsep sep'></div> """ + embed_image(img));

    img = plot_dist([dat[name]['frag_size_density'] for name in names],names,{'xlabel':'Fragment length','ylabel':"Number of fragments\n(density)",'xlog':1,'ylog':1,'xlim':[100,None],'handler':'smooth_dist'});
    fh.write(""" <div class='sep'></div> """ + embed_image(img));

    write_html_footer(fh)
    fh.close()
 
def set(dat,name,key,value):
    if name not in dat: dat[name] = {}
    if key not in dat[name]: dat[name][key] = value

def addto(dat,name,key,value):
    if name not in dat: dat[name] = {}
    if key not in dat[name]: dat[name][key] = 0
    dat[name][key] += value

def appendto(dat,name,key,xval,yval):
    if name not in dat: dat[name] = {}
    if key not in dat[name]: 
        dat[name][key] = {}
        dat[name][key]['xdat'] = []
        dat[name][key]['ydat'] = []
    dat[name][key]['xdat'].append(xval)
    dat[name][key]['ydat'].append(yval)

def parse_file(dat,line):
    x = re.split(r'\s+', line)
    name = x[0]
    file = x[1]
    fh = open(file,'r')
    for dat_line in fh:
        row = re.split(r'\t', dat_line.rstrip('\n'))
        if row[0]=='SN' and row[1]=='nraw_reads': addto(dat,name,'nraw_reads',int(row[2])); continue
        if row[0]=='SN' and row[1]=='nwr_reads_bx': addto(dat,name,'nwr_reads_bx',int(row[2])); continue
        if row[0]=='SN' and row[1]=='n_excluded': addto(dat,name,'n_excluded',int(row[2])); continue
        if row[0]=='SN' and row[1]=='n_all_reads': addto(dat,name,'n_all_reads',int(row[2])); continue
        if row[0]=='FRAG_NREADS': addto(dat,name,'n_reads_in_good_fragments',2*int(row[1])*int(row[3])); continue
        if row[0]=='FRAG_SIZE':
            appendto(dat,name,'frag_size',int(row[1]),int(row[3]))
            appendto(dat,name,'frag_size_density',int(row[1]),float(row[3])/(float(row[2])-float(row[1])))
            continue
        if row[0]=='LM' and row[1]=='genome_length': set(dat,name,'genome_length',int(row[2])); continue
        if row[0]=='LM' and row[1]=='max_frag_gap': set(dat,name,'max_frag_gap',row[2]); continue
        if row[0]=='LM' and row[1]=='min_readpairs_per_fragment': set(dat,name,'min_readpairs_per_fragment',row[2]); continue


#-------------------------

dat = {}
lines = open(inputs, 'rb').readlines()
for input in lines:
    parse_file(dat,input)

write_html(output,dat)

