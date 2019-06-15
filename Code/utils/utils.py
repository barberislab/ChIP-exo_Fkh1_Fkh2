import pandas as pd
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np


def _color_red_or_green(val):
    """Function used to color pandas dataframes. Red < 0, green >= 0."""
    if val < 0:
        color = 'red'
    else: 
        color = 'green'
    return 'color: %s' % color

def get_functional_enrichment(df, col, df_in=None):
    """df is the complete dataset. col the column to compare. Thresholds the SNR values
    to use as a cutoff. df_in, optional, is a subset of the indices to use instead of looking
    at the thresholds.
    By default, without df_in, the function returns enrichment w.r.t. the measured genes for all 4 of our
    experiments. With df_in, it returns enrichtment of the subset in df_in w.r.t. the entire df.
    """
    
    df_go_terms = pd.DataFrame()

    # genome-wide
    total_count = df[col].value_counts() / len(df) * 100
    df_go_terms['Genome-wide %'] = total_count
    
    
    if df_in is None:
        df_in = df.copy()
        
        # measured only
        df_measured = df_in[(df_in['Fkh1 log']>=0) | (df_in['Fkh1 stat']>=0) | (df_in['Fkh2 log']>=0) | (df_in['Fkh2 stat']>=0) ]
        measured_count = df_measured[col].value_counts() / len(df_measured) * 100
        df_go_terms['Measured %'] = measured_count

        for exp in ['Fkh1 log','Fkh1 stat','Fkh2 log','Fkh2 stat']:
            # filter out targets
            df_selected = df_in[df_in[exp]>=1]

            # calculate functional enrichment w.r.t. the measured genes
            selected_count = df_selected[col].value_counts() / len(df_selected) * 100
            perc_selected = selected_count - measured_count
            df_go_terms[exp + ' %'] = perc_selected   
    else:            
        df_selected = df.loc[df_in]        
        
        selected_count = df_selected[col].value_counts() / len(df_selected) * 100
        df_go_terms['Verified %'] = selected_count
        df_go_terms['Enrichment %'] =  selected_count - df_go_terms['Genome-wide %'] 
    
    df_go_terms = df_go_terms.round(decimals=2)
        
    df_go_terms = df_go_terms.fillna(0)

    return df_go_terms


def get_functional_enrichment_multi(d,col):
    """d is a dictionary of dataframes that are already filtered for targets. The keys are either 'all' or indicate the experiment.
    col indicates a column like: primary GO term.
    """
    
    df_go_terms = pd.DataFrame()

    # genome-wide
    total_count = d['all'][col].value_counts() / len(d['all']) * 100
    df_go_terms['Genome-wide %'] = total_count

    for exp in ['Fkh1 log','Fkh1 stat','Fkh2 log','Fkh2 stat']:
        # filter out targets
        df_selected = d[exp]

        # calculate functional enrichment w.r.t. the measured genes
        selected_count = df_selected[col].value_counts() / len(df_selected) * 100
        perc_selected = selected_count - total_count
        df_go_terms[exp + ' %'] = perc_selected   
    
    df_go_terms = df_go_terms.round(decimals=2)
        
    df_go_terms = df_go_terms.fillna(0)

    return df_go_terms


def draw_stackplot(d, thresholds, exp, ax, no_xlabel, pos_mult = (1,1)):

    t = thresholds # shorthand

    ########
    # Assign an x-axis position for each phase's stacked gene column 
    # (roughly the middle unless there are multiple groups in one phase)
    ########
    # Phase assignment from the publication
    # "G1(P) Pre-replicative late G1 -10 38"
    # "G1/S G1/S 38 55"
    # "S S 55 89"
    # "G2 G2 89 130"
    # "G2/M G2/M 130 145"
    # "M Mitosis 145 180"
    # "M/G1 early G1 180 195"
    # "G1 G1 195 290"
    phase_map = {'G1(P)':1,
                 'G1/S':2,
                 'S':3,
                 'G2':4,
                 'G2/M':5,
                 'M':6,
                 'M/G1':7, 
                 'G1':8}
    
    # have to make sure d is sorted ascending on SNR: lowest SNR first
    SNR_col = 'maxPeak_AD_12 '+exp
    d = d.sort_values(by=SNR_col,ascending=True)

    # the x values will be determined by the phase or sub-phase each gene peaks in
    x = d['Expression peak phase'].values
    x = [phase_map[g] for g in x]
    # the y values are based on the ranking in SNR
    y = d[SNR_col].values
    y = []
    start_height = 0.01
    step_size = 0.3
    empty_space = 0.05
    phase_stacker = {k:start_height for k in phase_map.values()}   
    
    for xi in x:
        y.append(phase_stacker[xi]+step_size)
        phase_stacker[xi] += step_size+empty_space
        
    names = d['Standard name'].values
    enz = d['is enzyme'].values
    
    colors_dict = {"Cell cycle":"#2ca02c" ,"Cell division":"#ffe119", 
                   "DNA replication":"#0080ff", "Signal transduction":"#cc33cc",
                   "Metabolism":"#ff7f0e","None":"#F8F8FF"}
    colors = [colors_dict[v] for v in d['Primary GO term'].values]
    ecolors = ['black'] * len(d)
    
    lw = []; ls = []; alpha = []; star = []; triangle = [];
    for i in range(len(d)):
        row = d.iloc[i]
        
        if (row['maxPeak_AD_12 '+exp] >= t[0] and row['GEM '+exp] >= t[1] and row['MACE '+exp] <= t[2]):
            lw.append(1.0)
            b = 'dashed'
        else: 
            lw.append(1.0)
            b = 'solid'
        
        ls.append(b)
        
        times_verified = sum([row['MacIsaac 2006 ' + exp[:4]],row['Venters 2011 ' + exp[:4]],row['Ostrow 2014 ' + exp[:4]]])
        if times_verified == 0: # Novel = triangle
            alpha.append(0.5)
            star.append(False)
            triangle.append(True)
        elif times_verified == 1: 
            alpha.append(0.5)
            star.append(False)
            triangle.append(False)
        elif times_verified == 2: 
            alpha.append(0.5)
            star.append(False)
            triangle.append(False)
        else:
            alpha.append(0.5) # 4x verified = star
            star.append(True)
            triangle.append(False)
        
    w = ['bold' if enz[i] else 'normal' for i in range(len(d))]

    
    ### ACTUALLY BUILD THE STACKPLOT
    for i, txt in enumerate(names):
            
        size = [0.8, 0.3]
        square = Rectangle((x[i]-size[0]/2,y[i]-size[1]/2), size[0], size[1], facecolor=colors[i],
                           edgecolor=ecolors[i], alpha=alpha[i], linewidth=lw[i], linestyle=ls[i]);
        ax.add_patch(square);
        
        # x position changes with length of gene name and whether its a (bold lettered) enzyme
        # y position changes only with enzyme status
        # the multipliers allow fine control in individual images
        if len(txt) == 4:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.18*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.08*(1+0.12*enz[i])), weight = w[i]);
        elif len(txt) == 5:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.25*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.08*(1+0.12*enz[i])), weight = w[i]);
        elif len(txt) == 6:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.28*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.08*(1+0.12*enz[i])), weight = w[i]);
        elif len(txt) == 7:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.34*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.08*(1+0.12*enz[i])), weight = w[i]);
        elif len(txt) > 7:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.34*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.08*(1+0.12*enz[i])), weight = w[i]);
        else:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.15*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.08*(1+0.12*enz[i])), weight = w[i]);
        
        if star[i]:
            ax.plot(x[i]+size[0]/2.+0.08, y[i]-0.01*(1+0.12*enz[i]),'k*');
        if triangle[i]:
            ax.plot(x[i]+size[0]/2.+0.08, y[i]-0.01*(1+0.12*enz[i]),'k^');



    ax.set_ylim([0,np.max(y)+size[1]*1.1])
    ax.set_xlim([0.5,8.55]) # a bit extra space on right for '**' signifiers
    ax.set_yticks([])
    ax.set_ylabel('SNR')
    
    if not no_xlabel:
        ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8])
        ax.set_xticklabels(['G1(P)','G1/S','S','G2','G2/M','M','M/G1', 'G1'])
        ax.set_xlabel("Cell cycle phase")
    else:
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off')

    return ax


def plot_rank_displacement_SNR_foldchange(df, A, B, strA, strB, filename, legend_loc):
    """Compares rank displacement and SNR fold change between columns A and B in df. 
    strA and strB are used in the xlabels and legend to identify the columns. 
    filename indicates the image filename."""
    
    rank_disp = [] # sorted order (rank)
    val_disp = [] # SNR
    z_to_nz = []
    nz_to_z = []

    dfA = df.sort_values([A,'Standard name'],ascending=False)
    dfB = df.sort_values([B,'Standard name'],ascending=False)

    # loop for various number of genes
    # only look at the first x genes to find a trend if any exists
    for i,case in enumerate([100, len(dfA[dfA[A] >= 1]), len(dfA.index)]):
        dfA_sub = dfA.iloc[:case] 

        rank_disp.append({})
        val_disp.append({})
        z_to_nz.append([])
        nz_to_z.append([])

        for g in dfA_sub.index:
            # do not consider missing experimental values
            if not np.isnan(dfA_sub.loc[g][A]): 

                # extract rank and SNR 
                rA = dfA.index.get_loc(g)
                rB = dfB.index.get_loc(g)
                snrA = dfA.loc[g][A]
                snrB = dfB.loc[g][B]

                # only calculate index displacement and fold change for non-zeros in both datasets
                # to avoid infinities 
                if snrA == 0 and snrB != 0:
                    z_to_nz[i].append((g, snrA, snrB))
                elif snrB == 0 and snrA != 0:
                    nz_to_z[i].append((g, snrA, snrB))
                elif snrA == 0 and snrB == 0:
                    pass
                else:
                    rank_disp[i][g] = rA - rB
                    val_disp[i][g] = np.log10(snrB/snrA)
            else:
                pass


    print([len(nz_to_z[i]) for i in range(len(nz_to_z))],'genes went from non-zero to zero')    
    print([len(z_to_nz[i]) for i in range(len(z_to_nz))],'genes went from zero to non-zero')    

    hist_data = [list(d.values()) for d in rank_disp]
    box_data = [list(d.values()) for d in val_disp]
    lengths = [len(box_data[i]) for i in range(len(rank_disp))]

    fig = plt.figure()
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    ax2 = plt.subplot2grid((1, 2), (0, 1))

    ax1.boxplot(box_data)
    ax1.set_ylabel(r'$\log_{10}(\frac{SNR_{'+strB+'}}{SNR_{'+strA+'}})$')
    ax1.set_xticklabels(['Top 100','SNR$_{' + strA + '} \geq 1$ (' + str(lengths[1]) + ')','SNR$_{' + strA + '} > 0$ (' + str(lengths[-1]) + ')'])

    ax2.hist(hist_data,histtype='step', bins=150)
    limits = [min([min(hist_data[i]) for i in range(len(hist_data))]),max([max(hist_data[i]) for i in range(len(hist_data))])]
    if abs(limits[0]) > 2000 or abs(limits[1]) > 2000:
        seq = [n for n in range(round(limits[0]/2000.0)*2000,round(limits[1]/2000.0)*2000,2000)]
    else:
        seq = [n for n in range(round(limits[0]/1000.0)*1000,round(limits[1]/1000.0)*1000,1000)]
    
    ax2.set_xticks(seq)
    ax2.set_xlabel(r'Rank displacement: $rank_{' + strA + '} - rank_{' + strB + '}$')
    ax2.legend(['SNR$_{' + strA + '} > 0$ (' + str(lengths[-1]) + ')','SNR$_{' + strA + '} \geq 1$ (' + str(lengths[1]) + ')','Top 100'], loc=legend_loc)

    
    fig.set_size_inches(19,6)
    plt.show()
    
    if filename != '':
        fig.savefig(filename, bbox_inches='tight')
        

        
def draw_stackplot_cc(d, exp, ax, no_xlabel, pos_mult = (1,1)):

    ########
    # Assign an x-axis position for each phase's stacked gene column 
    # (roughly the middle unless there are multiple groups in one phase)
    ########
    # Phase assignment from the publication
    # "G1(P) Pre-replicative late G1 -10 38"
    # "G1/S G1/S 38 55"
    # "S S 55 89"
    # "G2 G2 89 130"
    # "G2/M G2/M 130 145"
    # "M Mitosis 145 180"
    # "M/G1 early G1 180 195"
    # "G1 G1 195 290"
    phase_map = {'G1(P)':1,
                 'G1/S':2,
                 'S':3,
                 'G2':4,
                 'G2/M':5,
                 'M':6,
                 'M/G1':7, 
                 'G1':8}
    
    # have to make sure d is sorted ascending on SNR: lowest SNR first
    d = d.sort_values(by=['Standard name'], ascending=False)

    # the x values will be determined by the phase or sub-phase each gene peaks in
    x = d['Expression peak phase'].values
    x = [phase_map[g] for g in x]

    y = []
    start_height = 0.01
    step_size = 0.3
    empty_space = 0.05
    phase_stacker = {k:start_height for k in phase_map.values()}   
    
    for xi in x:
        y.append(phase_stacker[xi]+step_size)
        phase_stacker[xi] += step_size+empty_space
        
    names = d['Standard name'].values
    enz = d['is enzyme'].values
    
    colors_dict = {"Cell cycle":"#2ca02c" ,"Cell division":"#ffe119", 
                   "DNA replication":"#0080ff", "Signal transduction":"#cc33cc",
                   "Metabolism":"#ff7f0e","None":"#F8F8FF"}
    colors = [colors_dict[v] for v in d['Primary GO term'].values]
    ecolors = ['black'] * len(d)
    
    lw = []; ls = []; alpha = []; colorF1 = []; colorF2 = []
    alpha = [0.5]*len(d)
    lw = [1.0]*len(d)
    ls = ['solid']*len(d)
    
    for i in range(len(d)):
        row = d.iloc[i]
        s = row['Target of']

        if 'F1L' in s and 'F1S' in s:
            colorF1.append('#e52521')
        elif 'F1L' in s:
            colorF1.append('green')
        elif 'F1S' in s:
            colorF1.append('yellow')
        else:
            colorF1.append('white')
            
        if 'F2L' in s and 'F2S' in s:
            colorF2.append('#e52521')
        elif 'F2L' in s:
            colorF2.append('green')
        elif 'F2S' in s:
            colorF2.append('yellow')
        else:
            colorF2.append('white')
        
    w = ['bold' if enz[i] else 'normal' for i in range(len(d))]

    
    ### ACTUALLY BUILD THE STACKPLOT
    for i, txt in enumerate(names):
            
        size = [0.75, 0.3]
        square = Rectangle((x[i]-size[0]/2,y[i]-size[1]/2), size[0], size[1], facecolor=colors[i],
                           edgecolor=ecolors[i], alpha=alpha[i], linewidth=lw[i], linestyle=ls[i]);
        ax.add_patch(square);
        
        squareF2 = Rectangle((x[i]+size[0]/1.85,y[i]-0.1), 0.1, 0.1 , facecolor=colorF2[i],
                           edgecolor=ecolors[i], alpha=1, linewidth=1, linestyle='solid');
        ax.add_patch(squareF2);
        
        squareF1 = Rectangle((x[i]+size[0]/1.85,y[i]+size[1]/100), 0.1, 0.1, facecolor=colorF1[i],
                           edgecolor=ecolors[i], alpha=1, linewidth=1, linestyle='solid');
        ax.add_patch(squareF1);
        
        if len(txt) == 4:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.21*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.07*(1+0.2*enz[i])), weight = w[i]);
        elif len(txt) == 5:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.25*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.07*(1+0.2*enz[i])), weight = w[i]);
        elif len(txt) == 6:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.28*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.07*(1+0.2*enz[i])), weight = w[i]);
        elif len(txt) == 7:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.36*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.07*(1+0.2*enz[i])), weight = w[i]);
        elif len(txt) > 7:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.39*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.07*(1+0.2*enz[i])), weight = w[i]);
        else:
            ax.annotate(txt, (x[i]-pos_mult[0]*0.15*(1+0.1*enz[i]),y[i]-pos_mult[1]*0.07*(1+0.2*enz[i])), weight = w[i]);



    ax.set_ylim([0,np.max(y)+size[1]*1.1])
    ax.set_xlim([0.5,8.55]) # a bit extra space on right for '**' signifiers
    ax.set_yticks([])
    ax.set_ylabel('')
    
    if not no_xlabel:
        ax.set_xticks([1, 2, 3, 4, 5, 6, 7, 8])
        ax.set_xticklabels(['G1(P)','G1/S','S','G2','G2/M','M','M/G1', 'G1'])
        ax.set_xlabel("Cell cycle phase")
    else:
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off')

    return ax


def identify_papers_showing_target(row, case):
    """
    return string that contains the papers showing this gene as a target.
    case identifies the transcription factor.
    """

    s = 'Mondeel'
    
    if row['MacIsaac 2006 ' + case[:4]]: # computational
        s += ', MacIsaac'
    if row['Venters 2011 ' + case[:4]]:
        s += ', Venters'
    if row['Ostrow 2014 ' + case[:4]]:
        s += ', Ostrow'
        
#    if len(s) > 0:
#        s = s[2:]
        
    return s