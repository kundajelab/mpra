from __future__ import print_function

import sys

import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
from scipy.stats import spearmanr, pearsonr

################################################################################
# scatter plots
# From David Kelley's basenji code: https://github.com/calico/basenji/blob/master/basenji/plots.py


def jointplot(vals1, vals2, out_pdf, 
              alpha=0.5, point_size=10, square=False, cor='pearsonr', 
              x_label=None, y_label=None, figsize=6, 
              sample=None, table=False,
              show=True, despine=True,
              dpi=300,
              axfont=18,
              tickfont=18,
              color='black',
              ratio=5,
              kde=True,
              bw = 'scott',
              axlim = None,
              hexbin = False,
              bincount = 100,
              title = None
              ):

    if table:
        out_txt = '%s.txt' % out_pdf[:-4]
        out_open = open(out_txt, 'w')
        for i in range(len(vals1)):
            print(vals1[i], vals2[i], file=out_open)
        out_open.close()

    if sample is not None and sample < len(vals1):
        indexes = np.random.choice(np.arange(0,len(vals1)), sample, replace=False)
        vals1 = vals1[indexes]
        vals2 = vals2[indexes]

    plt.figure(figsize=(figsize, figsize), dpi=dpi)

    if cor is None:
        cor_func = None
    elif cor.lower() in ['spearman','spearmanr']:
        cor_func = spearmanr
        corr, pval = spearmanr(vals1, vals2)
    elif cor.lower() in ['pearson','pearsonr']:
        cor_func = pearsonr
        corr, pval = pearsonr(vals1, vals2)
    else:
        cor_func = None
    
    if kde and hexbin:
        g = sns.jointplot(vals1, vals2, space=0, stat_func=None,
                          color = color,
                          kind = "hex",
                          joint_kws={'bins': 'log', 'gridsize': bincount},
                          size = figsize, ratio = ratio,
                          marginal_kws={'kde': True, 'hist': False, 'kde_kws': {'shade': True, 'bw': bw}}) 
    elif kde and not hexbin:
        g = sns.jointplot(vals1, vals2, space=0, stat_func=None,
                          color = color,
                          kind = "scatter",
                          joint_kws={'alpha':alpha, 's':point_size, 'edgecolor':'black'},
                          size = figsize, ratio = ratio,
                          marginal_kws={'kde': True, 'hist': False, 'kde_kws': {'shade': True, 'bw': bw}}) 
    else:
        g = sns.jointplot(vals1, vals2, color=color, space=0, stat_func=None, joint_kws={'alpha':alpha, 's':point_size, 'edgecolor':'black'},
                                          size = figsize, ratio = ratio)

    ax = g.ax_joint

    if square:
        if axlim == None:
            vmin, vmax = scatter_lims(vals1, vals2)
            ax.set_xlim(vmin,vmax)
            ax.set_ylim(vmin,vmax)

            ax.plot([vmin,vmax], [vmin,vmax], linestyle='--', color='black')
        else:
            ax.set_xlim(*axlim)
            ax.set_ylim(*axlim)

            ax.plot(axlim, axlim, linestyle='--', color='black') 

    else:
        xmin, xmax = scatter_lims(vals1)
        ax.set_xlim(xmin,xmax)
        ymin, ymax = scatter_lims(vals2)
        ax.set_ylim(ymin,ymax)

    if cor_func is spearmanr:
        if pval < 1e-99:
            ax.text(0.05, 0.89, 'Spearman R = %.3f\nP-value = 0' % (corr),
                ha = 'left', va = 'center',
                transform = ax.transAxes,
                fontsize=axfont)
        else:
            ax.text(0.05, 0.89, 'Spearman R = %.3f\nP-value = %.1E' % (corr, pval),
                    ha = 'left', va = 'center',
                    transform = ax.transAxes,
                    fontsize=axfont)
    if cor_func is pearsonr:
        if pval < 1e-99:
            ax.text(0.05, 0.89, 'Pearson R = %.3f\nP-value = 0' % (corr),
                ha = 'left', va = 'center',
                transform = ax.transAxes,
                fontsize=axfont)
        else:
            ax.text(0.05, 0.89, 'Pearson R = %.3f\nP-value = %.1E' % (corr, pval),
                    ha = 'left', va = 'center',
                    transform = ax.transAxes,
                    fontsize=axfont)

    if y_label is not None:
        ax.set_ylabel(y_label, fontsize = axfont, labelpad = 10)
    if x_label is not None:
        ax.set_xlabel(x_label, fontsize = axfont, labelpad = 10)
    if title is not None:
        g.fig.suptitle(title, fontsize = axfont)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(tickfont)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(tickfont) 
    
    # ax.grid(True, linestyle=':')
    # plt.tight_layout(w_pad=0, h_pad=0)
 
    if despine:
        sns.despine(offset=5, trim=True, ax=ax)
    
    plt.savefig(out_pdf, bbox_inches='tight', dpi=dpi)
    if show == True:
        plt.show()
    plt.close()


def regplot(vals1, vals2, out_pdf, poly_order=1, alpha=0.5, point_size=10, cor='pearsonr', print_sig=False, show=True, square=True, x_label=None, y_label=None, title=None, figsize=(6,6), sample=None, table=False):

    if table:
        out_txt = '%s.txt' % out_pdf[:-4]
        out_open = open(out_txt, 'w')
        for i in range(len(vals1)):
            print(vals1[i], vals2[i], file=out_open)
        out_open.close()

    if sample is not None and sample < len(vals1):
        indexes = np.random.choice(np.arange(0,len(vals1)), sample, replace=False)
        vals1 = vals1[indexes]
        vals2 = vals2[indexes]

    plt.figure(figsize=figsize)

    gold = sns.color_palette('husl',8)[1]
    ax = sns.regplot(vals1, vals2, color='black', order=poly_order, scatter_kws={'color':'black', 's':point_size, 'alpha':alpha}, line_kws={'color':gold})

    if square:
        xmin, xmax = scatter_lims(vals1, vals2)
        ymin, ymax = xmin, xmax
    else:
        xmin, xmax = scatter_lims(vals1)
        ymin, ymax = scatter_lims(vals2)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    if x_label is not None:
        ax.set_xlabel(x_label)
    if y_label is not None:
        ax.set_ylabel(y_label)

    if title is not None:
        plt.title(title)

    if cor is None:
        corr = None
    elif cor.lower() in ['spearman','spearmanr']:
        corr, csig = spearmanr(vals1, vals2)
        corr_str = 'SpearmanR: %.3f' % corr
    elif cor.lower() in ['pearson','pearsonr']:
        corr, csig = pearsonr(vals1, vals2)
        corr_str = 'PearsonR: %.3f' % corr
    else:
        corr = None

    if print_sig:
        if csig > .001:
            corr_str += '\n p %.3f' % csig
        else:
            corr_str += '\n p %.1e' % csig

    if corr is not None:
        xlim_eps = (xmax-xmin) * .03
        ylim_eps = (ymax-ymin) * .05

        ax.text(xmin+xlim_eps, ymax-3*ylim_eps, corr_str, horizontalalignment='left', fontsize=12)

    # ax.grid(True, linestyle=':')
    sns.despine()

    # plt.tight_layout(w_pad=0, h_pad=0)

    plt.savefig(out_pdf)
    if show == True:
        plt.show()
    plt.close()


def scatter_lims(vals1, vals2=None, buffer=.05):
    if vals2 is not None:
        vals = np.concatenate((vals1, vals2))
    else:
        vals = vals1
    vmin = np.nanmin(vals)
    vmax = np.nanmax(vals)

    buf = .05*(vmax-vmin)

    if vmin == 0:
        vmin -= buf/2
    else:
        vmin -= buf
    vmax += buf

    return vmin, vmax
