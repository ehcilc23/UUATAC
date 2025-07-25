'''Motif Visualization in NvTK.

Currently, this module only support DNA MOTIF Visualization.

Protein Sequence Motif Visualization was under development.
'''

# Modified motif visualization functions from DeepOmic
# Jiaqili@zju.edu.cn

__all__ = ["filter_heatmap", "plot_filter_heatmap", "plot_filter_logos", "seq_logo", "plot_seq_logo"]

import os, logging
import numpy as np
import pandas as pd

from PIL import Image

import matplotlib
import matplotlib.image as mpimg
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .Motif import normalize_pwm

# def imresize(arr, shape):
#     return np.array(Image.fromarray(arr).resize((shape[1], shape[0])))

def filter_heatmap(pwm, output_fname=None, save=False, fig_size=(10, 7), 
                    norm=True, cmap='hot_r', cbar_norm=True):
    """visualize pwm represented Filter in heatmap."""
    pwm_dim, pwm_len = pwm.shape

    plt.figure(figsize=fig_size)
    if norm:
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    else:
        norm = None
    cmap_reversed = matplotlib.cm.get_cmap(cmap)
    im = plt.imshow(pwm, cmap=cmap_reversed, norm=norm, aspect="auto")

    #plt.axis('off')
    ax = plt.gca()
    ax.set_xticks(np.arange(-.5, pwm_len, 1.), minor=True)
    ax.set_yticks(np.arange(-.5, pwm_dim, 1.), minor=True)
    ax.grid(which='minor', color='k', linestyle='-', linewidth=2)
    plt.xticks(list(range(pwm_len)))
    if pwm.shape[0] == 4:
        plt.yticks([0, 1, 2, 3], ['A', 'C', 'G', 'T'], fontsize=16)
    else:
        plt.yticks(list(range(pwm_dim)), list(range(pwm_dim)), fontsize=16)

    #cbar = plt.colorbar()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=16)
    if cbar_norm:
        cbar.set_ticks([0.0, 0.5, 1.0])

    if save:
        plt.savefig(output_fname, format="pdf")
    plt.show()
    plt.close()


def plot_filter_heatmap(W, factor=5, fig_size=(10,7), save=True):
    """visualize pwm represented Filter in heatmap."""
    for idx, pwm in enumerate(W):
        output_fname = "Motif_" + str(idx) +".pdf"
        pwm = normalize_pwm(pwm, factor=factor)
        filter_heatmap(pwm, output_fname=output_fname, save=save, fig_size=fig_size)


def plot_filter_logos(W, figsize=(10,7), height=25, nt_width=10, norm=0, alphabet='dna', norm_factor=3, num_rows=None):
    """visualize pwm represented filter Logo."""
    # W = np.squeeze(W.transpose([3, 2, 0, 1]))
    num_filters = W.shape[0]
    if not num_rows:
        num_rows = int(np.ceil(np.sqrt(num_filters)))
        num_cols = num_rows
    else:
        num_cols = int(np.ceil(num_filters//num_rows))
    grid = matplotlib.gridspec.GridSpec(num_rows, num_cols)
    grid.update(wspace=0.2, hspace=0.2, left=0.1, right=0.2, bottom=0.1, top=0.2)
    fig = plt.figure(figsize=figsize);
    if norm:
        MAX = np.max(W)
    else:
        MAX = None

    for i in range(num_filters):
        plt.subplot(grid[i]);
        if norm_factor:
            W_norm = normalize_pwm(W[i], factor=norm_factor, max=MAX)
        else:
            W_norm = W[i]
        logo = seq_logo(W_norm, height=height, nt_width=nt_width, norm=0, alphabet=alphabet)
        plot_seq_logo(logo, nt_width=nt_width, step_multiple=None)
        #if np.mod(i, num_rows) != 0:
        plt.yticks([])
    return fig


# help functions
def plot_seq_logo(logo, nt_width=None, step_multiple=None):
    """visualize sequence logo."""
    plt.imshow(logo, interpolation='none')
    if nt_width:
        num_nt = logo.shape[1]/nt_width
        if step_multiple:
            step_size = int(num_nt/(step_multiple+1))
            nt_range = range(step_size, step_size*step_multiple)
            plt.xticks([step_size*nt_width, step_size*2*nt_width, step_size*3*nt_width, step_size*4*nt_width],
                        [str(step_size), str(step_size*2), str(step_size*3), str(step_size*4)])
        else:
            plt.xticks([])
        plt.yticks([0, 50], ['2.0','0.0'])
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('none')
    else:
        plt.imshow(logo, interpolation='none')
        plt.axis('off');


def load_alphabet(char_path, alphabet, colormap='standard'):

    def load_char(char_path, char, color):
        colors = {}
        colors['green'] = [10, 151, 21]
        colors['red'] = [204, 0, 0]
        colors['orange'] = [255, 153, 51]
        colors['blue'] = [0, 0, 204]
        colors['cyan'] = [153, 204, 255]
        colors['purple'] = [178, 102, 255]
        colors['grey'] = [160, 160, 160]
        colors['black'] = [0, 0, 0]

        img = mpimg.imread(os.path.join(char_path, char+'.eps'))
        img = np.mean(img, axis=2)
        x_index, y_index = np.where(img != 255)
        y = np.ones((img.shape[0], img.shape[1], 3))*255
        for i in range(3):
            y[x_index, y_index, i] = colors[color][i]
        return y.astype(np.uint8)


    colors = ['green', 'blue', 'orange', 'red']
    if alphabet == 'dna':
        letters = 'ACGT'
        if colormap == 'standard':
            colors = ['green', 'blue', 'orange', 'red']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))

    elif alphabet == 'rna':
        letters = 'ACGU'
        if colormap == 'standard':
            colors = ['green', 'blue', 'orange', 'red']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))


    elif alphabet == 'structure': # structural profile

        letters = 'PHIME'
        if colormap == 'standard':
            colors = ['blue', 'green', 'orange', 'red', 'cyan']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))

    elif alphabet == 'pu': # structural profile

        letters = 'PU'
        if colormap == 'standard':
            colors = ['cyan', 'purple']
        elif colormap == 'bw':
            colors = ['black', 'grey']
        chars = []
        for i, char in enumerate(letters):
            chars.append(load_char(char_path, char, colors[i]))

    return chars


def seq_logo(pwm, height=30, nt_width=10, norm=0, alphabet='dna', colormap='standard'):
    """generate sequence logo from pwm."""
    def get_nt_height(pwm, height, norm):

        def entropy(p):
            s = 0
            for i in range(len(p)):
                if p[i] > 0:
                    s -= p[i]*np.log2(p[i])
            return s

        num_nt, num_seq = pwm.shape
        heights = np.zeros((num_nt,num_seq));
        for i in range(num_seq):
            if norm == 1:
                total_height = height
            else:
                total_height = (np.log2(num_nt) - entropy(pwm[:, i]))*height;
            if alphabet == 'pu':
                heights[:,i] = np.floor(pwm[:,i]*np.minimum(total_height, height));
            else:
                heights[:,i] = np.floor(pwm[:,i]*np.minimum(total_height, height*2));

        return heights.astype(int)


    # get the alphabet images of each nucleotide
    package_directory = os.path.dirname(os.path.abspath(__file__))
    char_path = os.path.join(package_directory,'chars')
    chars = load_alphabet(char_path, alphabet, colormap)

    # get the heights of each nucleotide
    heights = get_nt_height(pwm, height, norm)

    # resize nucleotide images for each base of sequence and stack
    num_nt, num_seq = pwm.shape
    width = np.ceil(nt_width*num_seq).astype(int)

    if alphabet == 'pu':
        max_height = height
    else:
        max_height = height*2
    #total_height = np.sum(heights,axis=0) # np.minimum(np.sum(heights,axis=0), max_height)
    logo = np.ones((max_height, width, 3)).astype(int)*255;
    for i in range(num_seq):
        nt_height = np.sort(heights[:,i]);
        index = np.argsort(heights[:,i])
        remaining_height = np.sum(heights[:,i]);
        offset = max_height-remaining_height

        for j in range(num_nt):
            if nt_height[j] > 0:
                # resized dimensions of image
                nt_img = np.array(Image.fromarray(chars[index[j]]).resize((nt_width, nt_height[j])))
                # nt_img = imresize(chars[index[j]], (nt_height[j], nt_width))
                # determine location of image
                height_range = range(remaining_height-nt_height[j], remaining_height)
                width_range = range(i*nt_width, i*nt_width+nt_width)

                # 'annoying' way to broadcast resized nucleotide image
                if height_range:
                    for k in range(3):
                        for m in range(len(width_range)):
                            logo[height_range+offset, width_range[m],k] = nt_img[:,m,k];

                remaining_height -= nt_height[j]

    return logo.astype(np.uint8)

