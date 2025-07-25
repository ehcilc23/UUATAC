'''Motif analysis in NvTK.

Currently, this module only support DNA MOTIF analysis.
'''

import logging
import numpy as np


def trim_ic(motif, cutoff=0.4, background=0.25):
    """Trim motif based on IC(Bernouli)"""
    H = (motif * np.log2(motif / background + 1e-6)).sum(0)
    where = np.where(H > cutoff)[0]
    motif = motif[:, where.min():where.max()+1]
    return motif


def calc_motif_IC(motif, background=0.25):
    """Motif IC Bernouli"""
    H = (motif * np.log2(motif / background + 1e-6)).sum()
    logging.info("Motif IC(Bernouli): %.4f" % H)
    return H


def info_content(pwm, bg=0.5):
    """Motif IC Bernouli"""
    pseudoc = 1e-6
    bg_pwm = [1-bg, bg, bg, 1-bg]
    
    ic = 0 
    for i in range(pwm.shape[0]):
        for j in range(4):
            ic += -bg_pwm[j]*np.log2(bg_pwm[j]) + pwm[i][j]*np.log2(pseudoc + pwm[i][j])
    return ic


def calc_motif_entropy(motif, background=0.25):
    '''Motif Entropy'''
    H = -(motif * np.log2(motif / background + 1e-6)).sum()
    logging.info("Motif Entropy: %.4f" % H)
    return H


def calc_motif_frequency(motif_IC):
    '''Motif Frequency'''
    f = np.power(2, -(motif_IC - 1))
    logging.info("Motif Frequency: %.4f" % f)
    return f


def calc_frequency_W(W, background=0.25):
    '''Calculate motif Frequency in pwms'''
    motif_frequency_l, motif_IC_l = [], []
    for pwm in W:
        pwm = normalize_pwm(pwm)
        motif_IC = calc_motif_IC(pwm)
        motif_freq = calc_motif_frequency(motif_IC)
        motif_IC_l.append(motif_IC); motif_frequency_l.append(motif_freq)
    return motif_frequency_l, motif_IC_l


def normalize_pwm(pwm, factor=None, max=None):
    '''Normalize pwm'''
    if not max:
        max = np.max(np.abs(pwm))
    pwm = pwm/max
    if factor:
        pwm = np.exp(pwm*factor)
    norm = np.outer(np.ones(pwm.shape[0]), np.sum(np.abs(pwm), axis=0))
    pwm = pwm/norm
    pwm[np.isnan(pwm)] = 0.2 # background
    return pwm

threshold_nrom = 0.9

def meme_generate(W, output_file='meme.txt', prefix='Motif_'):
    '''Generate meme file for pwms'''
    # background frequency
    nt_freqs = [1./4 for i in range(4)]

    # open file for writing
    f = open(output_file, 'w')

    # print intro material
    f.write('MEME version 4\n')
    f.write('\n')
    f.write('ALPHABET= ACGT\n')
    f.write('\n')
    f.write('strands: + -\n')
    f.write('\n')
    f.write('Background letter frequencies:\n')
    f.write('A %.4f C %.4f G %.4f T %.4f \n' % tuple(nt_freqs))
    f.write('\n')

    for j in range(len(W)):
        pwm = normalize_pwm(W[j])
        f.write('MOTIF %s%d %d\n' % (prefix, j, j))
        f.write('\n')
        f.write('letter-probability matrix: alength= 4 w= %d nsites= %d E= 0\n' % (pwm.shape[1], pwm.shape[1]))
        for i in range(pwm.shape[1]):
            f.write('  %.4f\t  %.4f\t  %.4f\t  %.4f\t\n' % tuple(pwm[:,i]))
        f.write('\n')

    f.close()

