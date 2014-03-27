import os, sys
import numpy as N
import matplotlib.pylab as plt
import get_corr_allscals as corscal
import combination as comb
import scipy.stats.stats as st

"""
to find groups
"""

pccoefArr, allArr, cc, cormat = corscal.do_pca()

nsignif = 4
classify_id = N.zeros(len(cc))
limit = [-0.4, 0.4]

## consider the sign of the correlation coefficients
for i in xrange(nsignif):
	ipos = N.where(cormat[i, :] >= limit[1])[0]
	ineg = N.where(cormat[i, :] <= limit[0])[0]
	print "PC %d" % i
	print ">={:f}:{:}".format(limit[1], ipos)
	print "<={:f}:{:}".format(limit[0], ineg)
	izero = N.where(N.abs(cormat[i, :]) < limit[1])[0]
	print "limit[0]< {:} < limit[1]".format(izero)
	classify_id[ipos] += 3**i*2
	classify_id[ineg] += 3**i
	classify_id[izero] += 0

print classify_id

range = N.arange(min(classify_id), max(classify_id)+1)
ngrp = 0
nmaxele = 0
for ii in range:
	id = N.where(classify_id == ii)[0]
	print '%s%d%s' % ("-"*5, ii, "-"*5)
	for i in id:
		print "{:d}:{:s}: {:}".format(i, cc[i], cormat[:nsignif, i])
	if len(id) != 0:
		ngrp += 1
		if nmaxele < len(id):
			nmaxele = len(id)

classify_id_abs = N.zeros(len(cc))
## taking absolute correlation coefficient
for i in xrange(nsignif):
	ibig = N.where(N.abs(cormat[i, :]) >= limit[1])[0]
	print "PC %d" % i
	print "absolute value >={:.2f}: {:}".format(limit[1], ibig)
	izero = N.where(N.abs(cormat[i, :]) < limit[1])[0]
	print "absolute value < {:.2f}: {:}".format(limit[1], izero)
	classify_id_abs[ibig] += 2**i
	classify_id_abs[izero] += 0

print classify_id_abs

range = N.arange(min(classify_id_abs), max(classify_id_abs)+1)
ngrp = 0
nmaxele = 0
for ii in range:
	id = N.where(classify_id_abs == ii)[0]
	print '%s%d%s' % ("-"*5, ii, "-"*5)
	for i in id:
		print "{:d}:{:s}: {:}".format(i, cc[i], cormat[:nsignif, i])
	if len(id) != 0:
		ngrp += 1
		if nmaxele < len(id):
			nmaxele = len(id)


## for each PC
for i in xrange(nsignif):
	ii = N.where(N.abs(cormat[i, :]) >= limit[1])[0]
	print "%s with PC %d" % ("*"*5, i)
	if len(ii) > 0:
		for jj in ii:
			print "%s, %.3f" % (cc[jj], cormat[i, jj])

"""
## to do PCA of the array which pertains the props from each group
n_grp_element = N.zeros(ngrp, dtype='int')
grp_element = N.zeros((ngrp, nmaxele), dtype='int')
ngrp = 0
for ii in range:
	id = N.where(classify_id == ii)[0]
	if len(id) != 0:
		print '%s%d%s' % ("-"*5, ii, "-"*5)
		n_grp_element[ngrp] = len(id)
		grp_element[ngrp, :len(id)] = id
		ngrp += 1


comb_list = N.zeros((N.prod(n_grp_element), ngrp), dtype='int')
ii = 0
## following only works when nsignif = 4
for iele0 in grp_element[0, :n_grp_element[0]]:
	for iele1 in grp_element[1, :n_grp_element[1]]:
		for iele2 in grp_element[2, :n_grp_element[2]]:
			for iele3 in grp_element[3, :n_grp_element[3]]:
				for iele4 in grp_element[4, :n_grp_element[4]]:
					for iele5 in grp_element[5, :n_grp_element[5]]:
						for iele6 in grp_element[6, :n_grp_element[6]]:
							for iele7 in grp_element[7, :n_grp_element[7]]:
								for iele8 in grp_element[8, :n_grp_element[8]]:
									for iele9 in grp_element[9, :n_grp_element[9]]:
										comb_list[ii, :] = [iele0, iele1, iele2, iele3, iele4, iele5, iele6, iele7, iele8, iele9]
										ii += 1	
ccArr = N.array(cc)
eval0_max = 0.
pccorrfnname = N.corrcoef
for i in xrange(N.prod(n_grp_element)):
	evals, evecs = corscal.get_pcs(allArr[:, comb_list[i, :]], ccArr[comb_list[i, :]], pccorrfnname)
	evals_norm = evals/N.sum(evals)
	if evals_norm[0] > eval0_max:
		print "eval0: %f" % evals_norm[0]
		print "{:}".format(ccArr[comb_list[i, :]])
		eval0_max = evals_norm[0]
		evals_max = evals
		evecs_max = evecs
		selected = allArr[:, comb_list[i, :]]
		selected_cc = ccArr[comb_list[i, :]]

nlos = 96
corcoeffnname = st.pearsonr
pccoefArr = corscal.project_on_pc(evecs_max, selected, selected_cc)
cormat = corscal.calc_corr(pccoefArr, selected, corcoeffnname, 0)
figname = os.path.join("plot", "selectedprops%d_pcvalFromfile%d_%scorrwithpcFrom%s" % (ngrp, 16, "cor", "pea"))
corscal.plot_corr(cormat, selected_cc, evals_max/N.sum(evals_max), figname)
"""
