"""
to compare the difference between mine and Joanne's
"""

import os, sys
import numpy as N
import scipy.stats.stats as st
import get_corr_allscals as corscal
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker

resabspath = os.path.abspath("/Users/yookyung/Research")
filpath = os.path.join(resabspath, "filaments", "observable")

sys.path.append(os.path.join(filpath, "observedmass"))
import observedmass_notsave as obsmass

sys.path.append(os.path.join(resabspath, "mass_anaylsis"))

import reading_fn as rfn

sys.path.append(os.path.join(resabspath, "mass_anaylsis", "corr_scatter_props"))
import read_clus_props as props

def readdata1():
	nodeclus = N.loadtxt(os.path.join(resabspath, "mass_analysis", "outlier_analysis", "nodenum_ffinder_fiducial.dat"))
	esum_yn = N.loadtxt("esum.txt", usecols=[1])
	esum_jc = N.loadtxt("prodtest.dat", usecols=[1])
	observed = N.loadtxt(os.path.join(filpath, "observed_mass_MGe1.00e+14_fullclus.dat"))
	return nodeclus, esum_yn, esum_jc, observed

def esum_compare():
	nlos = 96
	nodeclus, esum_yn, esum_jc, observed = readdata1()
	idiffArr = N.where(N.abs(esum_yn - esum_jc) > 0.01)[0] ## nodenum

	observed1 = N.zeros((len(idiffArr)*nlos, observed.shape[1]))

	for i, idiff in enumerate(idiffArr):
		observed1[i*nlos:(i+1)*nlos,:] = observed[nodeclus[idiff]*nlos:(nodeclus[idiff]+1)*nlos, :]

	for i in xrange(len(idiffArr)):
		print "node %d" % idiffArr[i]
		arr = observed1[i*nlos:(i+1)*nlos,:]
		ii = N.where(arr[:, 5] == 0)[0]
		arr1 = arr[ii, [6, 7, 8, 10, 11]]
		jj = N.where(N.all(arr1 > 0., axis=1))[0]
		arr1 = arr1[jj, :]
		covmat = N.cov(arr1.T)
		eval, evec = N.linalg.eig(arr1)
		print eval, evec



def subtract_average(props):
	for icol in xrange(props.shape[1]):
		props[:, icol] -= N.average(props[:, icol])
	return props


def compare_props_10(nodesonly=True):
	props_jc_10 = N.loadtxt("props_meansubs_10_jc.dat")
	rfilename, cc, props_yn, plt_dir = corscal.get_all(16, nodesonly=nodesonly)
	print rfilename
	print N.array(cc, dtype='str').T

	props_yn_sub = subtract_average(props_yn)
	props_yn_10 = props_yn_sub[:10, :]

	rat = props_yn_10/props_jc_10

	for icol in xrange(rat.shape[1]):
		print "%s" % cc[icol]
		ii = N.argsort(rat[:, icol])	
		print N.column_stack((ii, rat[ii, icol]))
	
	print "%s" % ("="*10)

	for icol in xrange(rat.shape[1]):
		print "%s" % cc[icol]
		ii = N.where(N.abs(rat[:, icol] - 1.) > 0.01)[0]
		print N.column_stack((ii, rat[ii, icol]))

	return props_yn, props_yn_10, props_jc_10
	

def compare_mass():
	observed = rfn.read_data()

	scatArr_yn, cc = props.get_scat_size_ave(5, 96)
	scatArr_jc = N.loadtxt(os.path.join(resabspath, "mass_analysis", "joanne", "deltam.dat"))

	fig = plt.figure(1)
	plt.clf()
	for icol in xrange(scatArr_yn.shape[1]):
		ax = fig.add_subplot(2, 3, icol)
		ax.plot(scatArr_yn[:, icol], scatArr_jc[:, 1+icol], 'ro')
		ax.set_title(cc[icol])

	limit = 0.03
	diff = scatArr_yn - scatArr_jc[:, 1:]
	ii = N.where(N.any(N.abs(diff) > limit, axis=1))[0]
	arrAll = obsmass.get_mass()

	if len(ii) == 0:
		print "no difference in %f" % limit
	else:
		for i in ii:
			print "%s%d%s" % ("-"*5, i, "-"*5)
			print "%d %.4f %.4f %.4f %.4f %.4f\n" % (i, diff[i, 0], diff[i, 1], diff[i, 2], diff[i, 3], diff[i, 4])
			print "%d %.4e %.4e %.4e %.4e %.4e\n" % (i, scatArr_yn[i, 0], scatArr_yn[i, 1], scatArr_yn[i, 2], scatArr_yn[i, 3], scatArr_yn[i, 4])
			print "%d %.4e %.4e %.4e %.4e %.4e\n" % (i, scatArr_jc[i, 1], scatArr_jc[i, 2], scatArr_jc[i, 3], scatArr_jc[i, 4], scatArr_jc[i, 5])
			for ilos in xrange(96):
				print arrAll[i*96+ilos, :]


def compare_allprops(props1, props2, cc):
	fig1 = plt.figure(1)
	plt.clf()
	fig2 = plt.figure(2)
	plt.clf()
	diff = props1 - props2

	for i in xrange(props1.shape[1]):
		print "%s%s%s" % ("+"*5, cc[i], "+"*5)
		ax1 = fig1.add_subplot(4, 5, i+1)
		xx = props1[:, i]
		yy = props2[:, i]
		ax1.plot(xx, yy, 'ro', ms=2)
		ax1.xaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax1.yaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax1.text(N.min(xx), N.min(yy)+0.9*(N.max(yy)-N.min(yy)), cc[i])

		ax2 = fig2.add_subplot(4, 5, i+1)
		ax2.plot(N.arange(len(xx)), diff[:, i], 'ro', ms=2)
		ax2.xaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax2.yaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax2.set_ylim(-0.01, 0.01)
		idiff = N.where(N.abs(diff[:, i]) > 0.01)[0]
		print N.column_stack((idiff, xx[idiff], yy[idiff]))


def plot_ratio_cormats():
	nprops = 19
	cormat_yn = N.loadtxt("cormat_to_do_pca.dat")
	cormat_jc = N.loadtxt("corrcoeff_jc.dat", usecols=[2])
	cormat_jc = cormat_jc.reshape(nprops, nprops)
	rat = cormat_yn/cormat_jc

	fig = plt.figure(1)
	plt.clf()
	plt.subplots_adjust(hspace=0.3, wspace=0.3)
	fontsize=8
	for ip in xrange(nprops): 
		ax = fig.add_subplot(4, 5, ip+1)
		ax.plot(N.arange(nprops), rat[ip, :], 'ro', ms=2.5)
		ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)

		ii = N.where(N.abs(rat[ip, :]-1.) > 0.01)[0]
		print "-"*10
		print "{:d}: {:}".format(ip, ii)
		print cormat_yn[ip, ii]
		print cormat_jc[ip, ii]

def sort_by_evals(evals, evecs):
	ievals = N.argsort(-evals)
	evals = evals[ievals]
	evecs = evecs[:, ievals]
	return evals, evecs

def get_evecs_evals_JCprops(allprops, cormat, cc):
	evals1, evecs1 = corscal.get_pcs(allprops, cc, N.corrcoef)
	evals2, evecs2 = N.linalg.eig(cormat)
	evals2, evecs2 = sort_by_evals(evals2, evecs2)
	return evals1, evecs1, evals2, evecs2

def do_steps_to_do_pca(evecs, allprops, cc, figname):
	nlos = 96
	pccoefArr = corscal.project_on_pc(evecs, allprops, cc, nlos)
	cormat = corscal.calc_corr(pccoefArr, allprops, st.pearsonr, 0)
	corscal.plot_corr(cormat, cc, figname)
	return pccoefArr, cormat


def do_steps_to_do_pca2(evecs, allprops_org, cc, figname):
	nlos = 96
	nprops = len(cc)
	allprops = allprops_org.copy()
	pccoefArr = N.empty_like(allprops)
	for iprops in xrange(allprops.shape[1]):
		allprops[:, iprops] -= N.average(allprops[:, iprops])
	for iprops in xrange(allprops.shape[1]):
		pccoefArr[:, iprops] = N.sum(evecs[:, iprops]*allprops, axis=1)
	cormat = N.zeros((nprops, nprops))
	for iprops in xrange(nprops):
		for jprops in xrange(nprops):
			cormat[iprops, jprops] = N.corrcoef(pccoefArr[:, iprops], allprops[:, jprops])[0, 1]
	corscal.plot_corr(cormat, cc, figname)
	return pccoefArr, cormat

	
pccoefArr_all, props_all_yn, corrcoefmat_all_yn = corscal.do_pca(nodesonly=False)
pccoefArr_nodes, props_nodes_yn, corrcoefmat_nodes_yn = corscal.do_pca(nodesonly=True)
props_all_jc = N.loadtxt(os.path.join(resabspath, "mass_analysis", "joanne", "cluster_stats.dat"))
inode = N.loadtxt(os.path.join(resabspath, "mass_analysis", "outlier_analysis", "nodenum_ffinder_fiducial.dat"), dtype='int')
props_nodes_jc = props_all_jc[inode, :]	
rfilename, cc, props, plot_dir = corscal.get_all(16, nodesonly=False) ## just to get cc

#compare_allprops(props_all_yn, props_all_jc, cc)
#compare_allprops(props_nodes_yn, props_nodes_jc, cc)
#plot_ratio_cormats()


#### to proceed the steps using joanne's properties or correlation matrix

corcoeff_jc = N.loadtxt("corrcoeff_jc.dat", usecols=[2])
cormat_forPC_jc = corcoeff_jc.reshape(19, 19)
evals1_jc, evecs1_jc, evals2_jc, evecs2_jc = get_evecs_evals_JCprops(props_nodes_jc, cormat_forPC_jc, cc)
pccoefArr1_jc, corrcoefmat1_jc = do_steps_to_do_pca(evecs1_jc, props_nodes_jc, cc, "jc_fromprops.eps")
pccoefArr2_jc, corrcoefmat2_jc = do_steps_to_do_pca(evecs2_jc, props_nodes_jc, cc, "jc_fromcormat.eps")

evals1_yn, evecs1_yn = corscal.get_pcs(props_nodes_yn, cc, N.corrcoef)
pccoefArr1_yn, corrcoefmat1_yn = do_steps_to_do_pca2(evecs1_yn, props_nodes_yn, cc, "yn_testbynotshifting.eps")
