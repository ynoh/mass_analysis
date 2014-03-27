"""
2011-10
1. read scalar properties using props
2. get the pc's all scalar properties
3. calcluate the correlation between scalar props and their projected values on pcs

2012-02
changing the eigenvector sign to fix projected delta M be positive
"""

import os, sys
import numpy as N
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
sys.path.append("..")
sys.path.append(os.path.join("..", ".."))
import reading_fn as rfn
import read_clus_props as props
import scipy.stats.stats as st
import combination


def change_pc_sign(arr):
	"""
	delta M_red be positive
	each column is corresponding to the eigenvector
	"""
	for i in xrange(arr.shape[1]):
		if arr[0, i] < 0.:
			arr[:, i] *= -1.

	return arr
	

def get_pcs(arr, arrname, corrfn):
	"""
	corrfn = N.corrcoef or N.cov
	"""
	#print "array shape to get pcs={:}".format(arr.shape)
	#arr = arr.T
	cormat = corrfn(arr.T)
	evals, evecs = N.linalg.eig(cormat)
	#N.savetxt("cormat_to_do_pca.dat", cormat)
	ievals = N.argsort(-evals)
	evals = evals[ievals]
	evecs = evecs[:, ievals]

	#print "{:s} cormat: {:}".format("-"*3, cormat)
	#print "{:s} evals: {:}".format("-"*3, evals)
	#print "{:s} normalized evals:{:}".format("-"*3, evals/N.sum(evals))
	#print "{:s} 1st, 2nd, 3rd pcs of each props:".format("-"*3)
	#for i in xrange(len(arrname)):
	#	print "{:s}: {:.3f} {:.3f} {:.3f}".format(arrname[i], evecs[i, 0], evecs[i, 1], evecs[i, 2]) 
	return evals, evecs

def project_on_pc(evecs, mprops_org, cc):
	"""
	evecs = n_mprops x n_mprops
	mprops_org = nclus x n_mprops : before subtracting the average
	"""
	mprops = mprops_org.copy()
	n_mprops = mprops.shape[1]
	print "nprops:%d in project on pc" % n_mprops
	if evecs.shape[1] != n_mprops:
		print "pc and props are not consistent"
		raise

	print "{:s}:{:}".format("="*5, mprops[0, :])
	## subtract the average and divided by std
	for iprops in xrange(n_mprops):
		print "process: %s" % cc[iprops]
		print N.average(mprops[:, iprops])
		mprops[:, iprops] = (mprops[:, iprops] - N.average(mprops[:, iprops]))/N.std(mprops[:, iprops], ddof=1)
		print N.min(mprops[:, iprops]), N.max(mprops[:, iprops])

	print "{:s}:{:}".format("="*5, mprops[0, :])
	pccoefArr = N.empty_like(mprops)
	for iprops in xrange(n_mprops):
		pccoefArr[:, iprops] = N.sum(evecs[:, iprops]*mprops, axis=1)

	return pccoefArr


def calc_corr(pccoefArr, mprops, corrfn, iget_cor):
	n_mprops = mprops.shape[1]
	cormat = N.empty((n_mprops, n_mprops))
	for i_mprops in xrange(n_mprops):
		for j_mprops in xrange(n_mprops):
			corrcoef = corrfn(pccoefArr[:, i_mprops], mprops[:, j_mprops])[iget_cor]
			cormat[i_mprops, j_mprops] = corrcoef
	print "%s correlation coeff with PC0 coeff. and props" % ("-"*3)
	print cormat[0, :]
	return cormat


def modified_cov(arr1, arr2):
	covmat = N.cov(arr1, arr2)
	covmat = covmat.reshape(4)
	return covmat

def plot_corr2(cormat, cc, evals, figname, xlabel="props"):
	"""
	cormat = (pccoef) x (mprops)
	change the plotting shape
	"""
	fig = plt.figure(2)
	fig.clf()
	n_mprops = cormat.shape[1]

	#mec_arr = ['blue', 'red', 'orange', 'cyan', 'green', 'white', 'white', 'white', 'white', 'white']
	#mew_arr = ['blue', 'red', 'orange', 'cyan', 'green', 'blue', 'red', 'orange', 'cyan', 'green'] 
	if xlabel == "props":
		nplot = 4
		nrow = 2
		ncol = 2
		for i_mprops in xrange(nplot):
			ax = fig.add_subplot(nrow, ncol, i_mprops+1)
			ax.plot(N.arange(n_mprops), cormat[i_mprops, :], 'ro')
			ax.axhline(y=-0.4, xmin=-.5, xmax=n_mprops, linestyle=':')
			ax.axhline(y=0.4, xmin=-.5, xmax=n_mprops, linestyle=':')
			ax.set_ylim(-1., 1.)
			ax.set_xlim(-.5, n_mprops-.5)
			plt.subplots_adjust(hspace=0., bottom=0.2, wspace=0.2)
			if i_mprops >= (nrow-1)*ncol:
				ax.set_xticklabels(cc, rotation='vertical', fontsize=12)
				ax.set_xticks(N.arange(n_mprops))
			else:
				ax.set_xticklabels("")
				ax.set_xticks(N.arange(n_mprops))
			#if i_mprops % nrow == 0:
			ax.set_ylabel(r"$Cor$ w/ PC%d(%.2f)" % (i_mprops, evals[i_mprops]))
			if len(figname) > 0:
				plt.savefig("%s_xaxis%s.eps" % (figname, xlabel))

	else:
		mfc_arr = ['b', 'g', 'r', 'c', 'm', 'w', 'w', 'w', 'w', 'w']
		mec_arr = ['b', 'g', 'r', 'c', 'm', 'b', 'g', 'r', 'c', 'm'] 
		ncolor = len(mfc_arr)
		for i_mprops in xrange(n_mprops):
			if i_mprops < ncolor:
				plt.plot(N.arange(n_mprops), cormat[:, i_mprops], mfc=mfc_arr[i_mprops], mec=mec_arr[i_mprops], marker = 'o', label=cc[i_mprops], linestyle="None")
			else:
				plt.plot(N.arange(n_mprops), cormat[:, i_mprops - ncolor], mfc=mfc_arr[i_mprops - ncolor], mec=mec_arr[i_mprops - ncolor], marker = 's', label=cc[i_mprops - ncolor], linestyle="None")
			plt.ylim(-1., 1.)
			plt.xlabel("PC")
			plt.ylabel(r"$Cor$")
			plt.legend(loc=4, numpoints=1)
			plt.xlim(-0.5, n_mprops+2)
			if len(figname) > 0:
				plt.savefig("%s_xaxis%s.eps" % (figname, xlabel))

def plot_corr(cormat, cc, evals, figname, xlabel="props"):
	"""
	cormat = (pccoef) x (mprops)
	"""
	fig = plt.figure(2)
	fig.clf()
	n_mprops = cormat.shape[1]
	fs = dict(lb1 = 16, lb2 = 17, tc = 13)
	dy = 0.6
	dyy = 0.2
	plt.rcParams["axes.unicode_minus"] = False

	#mec_arr = ['blue', 'red', 'orange', 'cyan', 'green', 'white', 'white', 'white', 'white', 'white']
	#mew_arr = ['blue', 'red', 'orange', 'cyan', 'green', 'blue', 'red', 'orange', 'cyan', 'green'] 
	if xlabel == "props":
		nplot = 4
		for i_mprops in xrange(nplot):
			ax = fig.add_subplot(nplot, 1, i_mprops+1)
			plt.subplots_adjust(hspace=0, bottom=0.22, left=0.08, right=0.9)
			ax.plot(N.arange(n_mprops), cormat[i_mprops, :], 'ro', ms=3.8)
			ax.axhline(y=-0.4, xmin=-.5, xmax=n_mprops, linestyle=':')
			ax.axhline(y=0.4, xmin=-.5, xmax=n_mprops, linestyle=':')
			ax.set_ylim(-1., 1.)
			ax.set_xlim(-.5, n_mprops-.5)
			ax.set_xticks(N.arange(n_mprops))
			if i_mprops == nplot-1:
				ax.set_xticklabels(cc, rotation=90, fontsize=fs['lb1'], position=[0.5, 0.0])
				#verticalalignment='top', horizontalalignment='right')
				#ax.set_xticklabels(cc, rotation='vertical', fontsize=fs['lb1'])
			else:
				plt.setp(ax.get_xticklabels(), visible=False)
			#ax.set_ylabel(r"$\xi$ w PC%d(%.2f)" % (i_mprops, evals[i_mprops]), fontsize=10)
			ax.set_ylabel(r"$PC_{%d}(%.2f)$" % (i_mprops, evals[i_mprops]), fontsize=fs['lb1'])
			ax.yaxis.set_label_coords(1.05, 0.5)
			ax.yaxis.set_major_locator(plt_ticker.MultipleLocator(dy))
			ax.yaxis.set_minor_locator(plt_ticker.MultipleLocator(dyy))
			plt.setp(ax.get_yticklabels(), fontsize=fs['tc'])
			if i_mprops == (nplot)/2:
				ax2 = fig.add_subplot(nplot, 1, i_mprops+1, sharex=ax, sharey=ax, axisbg="none")
				ax2.set_ylabel(r"$Cor$", fontsize=fs['lb2'])
				ax2.yaxis.set_label_coords(-0.06, 0.95)
				ax2.yaxis.set_major_locator(plt_ticker.MultipleLocator(dy))
				ax2.yaxis.set_minor_locator(plt_ticker.MultipleLocator(dyy))
				plt.setp(ax2.get_xticklabels(), visible=False)
				plt.setp(ax2.get_yticklabels(), visible=False)
			if len(figname) > 0:
				plt.savefig("%s_xaxis%s.eps" % (figname, xlabel))

	else:
		mfc_arr = ['b', 'g', 'r', 'c', 'm', 'w', 'w', 'w', 'w', 'w']
		mec_arr = ['b', 'g', 'r', 'c', 'm', 'b', 'g', 'r', 'c', 'm'] 
		ncolor = len(mfc_arr)
		for i_mprops in xrange(n_mprops):
			if i_mprops < ncolor:
				plt.plot(N.arange(n_mprops), cormat[:, i_mprops], mfc=mfc_arr[i_mprops], mec=mec_arr[i_mprops], marker = 'o', label=cc[i_mprops], linestyle="None")
			else:
				plt.plot(N.arange(n_mprops), cormat[:, i_mprops - ncolor], mfc=mfc_arr[i_mprops - ncolor], mec=mec_arr[i_mprops - ncolor], marker = 's', label=cc[i_mprops - ncolor], linestyle="None")
			plt.ylim(-1., 1.)
			plt.xlabel("PC")
			plt.ylabel(r"$Cor$")
			plt.legend(loc=4, numpoints=1)
			plt.xlim(-0.5, n_mprops+2)
			if len(figname) > 0:
				plt.savefig("%s_xaxis%s.eps" % (figname, xlabel))


def plot_org_scatter(mprops, cc, figname):
	fig = plt.figure(3)

	plt.clf()
	n_mprops = mprops.shape[1]	
	fontsize=6
	labelx = -0.2
	param = dict(fontsize='xx-small')
	plotadjust = dict(hspace=0., wspace=0., left=0.1, right=0.95, bottom=0.15, top=0.95)
	width = (plotadjust['right'] - plotadjust['left'])/(n_mprops*(1. + plotadjust['wspace']))
	height = (plotadjust['top'] - plotadjust['bottom'])/(n_mprops*(1. + plotadjust['hspace']))


	for j_mprops in xrange(n_mprops):
		#for j_mprops in xrange(0, n_mprops-i_mprops):
		for i_mprops in xrange(j_mprops, n_mprops):
			pos = [plotadjust['right'] - i_mprops*(1. + plotadjust['wspace'])*width, plotadjust['bottom'] + j_mprops*(1. + plotadjust['hspace'])*height, width, height]
			#pos = [plotadjust['left'] + i_mprops*(1. + plotadjust['wspace'])*width, plotadjust['bottom'] + j_mprops*(1. + plotadjust['hspace'])*height, width, height]
			#ax = fig.add_subplot(n_mprops, n_mprops, i_mprops*n_mprops + j_mprops + 1)
			#plt.subplots_adjust(hspace=0.2, wspace=0.2, left=0.08, right=0.95, bottom=0.08, top=0.9)
			ax = fig.add_axes(pos)
			ax.plot(mprops[:, i_mprops], mprops[:, j_mprops], 'ro', ms=1.5)
			xr = [N.min(mprops[:, i_mprops]), N.max(mprops[:, i_mprops])]
			yr = [N.min(mprops[:, j_mprops]), N.max(mprops[:, j_mprops])]
			sx = []
			if j_mprops == 0:
				ax.set_xlabel(cc[i_mprops], **param)
				#ax.set_xticks(xr)
				#ax.set_xticklabels(["%.2f" % xr[0], "%.2f" % xr[1]])
				ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(3))
				for xx in N.linspace(xr[0], xr[1], 2):
					sx.append('{:.2f}'.format(xx))
				ax.set_xticklabels(sx)
				for tick in ax.xaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
			else:
				ax.set_xticklabels([])
			sy = []
			if i_mprops == n_mprops-1:
				ax.set_ylabel(cc[j_mprops], **param)
				for yy in N.linspace(yr[0], yr[1], 2):
					sy.append('{:.2f}'.format(yy))
				#print sy
				ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(3))
				ax.set_yticklabels(sy)
				for tick in ax.yaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				#print cc[j_mprops]
			else:
				ax.set_yticklabels([])
				#ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(2))
				#ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(2))
				#ax.set_xticks(["", ""])
				#ax.set_yticks(["", ""])


	plt.savefig(figname)



def plot_proj_scatter(pccoefs, eval, figname):
	"""
	plot scatters using the values projected on pc
	"""
	fig = plt.figure(4)
	plt.clf()

	eval /= N.sum(eval)

	param = dict(fontsize='small')
	
	range = N.max(N.abs(pccoefs))
	labelx = -0.15

	nplt = 4
	ncol = 3
	nrow = 2
	iplt = 1
	for iprops in xrange(nplt):
		for jprops in xrange(iprops+1, nplt):
			ax = fig.add_subplot(nrow, ncol, iplt)
			plt.subplots_adjust(hspace=0.2, wspace=0.25, left=0.08, right=0.95, bottom=0.08, top=0.9)
			ax.plot(pccoefs[:, iprops], pccoefs[:, jprops], 'ro', ms=3.)
			plt.axis('equal')
			#ax.set_xlim(-range, range)
			#ax.set_ylim(-range, range)
			ax.set_xlabel('PC%d(%.2f)' % (iprops, eval[iprops]), **param)
			ax.set_ylabel('PC%d(%.2f)' % (jprops, eval[jprops]), **param)
			ax.yaxis.set_label_coords(labelx, 0.5)
			iplt += 1

	plt.savefig(figname)



def plot_orgvsproj_scatter(ipc, pccoefs, mprops, corrcoefs, eval, cc, figname):
	"""
	plot mprops vs one of projected mprops
	"""
	fig = plt.figure(4)
	plt.clf()

	eval /= N.sum(eval)
	nprops = len(cc)

	fig = plt.figure(1)
	plt.clf()
	param = dict(fontsize='xx-small')
	
	labelx = -0.12

	ncol = 3
	nrow = N.ceil(1.*nprops/ncol) 
	iplt = 1
	fontsize = 7
	for iprops in xrange(nprops):
		ax = fig.add_subplot(nrow, ncol, iplt)
		plt.subplots_adjust(hspace=0., wspace=0.2, left=0.08, right=0.95, bottom=0.08, top=0.9)
		ax.plot(pccoefs[:, ipc], mprops[:, iprops], 'ro', ms=3.)
		#plt.axis('equal')
		#ax.set_xlim(-range, range)
		#ax.set_ylim(-range, range)
		ax.set_xlabel('PC%d(%.2f)' % (ipc, eval[ipc]), **param)
		ax.set_ylabel('%s' % cc[iprops], **param)
		ax.yaxis.set_label_coords(labelx, 0.5)
		ax.text(0.8*N.max(pccoefs[:, ipc]), 0.8*N.max(mprops[:, iprops]), "%.3f" % corrcoefs[ipc, iprops], **param)
		iplt += 1
		ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(4))
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)

	plt.savefig(figname)


def find_each_prop(prop_to_pick, propArr, ccList):
	ii = []
	print "choose %s" % prop_to_pick
	for i, cc in enumerate(ccList):
		if cc.find(prop_to_pick) != -1:
			ii.append(i)
	return N.array(ii, dtype='int')
		


def plot_scatArr(scatArr, name, figname, take_abs_scat):
	"""
	plot histogram of (<M_obs> - M_true)/ M_true
	"""
	fig = plt.figure(10)
	plt.clf()
	maxval = N.max((N.abs(N.min(scatArr)), N.max(scatArr)))
	if take_abs_scat:
		histbin = N.linspace(0, maxval, 21)
	else:
		histbin = N.linspace(-maxval, maxval, 21)
	nmax = 100.
	for iobs in xrange(len(name)):
		ax = fig.add_subplot(3, 2, iobs+1)
		plt.subplots_adjust(hspace=0.4)
		h, bins, patches = plt.hist(scatArr[:, iobs], bins=histbin, histtype='step')
		ax.set_xlabel(name[iobs])
		if iobs % 2 == 0:
			ax.set_ylabel(r"$N_{cluster}$")
		if N.max(h) > nmax:
			ax.set_ylim((0., N.max(h)))
		else:
			ax.set_ylim((0., nmax))
	plt.savefig(figname)


def comb(n, k):
	if k > n/2:
		kk = n - k
	else:
		kk = k
	numer = N.prod(N.arange(n, n - kk, -1))
	denom = N.prod(N.arange(kk, 1, -1))
	return numer/denom

def construct_small_array0(allArr, ccList, npick=9):
	"""
	npick counts everything
	"""
	ccArr = N.array(ccList)
	nclus = len(allArr)
	parArr = N.empty((nclus, npick)) ## will contain all the chosen properties

	## sum scatter size
	icols = []
	for ii, cc in enumerate(ccList):
		if cc.find("Delta") != -1:
			icols.append(ii)
	parArr[:, 0] = N.sum(allArr[:, icols], axis=1)
	cc0 = [r"$\sigma_i \Delta M_i$"]
	
	allArr2 = N.delete(allArr, icols, 1)  ## delete all mass scatter
	ccArr2 = N.delete(ccArr, icols)
	ncol = allArr2.shape[1] 
	jcols = combination.Combination(ncol, npick-1, 0, comb(ncol, npick-1)) ## lists of columns to pick
	jcols = N.array(jcols)
	print jcols.shape

	## this routine always keep the sum of the scatter sizes
	evalmax = 0.
	jcolmax = jcols[0, :]
	for jj in jcols:
		jj = N.array(jj, dtype='int') - 1
		parArr[:, 1:] = allArr2[:, jj]
		ccparList = cc0 + list(ccArr2[jj])
		evals, evecs = get_pcs(parArr, ccparList, N.corrcoef)
		normedeval0 = evals[0]/N.sum(evals)
		evalmax = max(normedeval0, evalmax)
		if evalmax == normedeval0:
			jcolmax = N.append(jj+1, 0)
			jcolmax = N.roll(jcolmax, 1)
			parArrmax = parArr.copy()
			ccparListmax = ccparList
			#print "{:}: {:.2f}, {:.2f}".format(jcolmax, normedeval0, evalmax)
	#pccoefArr = project_on_pc(evecs, allArr, cc, nlos)
	#cormat = calc_corr(pccoefArr, allArr, corcoeffn, 0)

	return jcolmax, ccparListmax, parArrmax

def find_props_cols(ccList, str):
	cols = []
	for ii, cc in enumerate(ccList):
		if cc.find(str) != -1:
			cols.append(ii)
	return cols

def construct_small_array(allArr, ccList, npick=9):
	"""
	npick counts everything
	pick only one of PC props
	"""
	ccArr = N.array(ccList)
	nclus = len(allArr)
	parArr = N.empty((nclus, npick)) ## will contain all the chosen properties

	## sum scatter size
	sc_cols = find_props_cols(ccList, "Delta")
	print "sc cols:{:}".format(sc_cols)
	parArr[:, 0] = N.sum(allArr[:, sc_cols], axis=1)
	cc0 = [r"$\sigma_i \Delta M_i$"]

	## find pc related array cols
	pc_cols = find_props_cols(ccList, "PC")
	print "pc cols:{:}".format(pc_cols)
	pcpropArr = allArr[:, pc_cols]
	pcpropccList = list(ccArr[pc_cols])
	#print pcpropccList

	## number of fixed props
	nfix = 2

	## columns to delete
	icols = N.array(sc_cols)
	icols = N.append(icols, pc_cols)
	
	## remaining props	
	others = N.delete(allArr, icols, 1)  ## delete all mass scatter
	ccArr2 = N.delete(ccArr, icols)
	ncol = others.shape[1] 
	jcols = combination.Combination(ncol, npick-nfix, 0, comb(ncol, npick-nfix)) ## lists of columns to pick
	jcols = N.array(jcols)
	print jcols.shape

	## this routine always keep the sum of the scatter sizes
	evalmax = 0.
	jcolmax = jcols[0, :] ## initialize
	for jj in jcols:
		jj = N.array(jj, dtype='int') - 1 ## combination returns numbers starting 1
		parArr[:, nfix:] = others[:, jj]
		for ipcprop in xrange(pcpropArr.shape[1]):
			ccparArr = N.array(cc0)
			parArr[:, 1] = pcpropArr[:, ipcprop]
			ccparArr = N.append(ccparArr, pcpropccList[ipcprop])
			ccparArr = N.append(ccparArr, ccArr2[jj])
			#print ccparList
			evals, evecs = get_pcs(parArr, ccparArr, N.corrcoef)
			normedeval0 = evals[0]/N.sum(evals)
			evalmax = max(normedeval0, evalmax)
			if evalmax == normedeval0:
				jcolmax = N.append(jj+1, 0)
				jcolmax = N.roll(jcolmax, 1)
				parArrmax = parArr.copy()
				ccparListmax = ccparArr
				#print "{:}: {:.2f}, {:.2f}".format(jcolmax, normedeval0, evalmax)
	#pccoefArr = project_on_pc(evecs, allArr, cc, nlos)
	#cormat = calc_corr(pccoefArr, allArr, corcoeffn, 0)

	return jcolmax, ccparListmax, parArrmax

		

def do_pca(ifile=16, nodesonly=True, mass_limit_env=5.e13, pccorrfn=N.corrcoef, corcoeffn=st.pearsonr, take_abs_scat=False, fix_sign=True, plot_shift=False, mtype_like="mo2"):
	"""
	modified 2012/02/08
	fix sign: fix delta M coeff to be positive
	pccorrfn: which corrfn will be used to calculate the PC's
	corcoeffn: which corrfn will be used to calculate the corrln between projected values on PCs and the scalar properties
	plot_shift: plot shif
	mtype_like related prop is taken out (12/04/??)
	"""
	nobs = 5
	nlos = 96
	nbox = 10
	filenameArr = rfn.get_filename(nbox)
	rfilename = filenameArr[ifile]
	ccList, allArr = props.get_array_for_all_scals(ifile, nodesonly=nodesonly, mass_limit_env=mass_limit_env, take_abs_scat=take_abs_scat, mtype_like=mtype_like)
	## allArr contains cluster information and scalar properties
	clus_info = allArr[:, :3]
	allArr = allArr[:, 3:]
	
	if pccorrfn == N.cov:
		pccorrfnname = "cov"
	else:
		pccorrfnname = "cor"

	if corcoeffn == st.pearsonr:
		corcoeffnname = "pea"
	elif corcoeffn == st.spearman:
		corcoeffnname = "spe"
	elif corcoeffn == "modified_cov":
		corcoeffnname = "cov"

	## take absolute value of Delta M_obs?
	if take_abs_scat:
		plot_dir = "plot_abs_test"
	else:
		plot_dir = "plot"

	if plot_shift:
		iprops = find_each_prop("Delta", allArr, ccList)
		plot_scatArr(allArr[:, iprops], ccList[iprops], os.path.join(plot_dir, "shiftsize.eps"), take_abs_scat)

	evals, evecs = get_pcs(allArr, ccList, pccorrfn)
	print "after getting evals and evecs allArr shape{:}".format(allArr.shape)
	if fix_sign:
		evecs = change_pc_sign(evecs)
	pccoefArr = project_on_pc(evecs, allArr, ccList)
	cormat = calc_corr(pccoefArr, allArr, corcoeffn, 0)

	#figname = os.path.join(plot_dir, "allscalprops%d_pcvalFrom%s_%scorrwithpcFrom%s" % (len(ccList), rfilename, corcoeffnname, pccorrfnname))
	## filename is renamed
	figname = os.path.join(plot_dir, "%s_scal%d-pcval-%s_w_pc-%s_node%s_mlike%s" % (corcoeffnname, len(ccList), rfilename, pccorrfnname, nodesonly, mtype_like))
	plot_corr(cormat, ccList, evals, figname)
	#plot_corr(cormat, ccList, evals/N.sum(evals), figname, xlabel="PC")
	#figname = os.path.join(plot_dir, "scatter_scal%d_pcvalFrom%s_node%s.eps" % (len(ccList), rfilename, nodesonly))
	#plot_org_scatter(allArr, ccList, figname)
	#figname = os.path.join(plot_dir, "proj_scat_scal%d-pcval-%s_%s_node%s.eps" % (len(ccList), rfilename, pccorrfnname, nodesonly))
	#plot_proj_scatter(pccoefArr, evals, figname)
	#figname = os.path.join(plot_dir, "scatter_allscalprops_pcvalFrom%s_vs_pc0From%s_coefFrom%s.eps" % (rfilename, pccorrfnname, corcoeffnname))
	#plot_orgvsproj_scatter(0, pccoefArr, allArr, cormat, evals, ccList, figname)
	#figname = os.path.join(plot_dir, "scatter_allscalprops_pcvalFrom%s_vs_pc1From%s_coefFrom%s.eps" % (rfilename, pccorrfnname, corcoeffnname))
	#plot_orgvsproj_scatter(1, pccoefArr, allArr, cormat, evals, ccList, figname)
	#figname = os.path.join(plot_dir, "scatter_allscalprops_pcvalFrom%s_vs_pc2From%s_coefFrom%s.eps" % (rfilename, pccorrfnname, corcoeffnname))
	#plot_orgvsproj_scatter(2, pccoefArr, allArr, cormat, evals, ccList, figname)

	print "%s print for tex %s" % ("-"*3, "-"*3)
	for i in xrange(len(ccList)):
		print "{:s} & {:.2f}({:.2f}) & {:.2f}({:.2f}) & {:.2f}({:.2f}) & {:.2f}({:.2f})\\\\".format(ccList[i], evecs[i, 0], cormat[0, i], evecs[i, 1], cormat[1, i], evecs[i, 2], cormat[2, i], evecs[i, 3], cormat[3, i]) 
	print "%s %s" % ("-"*3, "-"*3)

	individual_mass = False ## doesn't work
	if individual_mass:
		for iobs in xrange(nobs):
			print "%s observable #%d" % ("*"*5, iobs)
			print "check scatter arr: {:}".format(scatArr[:3, iobs])
			parArr = N.column_stack((scatArr[:, iobs], dotprodArr, evalprops))
			evals_i, evecs_i = get_pcs(parArr, [cc_scat[iobs]]+cc_dotprod+cc_eval, pccorrfn)
			pccoefArr_i = project_on_pc(evecs_i, parArr)
			cormat_i = calc_corr(pccoefArr_i, parArr, corcoeffn, 0)
			figname = os.path.join(plot_dir, "mass%dallscalprops_pcvalFrom%s_%scorrwithpcFrom%s" % (iobs, rfilename, corcoeffnname, pccorrfnname))
			plot_corr(cormat_i, [cc_scat[iobs]] + cc_dotprod + cc_eval, figname)

	return pccoefArr, allArr, ccList, cormat

#def do_pca_small(npick=9, nodesonly=True, mass_limit_env=5.e13):
def do_pca_small(cc, allArr, plot_dir, rfilename, npick=9, make_subset=False):
	"""
	if make_subset is True, make a subset in this function.
	otherwise, place subset of cc as cc and subset of the scalar props as allArr.
	"""
	nobs = 5
	#cc, allArr = props.get_array_for_all_scals(16, nodesonly=nodesonly, mass_limit_env=mass_limit_env)
	
	## for everything
	pccorrfn = N.corrcoef
	#pccorrfn = N.cov
	corcoeffn = st.pearsonr

	if pccorrfn == N.cov:
		pccorrfnname = "cov"
	else:
		pccorrfnname = "cor"

	if corcoeffn == st.pearsonr:
		corcoeffnname = "pea"
	elif corcoeffn == st.spearman:
		corcoeffnname = "spe"
	elif corcoeffn == "modified_cov":
		corcoeffnname = "cov"

	if make_subset:
		jcols, ccparList, parArr = construct_small_array(allArr, cc, npick=npick)
	else:
		ccparList = cc
		parArr = allArr	
		npick = len(cc)
	evals, evecs = get_pcs(parArr, ccparList, pccorrfn)
	print "after getting evals and evecs allArr shape{:}".format(parArr.shape)
	pccoefArr = project_on_pc(evecs, parArr, ccparList)
	cormat = calc_corr(pccoefArr, parArr, corcoeffn, 0)

	figname = os.path.join(plot_dir, "parscalprops%d_pcvalFrom%s_%scorrwithpcFrom%s" % (npick, rfilename, corcoeffnname, pccorrfnname))
	plot_corr(cormat, ccparList, evals/N.sum(evals), figname)
	#plot_corr(cormat, cc, figname, xlabel="PC")
	#figname = os.path.join(plot_dir, "scatter_allscalprops_pcvalFrom%s.eps" % rfilename)
	#plot_org_scatter(allArr, cc, figname)
	#figname = os.path.join(plot_dir, "scatter_allscalprops_pcvalFrom%s_proj_onpcFrom%s.eps" % (rfilename, pccorrfnname))
	#plot_proj_scatter(pccoefArr, evals, figname)
	if make_subset:
		return jcol, ccparList, parArr, evals, evecs, pccoefArr, cormat
	else:
		return evals, evecs, pccoefArr, cormat


"""
##some analysis after running do_pca, to put in the paper

## for figure 9
pccoefArr, allArr, ccList, cormat = do_pca()
evals, evecs = get_pcs(allArr, ccList, N.corrcoef)
print ["%.2f" % e for e in evals/N.sum(evals)]
## this cormat is the correlation matrix between pccoefs and the properties
## calculate the correlation matrix between the properties
cormat_props_nn = N.corrcoef(allArr.T)
cormat_props = cormat_props_nn.flatten()
#boundary = [-1., -0.5, -0.35, -0.2, 0., 0.2, 0.35, 0.5, 1.01] ## 1.01 to include 1.00
boundary = [-1., -0.5, -0.35, -0.2, 0.2, 0.35, 0.5, 1.01] ## 1.01 to include 1.00

nr = len(boundary)-1
nprops = allArr.shape[1]
ii = 0
for i in xrange(nr):
	#if boundary[i] < -0.2 or boundary[i+1] > 0.2:
	print "%.2f <= xi %.2f" % (boundary[i], boundary[i+1])
	ii = N.where(N.all((cormat_props >= boundary[i], cormat_props < boundary[i+1]), axis=0))[0]
	ix = ii/nprops
	iy = ii%nprops
	jj = 0
	for ixx, iyy in zip(ix, iy):
		if iyy > ixx:
			#print ixx, iyy
			print "%s(%d)-%s(%d): %.2f" % (ccList[ixx], ixx, ccList[iyy], iyy, cormat_props[ii[jj]])
		jj += 1

## for some discussion
nclus = len(allArr)
nlos = 96
nobs = len(evals)

observed, obsname = rfn.read_data()
maxMobs = N.zeros((nclus, nobs))

for iobs in xrange(nobs):
	aa = observed[:, iobs+2].reshape(nclus, nlos)


cormat_abs = N.abs(cormat)
print "%s max correlation with PC %s" % ("*"*5, "*"*5)
for ipc in xrange(cormat.shape[0]):
	ii = N.where(cormat_abs[ipc, :] == N.max(cormat_abs[ipc, :]))[0]
	print "PC%d: %s (%.3f)" % (ipc, ccList[ii], cormat[ipc, ii])

corbig = 0.5
print "%s correlation of prop with PC >= %.2f %s" % ("*"*5, corbig, "*"*5)
for iprop in xrange(cormat.shape[1]):
	jj = N.where(cormat_abs[:, iprop] >= corbig)[0]
	print "{:s}: {:}: {:}".format(ccList[iprop], jj, cormat[jj, iprop])
	print "the smallest PC: %d" % N.min(jj)

## to put in the do_pca_small()
jj = N.where(cormat_abs[0, :] >= 0.4)[0]
ccparList = N.array(ccList, dtype='str')
ccparList = list(ccparList[jj])
scalprops = allArr[:, jj]
do_pca_small(ccparList, scalprops, "plot", rfn.get_filename(10)[16], npick=9, make_subset=False)
"""
