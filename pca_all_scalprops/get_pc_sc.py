"""
do principal component analysis for the entire clusters using scalar properties
this is basically same as what Jeeson-Daniel et al did

add drawing the histogram of each property
"""

import os
import numpy as N
import scipy.stats.stats as st
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
from matplotlib import cm, colors
from matplotlib.backends.backend_pdf import PdfPages

def get_pcs(arr, arrname):
	arr = arr.T
	cormat = N.corrcoef(arr)
	evals, evecs = N.linalg.eig(cormat)
	ievals = N.argsort(-evals)
	evals = evals[ievals]
	evecs = evecs[:, ievals]

	print "{:s}evals:{:}".format("-"*5, evals)
	print "{:s}normalized evals:{:}".format("-"*5, evals/N.sum(evals))
	print "{:s}1st, 2nd, 3rd pcs of each prop:"
	print ["{:s}: {:3.2f} {:3.2f} {:3.2f}".format(arrname[i], evecs[i, 0], evecs[i, 1], evecs[i, 2]) for i in xrange(len(arrname))]

	return evals, evecs


def get_propsname(fname):
	ff = open(fname, "r")
	ff.readline()
	comments = ff.readline()
	ff.close()
	comments = comments.split()
	comments = N.array(comments)
	print comments
	return comments


def take_logmass(propsname, arr):
	"""
	take log for only masses
	"""
	icol = N.array([], dtype="int")
	icol = N.append(icol, N.where(propsname == "fof-mass")[0])
	icol = N.append(icol, N.where(propsname == "CFP-CFM")[0])
	icol = N.append(icol, N.where(propsname == "AMP-AHM")[0])
	icol = N.append(icol, N.where(propsname == "AFP-AFM")[0])
	icol = N.append(icol, N.where(propsname == "GRP-AGR")[0])

	if len(icol) > 0:
		print "{:} taking log of masses: {:}".format("-"*5, propsname[icol])
		arr[:, icol] = N.log10(arr[:, icol] + 1.e-7)
	return icol, arr


def take_logall(arr):
	"""
	take log for all properties
	"""
	print "%s %s %s" % ("-"*5, "taking log for everything", "-"*5)
	arr = N.log10(arr + 1.e-7)
	return arr


def calc_fluctuation(arr):
	"""
	to use the deviation instead of using the true value
	i.e. subtract the average and divide by the average
	"""
	## arr: (nclusxnprops) 
	arr = arr.T ## to write code simply
	print "%s %s %s" % ("-"*5, "scatter will be used", "-"*5)
	print "array shape:{:}".format(arr.shape)
	for ii, each_row in enumerate(arr):
		ave = N.average(each_row)
		print "%d: %.4f, %.4f" % (ii, ave, N.std(each_row))
		each_row = (each_row - ave)/N.std(each_row)
		arr[ii, :] = each_row

	return arr.T ## re-transpose


def get_node_clus(sprops, propsname):
	"""
	get node clusters: these clusters have planes or plane mass fraction
	"""
	icol = N.array([], dtype='int')
	icol = N.append(icol, N.where(N.any([propsname == 'CFP-CFM', propsname == 'FrCFP-CFM'], axis=0))[0])
	icol = N.append(icol, N.where(N.any([propsname == 'AFP-AFM', propsname == 'FrAFP-AFM'], axis=0))[0])
	icol = N.append(icol, N.where(N.any([propsname == 'AMP-AHM', propsname == 'FrAMP-AHM'], axis=0))[0])
	icol = N.append(icol, N.where(N.any([propsname == 'GRP-AGR', propsname == 'FrGRP-AGR'], axis=0))[0])
	print "{:} to find nodes, I used: {:}".format("-"*5, propsname[icol])
	inode = N.where(N.any(sprops[:, icol] != 0., axis=1))[0]
	print "%s n of nodes: %d" % ("-"*5, len(inode))

	return N.array(inode, dtype='int')
	

def proj_on_spropspc(pcArr, sprops):
	"""
	get the sprops in pc coordinates
	"""
	pccoefs = N.dot(sprops, pcArr)
	return pccoefs

	
def get_pixel_density(pc, obs, xr, yr, npix):
	xx = N.linspace(xr[0], xr[1], npix)
	yy = N.linspace(yr[0], yr[1], npix)
	#XX, YY = plt.meshgrid(xx, yy)
	dx = xx[1] - xx[0]
	dy = yy[1] - yy[0]
	## get the density in each pixel
	den = N.zeros((npix, npix))
	ixArr = N.array((pc - xr[0])/dx, dtype='int')
	iyArr = N.array((obs - yr[0])/dy, dtype='int')
	for ix, iy in zip(ixArr, iyArr):
		den[ix, iy] += 1
	#return XX, YY, den
	return den


def get_corrcoef_pccoef_sprops(pccoefs, sprops, corrfn):
	nsprops = sprops.shape[1]
	corrcoefArr = N.empty((nsprops, nsprops))
	for isprops in xrange(nsprops):
		for ipc in xrange(nsprops):
			corrcoef = corrfn(pccoefs[:, ipc], sprops[:, isprops])[0]
			corrcoefArr[isprops, ipc] = corrcoef
	return corrcoefArr

def plot_corr_pccoef_sprops_on_mesh(pccoefs, sprops, corrcoefArr, evals, propsname, figpp, npix=7):
	"""
	reference: matplotlib.sourceforge.net/examples/pylab_examples/multi_image.html
	get all images 
	"""
	nsprops = pccoefs.shape[1]
	fig = plt.figure(1)
	plt.clf()
	#camp = plt.cm.cool
	## plotting environment 
	coloraxes = [0.2, 0.04, 0.6, 0.01]
	cax = fig.add_axes(coloraxes)
	plotadjust = dict(hspace=0., wspace=0., left=0.1, right=0.95, bottom=0.12, top=0.97)
	width = (plotadjust['right'] - plotadjust['left'])/(nsprops*(1. + plotadjust['wspace']))
	height = (plotadjust['top'] - plotadjust['bottom'])/(nsprops*(1. + plotadjust['hspace']))
	print "%s width:%.2f, height:%.2f" % ("-"*5, width, height)
	ax = []
	images = []
	fontsize = 6
	labelx = -0.8
	nticks = 3
	#pcname = ('pc1', 'pc2', 'pc3', 'pc4', 'pc5')
	pcname = []
	for ii in xrange(nsprops):
		pcname.append("PC%d(%3.1f)" % (ii, evals[ii]))
	## intensity scale
	vmin = 1e5
	vmax = -1e5
	## xrange/yrange 
	xr = N.array([N.min(pccoefs, axis=0), N.max(pccoefs, axis=0)])
	yr = N.array([N.min(sprops, axis=0), N.max(sprops, axis=0)])

	for isprops in xrange(nsprops):
		for ipc in xrange(nsprops):
			## plot position
			pos = [plotadjust['left'] + ipc*(1. + plotadjust['wspace'])*width, plotadjust['bottom'] + isprops*(1. + plotadjust['hspace'])*height, width, height]
			a = fig.add_axes(pos)
			## plot ticks and labels
			a.yaxis.set_label_coords(labelx, 0.5)
			sx = []
			if isprops > 0:
				a.set_xticklabels([])
			else:
				a.xaxis.set_major_locator(plt_ticker.MaxNLocator(nticks))
				for xx in N.linspace(xr[0, ipc], xr[1, ipc], nticks):
					sx.append('{:.2f}'.format(xx))
				a.set_xticklabels(sx)
				for tick in a.xaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				a.set_xlabel(pcname[ipc], fontsize='xx-small')
			sy = []
			if ipc > 0:
				a.set_yticklabels([])
			else:
				for yy in N.linspace(yr[0, isprops], yr[1, isprops], nticks):
					sy.append('{:.2f}'.format(yy))
				#print sy
				a.yaxis.set_major_locator(plt_ticker.MaxNLocator(nticks))
				a.set_yticklabels(sy)
				for tick in a.yaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				a.set_ylabel(propsname[isprops][:6], fontsize='xx-small')

			## get the density
			den = get_pixel_density(pccoefs[:, ipc], sprops[:, isprops], xr[:, ipc], yr[:, isprops], npix)

			vmin = min(vmin, N.amin(den))
			vmax = max(vmax, N.amax(den))
			den_img = a.imshow(den, interpolation='nearest')
			#plt.setp(plt.gca(), 'ylim', reversed(plt.getp(plt.gca(), 'ylim')))
			a.set_ylim(a.get_ylim()[::-1])

			#cctext = '{:.2f}, {:.2f}, {:.2f}'.format(corrcoefArr[isprops, ipc], corrcoefArr[nobs+isprops, ipc], corrcoefArr[2*nobs+isprops, ipc])

			a.text(0.1, 0.1, "%.2f" % corrcoefArr[isprops, ipc], color='w', fontsize=6)
			#a.text(pos[0], pos[1], '{:d}-{:d}:{:.4f}'.format(ipc, isprops, corrcoefArr[isprops, ipc]), color='w')

			images.append(den_img)

			ax.append(a)

	class ImageFollower:
		"""
		reference: same as plot_corr_pc_n_obs_allclus
		update image in response to change in clim or camp on another image
		i.e. to get the same color scale for all images
		set the first image as the master, with all the others 
		"""
		def __init__(self, follower):
			self.follower = follower

		def __call__(self, leader):
			self.follower.set_cmap(leader.get_cmap())
			self.follower.set_clim(leader.get_clim())

	## set the first image as the master, with all the others
	norm = colors.Normalize(vmin=vmin, vmax=vmax)
	for i, im in enumerate(images):
		im.set_norm(norm)
		if i > 0:
			images[0].callbacksSM.connect('changed', ImageFollower(im))
	
	cbar = fig.colorbar(images[0], cax, orientation='horizontal') 
	cbar.ax.set_yticklabels(['%.0f'.format(N.min(den)), '%.0f'.format(N.max(den))], fontsize='xx-small')


	## following lines are neccessary only if we want to run this interactively and modify the color map
	## return the current axes to the first one
	plt.axes(ax[0]) 

	## because the current image must be in current axes
	plt.sci(images[0])

	figpp.savefig()


def plot_corr_pccoef_sprops_points(pccoefs, sprops, corrcoefArr, evals, propsname, figpp, npix=7):
	"""
	reference: matplotlib.sourceforge.net/examples/pylab_examples/multi_image.html
	get all images 
	"""
	nsprops = pccoefs.shape[1]
	fig = plt.figure(1)
	plt.clf()
	#camp = plt.cm.cool
	## plotting environment 
	plotadjust = dict(hspace=0., wspace=0., left=0.1, right=0.95, bottom=0.12, top=0.97)
	width = (plotadjust['right'] - plotadjust['left'])/(nsprops*(1. + plotadjust['wspace']))
	height = (plotadjust['top'] - plotadjust['bottom'])/(nsprops*(1. + plotadjust['hspace']))
	print "%s width:%.2f, height:%.2f" % ("-"*5, width, height)
	fontsize = 6
	labelx = -0.8
	nticks = 3
	#pcname = ('pc1', 'pc2', 'pc3', 'pc4', 'pc5')
	pcname = []
	for ii in xrange(nsprops):
		pcname.append("PC%d(%3.1f)" % (ii, evals[ii]))

	## xrange/yrange 
	xr = N.array([N.min(pccoefs, axis=0), N.max(pccoefs, axis=0)])
	yr = N.array([N.min(sprops, axis=0), N.max(sprops, axis=0)])

	for isprops in xrange(nsprops):
		for ipc in xrange(nsprops):
			## plot position
			pos = [plotadjust['left'] + ipc*(1. + plotadjust['wspace'])*width, plotadjust['bottom'] + isprops*(1. + plotadjust['hspace'])*height, width, height]
			a = fig.add_axes(pos)
			## plot ticks and labels
			a.yaxis.set_label_coords(labelx, 0.5)
			sx = []
			if isprops > 0:
				a.set_xticklabels([])
			else:
				a.xaxis.set_major_locator(plt_ticker.MaxNLocator(nticks))
				for xx in N.linspace(xr[0, ipc], xr[1, ipc], nticks):
					sx.append('{:.2f}'.format(xx))
				a.set_xticklabels(sx)
				for tick in a.xaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				a.set_xlabel(pcname[ipc], fontsize='xx-small')
			sy = []
			if ipc > 0:
				a.set_yticklabels([])
			else:
				for yy in N.linspace(yr[0, isprops], yr[1, isprops], nticks):
					sy.append('{:.2f}'.format(yy))
				#print sy
				a.yaxis.set_major_locator(plt_ticker.MaxNLocator(nticks))
				a.set_yticklabels(sy)
				for tick in a.yaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				a.set_ylabel(propsname[isprops][:6], fontsize='xx-small')

			plt.plot(pccoefs[:, ipc], sprops[:, isprops], 'ro', ms=1.2)

			a.text(0.1, 0.1, "%.2f" % corrcoefArr[isprops, ipc], color='w', fontsize=6)

	figpp.savefig()


def plot_props_hist(propsname, sprops, figname, nbin=20):
	
	figpp = PdfPages(figname)

	fig = plt.figure(2)
	nrow = 2
	ncol = 3

	for iprops, pname in enumerate(propsname):
		iplt = iprops
		if iprops % (nrow*ncol) == 0: 
			plt.clf()

		while iplt >= nrow*ncol:
			iplt -= nrow*ncol 

		ax = fig.add_subplot(nrow, ncol, iplt+1)
		h, bins, patches = plt.hist(sprops[:, iprops], bins=nbin)
		ax.set_xlabel(pname)
		if iplt % ncol == 0:
			ax.set_ylabel(r"$N_{cluster}$")

		ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(5))

		if iprops % (nrow*ncol) == nrow*ncol - 1:
			figpp.savefig()

	figpp.close()	


def choose_props(propsname, del_propArr):
	"""
	choose properties among scalar props 
	"""
	idelcol = N.array([], dtype='int')
	for del_prop in del_propArr:
		icol = N.where(propsname == del_prop)[0]
		idelcol = N.append(idelcol, icol)
	idelcol = N.sort(idelcol)
	print "{:} deleted columns: {:}".format("-"*5, idelcol)
	icolArr = N.arange(len(propsname))
	if len(idelcol) > 0:
		icolArr = N.delete(icolArr, idelcol)
	return icolArr


def get_stats_plots(evecs, evals, sprops, propsname, corrfn, figname, plotting_coeffs):
	"""
	calculate the PC coeffs, 
	calculate the corrln coeffs between PCs and the original properties and plot them
	"""
	pccoefs = proj_on_spropspc(evecs, sprops)
	corrcoefs = get_corrcoef_pccoef_sprops(pccoefs, sprops, corrfn)
	print "%s %s" % ("-"*5, "corrln coefficients")
	print ["{:s}: {:3.2f} {:3.2f} {:3.2f}".format(propsname[i], corrcoefs[i, 0], corrcoefs[i, 1], corrcoefs[i, 2]) for i in xrange(len(propsname))]

	if plotting_coeffs:
		figpp = PdfPages(os.path.join("corr_plots", "corr_pcs_%s_on_mesh.pdf" % figname))
		plot_corr_pccoef_sprops_on_mesh(pccoefs, sprops, corrcoefs, evals, propsname, figpp)
		figpp.close()

	return pccoefs, corrcoefs


def get_data(take_fraction, logall, logmass, standardized, nodeclustersonly, noneed):
	"""
	get the data and rearrange for use
	"""
	fname = os.path.join("..", "cluster_scalars_z0.1_frac{:}.dat".format(take_fraction))
	#fname = os.path.join("..", "cluster_scalars_z0.1_frac{:}_changed.dat".format(take_fraction))
	print "%s read %s" % ("-"*5, fname)
	propsname = get_propsname(fname)
	sprops = N.loadtxt(fname)

	## take only properties
	sprops = sprops[:, 2:]
	propsname = propsname[2:]

	## when delete the columns, makesure that nplane is the right number
	if len(noneed) > 0:
		icols = choose_props(propsname, noneed)
		print "chosen columns: {:}".format(icols)
		sprops = sprops[:, icols]
		propsname = propsname[icols]
		print "the original props I'm going to use"
		print ["{:s}: {:3.2f}".format(propsname[i], sprops[0, i]) for i in xrange(len(propsname))]

	## followings shouldn't be true at the same time
	if logall:
		sprops = take_logall(sprops)
	elif logmass:
		iplane, sprops = take_logmass(propsname, sprops)

	if standardized:
		sprops = calc_fluctuation(sprops)

	print "the converted props are"
	print ["{:s}: {:}".format(propsname[i], sprops[0, i]) for i in xrange(len(propsname))]

	if nodeclustersonly:
	## for only nodes
		print "%sconsider only node clusters%s" % ("="*5, "="*5)
		inode = get_node_clus(sprops, propsname)
		sprops = sprops[inode, :]
	else:
		print "%sconsider all clusters(include non-nodes)%s" % ("="*5, "="*5)

	return sprops, propsname


def do_pca(ifile=8, take_fraction=True, logall=False, logmass=False, standardized=True, plotting_coeffs=True, plotting_hist=True, nodeclustersonly=True, noneed=['FrGRP-AGR', 'FrAFP-AFM']):

	"""
	check the data which I printed out between plane mass or mass fraction in plane, 
	also which matrix is used to do pca
	"""
	## get data
	sprops, propsname = get_data(take_fraction, logall, logmass, standardized, nodeclustersonly, noneed)

	name_base = "standardized{:}_sprops_nodesonly_logall{:}_logmass{:}_nprops{:}".format(standardized, logall, logmass, len(propsname))

	## do pca
	evals, evecs = get_pcs(sprops, propsname)

	pccoefs, corrcoefs = get_stats_plots(evecs, evals, sprops, propsname, st.pearsonr, name_base, plotting_coeffs)

	if plotting_hist:
		figname = "hist_pcs_standardized{:}_sprops_nodesonly_logall{:}_logmass{:}_nprops{:}.pdf".format(standardized, logall, logmass, len(propsname))
		plot_props_hist(propsname, sprops, figname)

	return pccoefs, corrcoefs, sprops, evals, evecs

