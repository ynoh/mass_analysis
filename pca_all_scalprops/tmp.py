def plot_corr_pc_n_obs_allclus(pcall, obsall, nobs, npix, corrcoefArr, figname):
	"""
	reference: matplotlib.sourceforge.net/examples/pylab_examples/multi_image.html
	get all images 
	"""
	fig = plt.figure(1)
	#camp = plt.cm.cool
	## plotting environment 
	obsname = (r'$N_{red}$', r'$N_{phase}$', r'$Y(SZ)$', r'$WL$', r'$V_{ph}$')
	pcname = ('pc1', 'pc2', 'pc3', 'pc4', 'pc5')
	coloraxes = [0.2, 0.06, 0.6, 0.03]
	cax = fig.add_axes(coloraxes)
	plotadjust = dict(hspace=0., wspace=0., left=0.1, right=0.95, bottom=0.15, top=0.95)
	width = (plotadjust['right'] - plotadjust['left'])/(nobs*(1. + plotadjust['wspace']))
	height = (plotadjust['top'] - plotadjust['bottom'])/(nobs*(1. + plotadjust['hspace']))
	ax = []
	images = []
	fontsize = 6
	## intensity scale
	vmin = 1e5
	vmax = -1e5
	## xrange/yrange 
	xr = N.array([N.min(pcall, axis=0), N.max(pcall, axis=0)])
	yr = N.array([N.min(obsall, axis=0), N.max(obsall, axis=0)])
	nticks = 4

	for iobs in xrange(nobs):
		for ipc in xrange(nobs):
			## plot position
			pos = [plotadjust['left'] + ipc*(1. + plotadjust['wspace'])*width, plotadjust['bottom'] + iobs*(1. + plotadjust['hspace'])*height, width, height]
			a = fig.add_axes(pos)
			## plot ticks and labels
			sx = []
			if iobs > 0:
				a.set_xticklabels([])
			else:
				a.xaxis.set_major_locator(plt_ticker.MaxNLocator(nticks))
				for xx in N.linspace(xr[0, ipc], xr[1, ipc], nticks):
					sx.append('{:.2f}'.format(xx))
				a.set_xticklabels(sx)
				for tick in a.xaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				a.set_xlabel(pcname[ipc], fontsize='small')
			sy = []
			if ipc > 0:
				a.set_yticklabels([])
			else:
				for yy in N.linspace(yr[0, iobs], yr[1, iobs], nticks):
					sy.append('{:.2f}'.format(yy))
				print sy
				a.yaxis.set_major_locator(plt_ticker.MaxNLocator(nticks))
				a.set_yticklabels(sy)
				for tick in a.yaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				a.set_ylabel(obsname[iobs], fontsize='small')

			den = get_pixel_density(pcall[:, ipc], obsall[:, iobs], xr[:, ipc], yr[:, iobs], npix)
			cctext = '{:.2f}, {:.2f}, {:.2f}'.format(corrcoefArr[iobs, ipc], corrcoefArr[nobs+iobs, ipc], corrcoefArr[2*nobs+iobs, ipc])

			a.text(pos[0], pos[1], cctext, color='w', fontsize='xx-small')
			#a.text(pos[0], pos[1], '{:d}-{:d}:{:.4f}'.format(ipc, iobs, corrcoefArr[iobs, ipc]), color='w')

			vmin = min(vmin, N.amin(den))
			vmax = max(vmax, N.amax(den))
			images.append(a.imshow(den))

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
	#cbar.ax.set_yticklabels(['%.0f'.format(N.min(den)), '%.0f'.format(N.max(den))], fontsize='xx-small')


	## following lines are neccessary only if we want to run this interactively and modify the color map
	## return the current axes to the first one
	plt.axes(ax[0]) 

	## because the current image must be in current axes
	plt.sci(images[0])

	plt.savefig(figname)
	#plt.show()

