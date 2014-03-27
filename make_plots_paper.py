## to plot for paper figures

import os
import numpy as N
import matplotlib.pylab as plt
import matplotlib.ticker as plt_ticker
import reading_fn as rfn


def get_parameters():
	obsname = (r'$N_{red}$', r'$N_{phase}$', r'$Y(SZ)$', r'$WL$', r'$V_{Phase}$')
	nlos = 96
	nobs = 5
	return obsname, nlos, nobs
	

class MassScatterWidth:
	"""
	scatter width plot
	"""

	def __init__(self):
		self.adjustparam = dict(hspace=0.35, wspace=0.5, left=0.08, right=0.95, bottom=0.05, top=0.97)
		self.obsname, self.nlos, self.nobs = get_parameters()

		def plot_scatter_width
		


class MassCorrScatter:
	"""
	make scatter plot
	"""

	def __init__(self, take_log=False, use_scatter=True):
		self.adjustparam = dict(hspace=0.35, wspace=0.5, left=0.08, right=0.95, bottom=0.05, top=0.97)
		self.obsname, self.nlos, self.nobs = get_parameters()
		self.take_log = take_log
		self.use_scatter = use_scatter

	def get_valid_array_altogether(self, m_est, trueM):
		truefalse_all = N.ones(len(m_est))
		ikeep = N.nonzero(N.all(m_est > 1.e5, axis=1))[0]
		if len(ikeep) != len(m_est):
			print "invalid mass exists"
		if self.take_log:
			m_est = N.log10(m_est[ikeep, :])
			masstype = "log"
		elif self.use_scatter:
			m_est = (m_est[ikeep, :] - trueM)/trueM
			masstype = "dtr"
		else:
			m_est = m_est[ikeep, :]/1.e14
			masstype = "org"
		return m_est, masstype

	def plot_scatter(self, iclus):
		observed = rfn.read_data()
		fig = plt.figure(2)
		plt.clf()

		trueM = observed[iclus*self.nlos, 1]
		mass = observed[iclus*self.nlos:(iclus+1)*self.nlos, 2:]
		mass, masstype = self.get_valid_array_altogether(mass, trueM)

		ii = 0
		for iobs in xrange(self.nobs):
			for jobs in xrange(iobs+1, self.nobs, 1):

				ax = fig.add_subplot(4, 3, ii+1)
				ax.plot(mass[:, iobs], mass[:, jobs], 'bo', ms=2.5)

				plt.subplots_adjust(**self.adjustparam)
				ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(3))
				ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(3))

				#ax.set_xlabel("%s/<%s>" % (self.obsname[iobs], self.obsname[iobs])) 
				#ax.set_ylabel("%s/<%s>" % (self.obsname[jobs], self.obsname[jobs])) 

				fontsize=9
				for tick in ax.xaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				for tick in ax.yaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)

				ii += 1
				print self.obsname[iobs], self.obsname[jobs]

		plt.savefig(os.path.join("paper", "figure", "scatter_%s_clus%d.eps" % (masstype, iclus)), orientation='portrait', transparent=True)


class MassCorrCoefHist:
	"""
	plot histogram of correlation coeff/covariance 
	"""
	def __init__(self, datafilename = "dtr_orgxiArr"):
		self.datafilename  = datafilename
		self.nobs = 5
		self.adjustparam = dict(hspace=0.3, wspace=0.3, left=0.05, right=0.95, bottom=0.05, top=0.97)
		self.obsname = (r'$N_{red}$', r'$N_{phase}$', r'$Y(SZ)$', r'$WL$', r'$V_{Phase}$')

	def plot_corrcoef(self, icor):
		"""
		icor: which correlation coefficient will be plotted
		icor 0:pearsonr, 1:spearman, 2:kendalltau, 3:covariance
		"""
		xi_arr = N.loadtxt(os.path.join("corr_scatter", "scatterplot", "data", "%s.dat" % self.datafilename))
		ncorcoef = xi_arr.shape[1]/(self.nobs*(self.nobs-1)/2.)
		nbin = 10

		fig = plt.figure(1)
		fig.clf()

		ii = 0
		for iobs in xrange(self.nobs):
			for jobs in xrange(iobs+1, self.nobs, 1):
				eachval = xi_arr[:, ncorcoef*ii + icor]
				inotnan = N.where(N.isnan(eachval)==False)[0]
				eachval = eachval[inotnan]

				ax = fig.add_subplot(4, 3, ii+1)
				h, bins, patches = plt.hist(eachval, bins=nbin, facecolor='grey', edgecolor='k')
				plt.subplots_adjust(**self.adjustparam)
				ax.xaxis.set_major_locator(plt_ticker.MaxNLocator(3))
				ax.yaxis.set_major_locator(plt_ticker.MaxNLocator(4))
				
				if N.int(ii/3) == 2:
					ax.set_xlabel(r"$\xi$")
				if ii % 3 == 0:
					ax.set_ylabel(r"$N_{cluster}$")

				fontsize = 9
				for tick in ax.xaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)
				for tick in ax.yaxis.get_major_ticks():
					tick.label1.set_fontsize(fontsize)

				ii += 1
				print self.obsname[iobs], self.obsname[jobs]

		plt.savefig(os.path.join("paper", "figure", "hist_%s.eps" % self.datafilename), orientation='portrait', transparent=True)


if __name__ == "__main__":

	## scatter plot
	paperfig2 = MassCorrScatter()
	paperfig2.plot_scatter(0)

	paperfig3 = MassCorrCoefHist()
	paperfig3.plot_corrcoef(0)
