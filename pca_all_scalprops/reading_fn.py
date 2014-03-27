"""
functions frequently used to read data
read_principal_component is changed (2011/11/28)

converted mass function is added (2012/02/??)

m_oper becomes p_like following the notation in the paper
i.e. m_oper was nxn matrix but now it is nx1 matrix

2012/03/31 read_data: the columns of observed data is changed
		*_old  is removed

functions: the order is bottom to top
	get_filename(nbox):
	def get_estimated_mtrue(pcfilename = "alllos_sca_cov_1.00e+14", observed=[], nlos=96, nobs=5):
	get_mlh_op_entire_clus(observed, mtype, nlos=96, nobs=5, weight=True):
	get_mlh_operator(mobs_valid_arr, weight=True):
	get_moff_corrected_arr(observed_per_clus, m_ave_tot):
	project_on_pc(m_est_perclus, evecs, pccorrfnname, normalize=True):
	change_n_to_p(pc):
	convert_mass(observed, nlos, mtype):
	get_valid_array_altogether_alllos(observed, nobs, nlos, mtype="org", m_oper=[], m_off="none"):
	get_valid_array_altogether(m_est_perclus, m_true, mtype, m_oper=N.array([]), m_off="none"):
	get_m_off_allclus(observed, nlos, mtype):
	calc_maxlikelihoodmass(observed, nlos, mtype):
	calc_ave_true_mass(observed, nlos):
	get_obsname():
	read_data(takeout=True, setnan0=True):
	set_nan_to_zero(observed_mass, nobs):
	take_out_los_having_clus(observed_mass, flagval):
	read_principal_component(filename):

"""

import os, sys
import numpy as N
#sys.path.append(os.path.join("corr_scatter", "get_pc"))
#import get_evec


#res_abspath = os.path.abspath("/Users/yookyung/Research/")
res_abspath = os.path.abspath(os.path.join("..", "..", ".."))
obs_path = os.path.join(res_abspath, "filaments/observable/")
pc_path = os.path.join(res_abspath, "mass_analysis/corr_scatter/get_pc/eigen/")

def read_principal_component(filename):
	filenamehead = "eigenvalvecs"
	filename = "%s_%s.dat" % (filenamehead, filename)
	print "read %s" % filename
	data = N.loadtxt(os.path.join(pc_path, filename), comments='#')
	if filename.find("alllos") == -1:
		eval = data[:, 1]
		evec = data[:, 2:]
	else:
		eval = data[:, 0]
		evec = data[:, 1:]
	return eval, evec


def take_out_los_having_clus(observed_mass, flagval):
	"""
	set mass = 0 for that line-of-sight
	"""
	iflag1 = N.where(flagval > 0)[0]
	#print N.column_stack((iflag1, flagval[iflag1]))
	observed_mass[iflag1, :] = 0.
	return observed_mass


def set_nan_to_zero(observed_mass, nobs):
	"""
	change nan to zero to take out that value easily
	"""
	for iobs in xrange(nobs):
		inan = N.where(N.isnan(observed_mass[:, iobs]))[0]
		observed_mass[inan, iobs] = 0.
	
	return observed_mass

def read_data(takeout=True, setnan0=True):
	"""
	read data not from the saved file but calculate everytime
	if takeout = True, flagval = 1 should be taken out
	the line-of-sight where flagval = 1 has a big cluster along there so the mass is too big
	it reads original id
	column info: "mass id", "original id","nodeM", "nx", "ny", "nz", "flag", "RedRi", "PhaseRi", "SZ", "3SigVdsp", "BiviVdsp", "WeakL"
	2012-03-31: the columns of vel dsp and wl are swapped and the function returns obsname 
	"""
	sys.path.append(os.path.join(obs_path, "observedmass"))
	print obs_path
	import observedmass_notsave as calc_obs
	reload(calc_obs)
	observed = calc_obs.get_mass()
	
	## swap weak lensing & 3sig-clipping
	observed[:, N.array([9, 11])] = observed[:, N.array([11, 9])]

	## before delete keep flagval to take out the line-of-sight along which there is a cluster 
	flagval = observed[:, 5]

	## delete no need columns, including 3-sig clipping method
	observed = N.delete(observed, N.s_[2:6:], 1)
	observed = observed[:, :-1] # remove 3-sig clipping

	## 2012-03-31 to match the order with Joanne,
	## WL and vel dsp. columns are swapped
	observed[:, [-2, -1]] = observed[:, [-1, -2]]

	nobs = 5
	obsname = [r"$M_{N_{\rm red}}$", r"$M_{N_{\rm ph}}$", r"$M_{\rm SZ}$", r"$M_{\rm Vel}$", r"$M_{\rm WL}$"]

	## takeout == True, set all mass be 0 for the line-of-sight along which there is a cluster
	if takeout:
		print "the mass observed along the line-of-sight where a massive cluster is set to be 0."
		observed[:, 2:] = take_out_los_having_clus(observed[:, 2:], flagval)

	if setnan0:
		print "nan value will be set to 0."
		observed[:, 2:] = set_nan_to_zero(observed[:, 2:], nobs)

	print "Finished getting observed mass array: {:}".format(observed.shape)

	return observed, obsname

def get_obsname():
	return [r"$M_{N_{\rm red}}$", r"$M_{N_{\rm ph}}$", r"$M_{\rm SZ}$", r"$M_{\rm Vel}$", r"$M_{\rm WL}$"]



def calc_ave_true_mass(observed, nlos):
	nclus = len(observed)/nlos
	mtrue_ave = 0.
	nvalid_los = 0.
	for iclus in xrange(nclus):
		ikeep = N.nonzero(N.all(observed[iclus*nlos:(iclus+1)*nlos, 2:] > 1.e5, axis=1))[0]
		mtrue_ave += observed[iclus*nlos, 1]*len(ikeep)
		nvalid_los += len(ikeep)
	mtrue_ave /= 1.*nvalid_los
	return mtrue_ave


def get_m_off_allclus(observed, nlos, mtype):
	if mtype == "mlo":
		m_obs_ave_tot = N.average(observed[:, 2:], axis=0)
		m_true_ave_tot = calc_ave_true_mass(observed, nlos)
		print "average mobs:", m_obs_ave_tot
		print "true mass average", m_true_ave_tot
		m_off = m_obs_ave_tot - m_true_ave_tot
		m_off.reshape(1, len(m_off))
	else:
		m_off = "none"
	return m_off


def get_valid_array_altogether(m_est_perclus, m_true, mtype, m_oper=N.array([]), m_off="none"):
	"""
	m_est_perclus = (96 x nobs)
	drop invalid data's line-of-sight from all observations
	take log or not 
	m_est_perclus array itself is modified
	log: take log of mass
	org: use observed mass or the mass divided by 1.e14
	sca: (observed mass - FoF mass)/FoF mass
	obs: M^{each line of sight}_{iobs} - <M>_{for all measurements}^{each line of sight}
	mlh: use max likelihood mass, m_oper is necessary for this case
	if m_off should be m_obs_ave_tot - m_true_ave_tot vector for mtype=mlo case.
	2012/03/31: valid_array_altogether_ikeep is combined
	"""
	ikeep = N.nonzero(N.all(m_est_perclus > 1.e5, axis=1))[0]

	## drop invalid data row
	m_est_perclus = m_est_perclus[ikeep, :]

	## want to use log?
	if mtype == "org":
		pass
	elif mtype == "log":
		m_est_perclus = N.log10(m_est_perclus)
	elif mtype == "sca":
		m_est_perclus = (m_est_perclus - m_true)/m_true
	elif mtype == "obs":
		obs_ave = N.average(m_est_perclus, axis=1) ## 96x1
		#	m_est_perclus_perobs = (m_est_perclus_perobs - obs_ave)/obs_ave
		for i in xrange(m_est_perclus.shape[1]):
			m_est_perclus[:, i] = (m_est_perclus[:, i] - obs_ave)/obs_ave
	elif mtype == "mlh":
		mlhmass = N.zeros(len(m_est_perclus))
		for i, mrow in enumerate(m_est_perclus):
			mlhmass[i] = N.dot(mrow, m_oper)
		if m_true > 1.e15:
			print mlhmass
		for i in xrange(m_est_perclus.shape[1]):
			m_est_perclus[:, i] = (m_est_perclus[:, i] - mlhmass)/mlhmass
	elif mtype == "mlo": ## same as mlh but subtract offset
		mlhmass = N.zeros(len(m_est_perclus))
		for i, mrow in enumerate(m_est_perclus):
			mlhmass[i] = N.dot(mrow - m_off, m_oper)
		for i in xrange(m_est_perclus.shape[1]):
			m_est_perclus[:, i] = (m_est_perclus[:, i] - mlhmass)/mlhmass
	elif mtype == "ml2":
		m_oper = get_mlh_operator(m_est_perclus)
		mlhmass = N.zeros(len(m_est_perclus))
		for i, mrow in enumerate(m_est_perclus):
			mlhmass[i] = N.dot(mrow, m_oper)
		for i in xrange(m_est_perclus.shape[1]):
			m_est_perclus[:, i] = (m_est_perclus[:, i] - mlhmass)/mlhmass
	elif mtype == "mo2":
		m_oper = get_mlh_operator(m_est_perclus)
		m_ave = N.average(m_est_perclus, axis=0)
		m_off = m_ave - m_true
		mlhmass = N.zeros(len(m_est_perclus))
		for i, mrow in enumerate(m_est_perclus):
			mlhmass[i] = N.dot(mrow - m_off, m_oper)
		for i in xrange(m_est_perclus.shape[1]):
			m_est_perclus[:, i] = (m_est_perclus[:, i] - mlhmass)/mlhmass
	elif mtype == "mpc":
		print "this mtype(mpc) is valid for only all los analysis" 
	else:
		print "mtype: %s" % mtype
		print "this mtype is wrong or altogether function doesn't take care of this mtype"
		raise
	return m_est_perclus, ikeep


def get_valid_array_altogether_alllos(observed, nobs, nlos, mtype="org", m_oper=[], m_off="none", pcfilename="alllos_sca_cov_1.00e+14"):
	if mtype.find("mpc") == -1:
		mobs_valid_arr = N.empty((len(observed), nobs))
		nclus = len(observed)/nlos
		ikeep_arr = N.zeros(len(observed), dtype='int')
		ivalid_los_arr = N.zeros(nclus)
		ivalid_los = 0

		for iclus in xrange(nclus):
			mobs_valid, ikeep = get_valid_array_altogether(observed[iclus*nlos:(iclus+1)*nlos, 2:], observed[iclus*nlos, 1], mtype, m_oper=m_oper, m_off=m_off)
			mobs_valid_arr[ivalid_los:ivalid_los+len(mobs_valid), :] = mobs_valid
			ikeep_arr[ivalid_los:ivalid_los+len(mobs_valid)] = ikeep
			ivalid_los_arr[iclus] = len(mobs_valid)
			ivalid_los += len(mobs_valid)

		mobs_valid_arr = mobs_valid_arr[:ivalid_los, :]
		ikeep_arr = ikeep_arr[:ivalid_los]
	else:
		mtrue_est, mvalid, ivalid_los_arr, ikeep_arr = get_estimated_mtrue(pcfilename=pcfilename, observed=observed, nlos=nlos, nobs=nobs)
		mobs_valid_arr = N.zeros_like(mvalid)
		for iobs in xrange(nobs):
			mobs_valid_arr[:, iobs] = mvalid[:, iobs]/mtrue_est
	return mobs_valid_arr, ivalid_los_arr, ikeep_arr


def convert_mass(observed, nlos, mtype):
	"""
	convert the mass to the appropriate form for the purpose
	this is different from get_valid_arr_altogether
	this function handles all clusters and non-valid los
	non-valid los has still 0 for the observed mass
	invalid line-of-sight is not excluede
	"""
	fw1 = open(os.path.join(res_abspath, "mass_analysis", "mtype-%s.dat" % mtype), "w")
	m_conv = N.zeros_like(observed[:, 2:])
	nclus = len(observed)/nlos
	if mtype == "log":
		m_conv = N.log10(observed[:, 2:] + 1.e-7)
	elif mtype == "sca":
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:] 
			mtrue1 = observed[iclus*nlos, 1]
			m_conv[iclus*nlos:(iclus+1)*nlos, :] = (mobs1 - mtrue1)/mtrue1
	elif mtype == "obs":
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:] 
			mobs_ave = N.average(mobs1, axis=1) ## nlos x 1
			clusnum = N.zeros(nlos)
			clusnum.fill(iclus)
			truemass = N.zeros(nlos)
			truemass.fill(observed[iclus*nlos, 1])
			N.savetxt(fw1, N.column_stack((clusnum, truemass, mobs_ave)), fmt="%5d\t%.4e\t%.4e")
			for i in xrange(m_conv.shape[1]):
				m_conv[iclus*nlos:(iclus+1)*nlos, i]  = (mobs1[:, i] - mobs_ave)/(mobs_ave + 1.e-7) ## to prevent from dividing by 0
	elif mtype == "mlh":
		ikeep = N.nonzero(N.all(observed[:, 2:] > 1.e5, axis=1))[0]	
		m_oper = get_mlh_operator(observed[ikeep, 2:])
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:] 
			mlhmass = N.zeros(len(mobs1))
			for i, mrow in enumerate(mobs1):
				mlhmass[i] = N.dot(mrow, m_oper)
			#if iclus == 0:
			#	print mlhmass
			clusnum = N.zeros(nlos)
			clusnum.fill(iclus)
			truemass = N.zeros(nlos)
			truemass.fill(observed[iclus*nlos, 1])
			N.savetxt(fw1, N.column_stack((clusnum, truemass, mlhmass)), fmt="%5d\t%.4e\t%.4e")
			for i in xrange(mobs1.shape[1]):
				m_conv[iclus*nlos:(iclus+1)*nlos, i] = (mobs1[:, i] - mlhmass)/(mlhmass + 1.e-7)
	elif mtype == "ml2":
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:] 
			ikeep = N.nonzero(N.all(mobs1 > 1.e5, axis=1))[0]
			m_oper = get_mlh_operator(mobs1[ikeep, :]) ## covariance matrix should be calculated by only using the valid mass
			mlhmass = N.zeros(len(mobs1))
			for i, mrow in enumerate(mobs1):
				mlhmass[i] = N.dot(mrow, m_oper)
			clusnum = N.zeros(nlos)
			clusnum.fill(iclus)
			truemass = N.zeros(nlos)
			truemass.fill(observed[iclus*nlos, 1])
			N.savetxt(fw1, N.column_stack((clusnum, truemass, mlhmass)), fmt="%5d\t%.4e\t%.4e")
			#N.savetxt(fw1, mlhmass, fmt="%.4e")
			for i in xrange(mobs1.shape[1]):
				m_conv[iclus*nlos:(iclus+1)*nlos, i] = (mobs1[:, i] - mlhmass)/(mlhmass + 1.e-7)
	elif mtype == "mlo": ## same as mlh but subtract offset
		ikeep = N.nonzero(N.all(observed[:, 2:] > 1.e5, axis=1))[0]	
		m_oper = get_mlh_operator(observed[ikeep, 2:])
		m_off = get_m_off_allclus(observed, nlos, mtype)
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:]
			mlhmass = N.zeros(len(mobs1))
			for i, mrow in enumerate(mobs1):
				mlhmass[i] = N.dot(mrow - m_off, m_oper)
			clusnum = N.zeros(nlos)
			clusnum.fill(iclus)
			truemass = N.zeros(nlos)
			truemass.fill(observed[iclus*nlos, 1])
			N.savetxt(fw1, N.column_stack((clusnum, truemass, mlhmass)), fmt="%5d\t%.4e\t%.4e")
			for i in xrange(mobs1.shape[1]):
				m_conv[iclus*nlos:(iclus+1)*nlos, i] = (mobs1[:, i] - mlhmass)/(mlhmass + 1.e-7)
	elif mtype == "mo2":
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:]
			ikeep = N.nonzero(N.all(mobs1 > 1.e5, axis=1))[0]
			m_oper = get_mlh_operator(mobs1[ikeep, :])
			m_ave = N.average(mobs1[ikeep, :], axis=0)
			m_off = m_ave - observed[iclus*nlos, 1]
			mlhmass = N.zeros(len(mobs1))
			for i, mrow in enumerate(mobs1):
				mlhmass[i] = N.dot(mrow - m_off, m_oper)
			clusnum = N.zeros(nlos)
			clusnum.fill(iclus)
			truemass = N.zeros(nlos)
			truemass.fill(observed[iclus*nlos, 1])
			N.savetxt(fw1, N.column_stack((clusnum, truemass, mlhmass)), fmt="%5d\t%.4e\t%.4e")
			for i in xrange(mobs1.shape[1]):
				m_conv[iclus*nlos:(iclus+1)*nlos, i] = (mobs1[:, i] - mlhmass)/(mlhmass + 1.e-7)
	else:
		m_conv = observed[:, 2:].copy()
	fw1.close()
	return m_conv


def change_n_to_p(pc):
	if pc[0] < 0.:
		pc *= -1.
	return pc



def project_on_pc(m0, ev, pccorrfnname, normalize=True):
	"""
	Basically, if not normalized, it's the length on the direction of that pc
	if normalized, the angle between the observed mass and the pc direction
	copied from corr_losDotprops_masspc/get_corr_LosDotProps_masspc.py
	"""
	m_est_perclus = m0.copy()
	evecs = ev.copy()
	proj_mest = N.empty_like(m_est_perclus)
	nobs = m_est_perclus.shape[1]

	if pccorrfnname == "cov":
		for iobs in xrange(nobs):
			m_est_perclus[:, iobs] = m_est_perclus[:, iobs] - N.average(m_est_perclus[:, iobs])
			if N.isnan(N.any(proj_mest[:, iobs])):
				print "evecs[:, %d]" % iobs
				print evecs[:, iobs]
	else:
		for iobs in xrange(nobs):
			m_est_perclus[:, iobs] = (m_est_perclus[:, iobs] - N.average(m_est_perclus[:, iobs]))/N.std(m_est_perclus[:, iobs])

	for iobs in xrange(nobs):
		evecs[:, iobs] = change_n_to_p(evecs[:, iobs])

	## sum up to get the coefficients for the roated coordinate
	if normalize:
		for iobs in xrange(nobs):
			m_est_perclus[:, iobs] /= N.sqrt(N.sum(m_est_perclus**2, axis=1))

	for iobs in xrange(nobs):
		proj_mest[:, iobs] = N.sum(evecs[:, iobs]*m_est_perclus, axis=1)

	return proj_mest



def get_moff_corrected_arr(observed_per_clus, m_ave_tot):
	## this function is used for calc_maxlikelihoodmass
	## subtract offset from the observed mass
	mobs = observed_per_clus[:, 2:]
	ikeep = N.nonzero(N.all(mobs > 1.e5, axis=1))[0]
	m_true = observed_per_clus[0, 1]
	moff = m_ave_tot - m_true
	#print m_ave_tot, m_true
	#print moff
	mobs += moff
	mobs = mobs[ikeep, :]
	return mobs

	
def calc_maxlikelihoodmass(observed, nlos, mtype):
	## this uses the covariance matrix of the entire sample
	## to avoid round up error divide the mass by 1.e14
	## equation shows that dividing the covariance matrix by constant doesn't affect the final results
	## invalid line-of-sight is not excluede
	nclus = len(observed)/nlos
	mlhmass = N.zeros(len(observed))
	if mtype == "mlh":
		ikeep = N.nonzero(N.all(observed[:, 2:] > 1.e5, axis=1))[0]	
		m_oper = get_mlh_operator(observed[ikeep, 2:])
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:] 
			for i, mrow in enumerate(mobs1):
				mlhmass[iclus*nlos+i] = N.dot(mrow, m_oper)
	elif mtype == "ml2":
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:] 
			ikeep = N.nonzero(N.all(mobs1 > 1.e5, axis=1))[0]
			m_oper = get_mlh_operator(mobs1[ikeep, :]) ## covariance matrix should be calculated by only using the valid mass
			for i, mrow in enumerate(mobs1):
				mlhmass[iclus*nlos+i] = N.dot(mrow, m_oper)
	elif mtype == "mlo": ## same as mlh but subtract offset
		ikeep = N.nonzero(N.all(observed[:, 2:] > 1.e5, axis=1))[0]	
		m_oper = get_mlh_operator(observed[ikeep, 2:])
		m_off = get_m_off_allclus(observed, nlos, mtype)
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:]
			for i, mrow in enumerate(mobs1):
				mlhmass[iclus*nlos+i] = N.dot(mrow - m_off, m_oper)
	elif mtype == "mo2":
		for iclus in xrange(nclus):
			mobs1 = observed[iclus*nlos:(iclus+1)*nlos, 2:]
			ikeep = N.nonzero(N.all(mobs1 > 1.e5, axis=1))[0]
			m_oper = get_mlh_operator(mobs1[ikeep, :])
			m_ave = N.average(mobs1[ikeep, :], axis=0)
			m_off = m_ave - observed[iclus*nlos, 1]
			for i, mrow in enumerate(mobs1):
				mlhmass[iclus*nlos+i] = N.dot(mrow - m_off, m_oper)
	return mlhmass


def get_mlh_operator(mobs_valid_arr, weight=True):
	"""
	calculate
	\sum_i C^{-1}_ij/ \sum_i,j C^{-1}_ij
	C^{-1} is the inverse of covariance matrix
	this covariance matrix can be calculated using the entire clusters line-of-sight
	and each cluster's 
	"""
	#print "weight=%s" % weight
	#massArr = mobs_valid_arr/1.e14
	massArr = mobs_valid_arr
	covmat = N.cov(massArr.T)
	inv_covmat = N.linalg.inv(covmat)
	m_oper = N.sum(inv_covmat, axis=1)
	if weight:
		m_oper /= N.sum(inv_covmat)
	return m_oper


def get_mlh_op_entire_clus(observed, mtype, nlos=96, nobs=5, weight=True):
	if mtype == "mlh" or mtype == "mlo":
		mobs_valid_arr_all, ivalid_los_arr, ikeep_arr = get_valid_array_altogether_alllos(observed, nobs, nlos)
		m_oper = get_mlh_operator(mobs_valid_arr_all, weight=weight)
		return m_oper
	else:
		return []


def get_estimated_mtrue(pcfilename = "alllos_sca_cov_1.00e+14", observed=[], nlos=96, nobs=5):
	"""
	implemented from corr_scatter/calc_mtrue_estimate.py

	"""
	if len(observed) == 0:
		observed, obsname = rfn.read_data()
	nclus = len(observed)/nlos
	mvalid, ivalid, ikeep = get_valid_array_altogether_alllos(observed, nobs, nlos)
	evals, evecs = read_principal_component(pcfilename)

	## ikeep_id contains the index of observed array (c.f. ikeep points the index of each observable)
	ikeep_id = N.zeros_like(ikeep)
	id = 0
	for iclus in xrange(nclus):
		ikeep_id[id:id+ivalid[iclus]] = iclus*nlos + ikeep[id:id+ivalid[iclus]]
		id += ivalid[iclus]

	mtrue = observed[ikeep_id, 1]

	## mest to mtrue ratio  
	mestTOmtrue = N.zeros_like(mvalid)
	for iobs in xrange(nobs):
		mestTOmtrue[:, iobs] = mvalid[:, iobs]/mtrue

	mestTOmtrue_ave = N.average(mestTOmtrue, axis=0)

	denom = N.dot(mestTOmtrue_ave, evecs[:, -1])
	mtrue_est = N.zeros(len(mvalid))
	for i, mest in enumerate(mvalid):
		mtrue_est[i] = N.dot(evecs[:, -1], mest)

	mtrue_est /= denom

	return mtrue_est, mvalid, ivalid, ikeep

	

def get_filename(nbox):
	"""
	log: take log of mass
	org: use observed mass or the mass divided by 1.e14
	sca: (observed mass - FoF mass)/FoF mass
	obs: M^{each line of sight}_{iobs} - <M>_{for all measurements}^{each line of sight}
	mlh: max likelihood mass using the full sample cov. w/o correcting the offset
	ml2: max likelihood mass using the indv cluster cov. w/o correting the offset
	mlo: max likelihood mass using the full sample cov. w/ correcting the offset
	mo2: max likelihood mass using the indv cluster cov. w/ correcting the offset
	mpc: full los case.. using the mass calculated in get_estimated_mtrue - M_4 in the paper definition
	"""
	filenameArr = N.array([])
	for ll in ['log', 'org', 'sca', 'obs', 'mlh', 'ml2', 'mlo', 'mo2', 'mpc']:
		for outl in ['incl', 'excl']:
			for co in ['cov', 'pea', 'spe', 'ken']:
				filename = "%s_%soutliers_%s" % (ll, outl, co)
				if outl[0] == 'e':
					filename = filename + '_nbox{:}'.format(nbox) 
				filenameArr = N.append(filenameArr, filename)
	return filenameArr
