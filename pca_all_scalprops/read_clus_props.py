"""
created(08/??/11) to make one vec props file
editted(08/22/11) to make one scalar props file

editted(10/11/11) more pc realted scalar properties are added.
at least for scalar property, don't have to save the file but use those functions to read the props I want to use 
get_clus_props.py create the file for properties 
this code contains functions to read: some are changed

triaxiality definition is changed

for the vector properties, the valid lines of sight should be chosen in other codes where these vector properties are used
"""

import os
import numpy as N
import reading_fn as rfn

path_abs = rfn.res_abspath
path_fil = os.path.join(path_abs, "filaments")
path_scat = os.path.join(path_abs, "mass_analysis")
path_data = os.path.join("lcdm250A", "z0.10m3e+10rho5.0")


#### for everything

def get_cluster_ids_mass(filename=os.path.join(path_fil, "finder", "lcdm250A", "z0.10", "halombig_0.10_10.5.halo")):
	clus_info = N.loadtxt(filename, usecols=N.arange(3))
	clus = N.where(clus_info[:, 2] >= 1.e14)[0]
	clus_info = clus_info[clus, :]
	return clus_info

def get_clus_info():
	"""
	use get_cluster_ids_mass() function, save the information already
	because reading entire data takes a lot of time
	"""
	clus_info = N.loadtxt(os.path.join(path_scat, "corr_scatter_props", "clus_info.dat"))
	return clus_info

#### vector properties

def get_subgrp_vec(filename=os.path.join(path_fil, "joanne", "halosubs_11.30_0.9068_13.70.dat")):
	subgrp = N.loadtxt(filename, usecols=(0, 4, 21, 22, 23, 24, 25, 26))
	nclus = subgrp[-1, 0] + 1
	print "nclus in subgrp data: %d " % nclus
	posvel = N.zeros((nclus, 6))
	biggest = N.where(subgrp[:, 1] == 0)[0]
	print "# of biggest subgroup which is %d should be same as the # of clus which is %d" % (len(biggest), nclus) 
	## to check
	clus = N.array(subgrp[biggest, 0], dtype="int")
	if len(N.unique(clus)) != nclus: 
		print "re-write the code"
	else:
		posvel[clus, :3] = subgrp[biggest, 2:5]
		posvel[clus, 3:] = subgrp[biggest, 5:]
	cc = [r"$\hat{r}_{sub}$", r"$\hat{v}_{sub}$"]
	print "got largest subgroup direction from {:s}: {:}".format(filename, posvel.shape)
	return posvel, cc


def get_massivefil_vec(nclus, filename=os.path.join(path_fil, "massivefil_angle2", path_data, "massive_fil_direc_1.00e+14-1.00e+16_10.0_IncE1.dat")):
	fil = N.loadtxt(filename, usecols=(0, 4, 5, 6))
	pos = N.zeros((nclus, 3))
	node = N.array(fil[:, 0], dtype="int")
	pos[node, :] = fil[:, 1:]
	cc = [r"$\hat{r}_{fil}$"]
	print "got massive filament vec from {:s}: {:}".format(filename, pos.shape)
	return pos, cc


def get_normal_vec(filename, nclus):
	normal = N.loadtxt(filename)
	pos = N.zeros((nclus, 3))
	node = N.array(normal[:, 0], dtype="int")
	pos[node, :] = normal[:, 1:]
	if filename.find("bymass") != -1:
		cc = [r"$\hat{n}_{mass}$"]
	elif filename.find("byfil") != -1:
		cc = [r"$\hat{n}_{fil}$"]
	print "got {:s} from {:s}: array shape: {:}".format(cc, filename, pos.shape)
	return pos, cc


def get_clus_vec(filename=os.path.join(path_fil, "cluster_momentinertia", path_data, "cluster_eigenvalvec.dat")):
	clusevec = N.loadtxt(filename, usecols=(4, 5, 6))
	cc = [r"$\hat{l}_1$"]
	print "got cluster long axis from {:s}: {:}".format(filename, clusevec.shape)
	return clusevec, cc
	

def calc_mass_offset_sub(m_perclus, m_true):
	m_ave = N.average(m_perclus, axis=1)
	m_true_arr = N.zeros_like(m_ave)
	m_true_arr.fill(m_true)
	return m_ave - m_true_arr
	

def calc_each_los_mass_offset_from_true_mass(observed, nlos):
	nclus = len(observed)/nlos
	deltaM_arr = N.zeros(len(observed))
	for iclus in xrange(nclus):
		deltaM_arr[iclus*nlos:(iclus+1)*nlos] = calc_mass_offset_sub(observed[iclus*nlos:(iclus+1)*nlos, 2:], observed[iclus*nlos, 1])
	return deltaM_arr
		

def calc_rat_mlike_mtrue(observed, nlos, mtype):
	"""
	mtype: mlh, mlo, ml2, mo2
	returns nlos*nclus array
	"""
	nclus = len(observed)/nlos
	m_ratio = rfn.calc_maxlikelihoodmass(observed, nlos, mtype)
	for iclus in xrange(nclus):
		m_ratio[iclus*nlos:(iclus+1)*nlos] /= observed[iclus*nlos, 1]			
	#m_ratio -= 1.
	return m_ratio


def get_array_for_all_vecs(nodesonly=True):
	"""
	don't save this data any more: to save data, use get_clus_props.py 2011-11-14
	"""
	## get cluster id (by mass, fof) and cluster mass
	#clus_info = get_cluster_ids_mass()
	clus_info = get_clus_info()
	nclus = len(clus_info)
	print "n clus: %d" % nclus 

	ccList = []
	## large subgroup direction
	vec, cc = get_subgrp_vec()
	ccList += cc
	clus_info = N.column_stack((clus_info, vec))

	## massive filament direction
	vec, cc = get_massivefil_vec(nclus)
	ccList += cc
	clus_info = N.column_stack((clus_info, vec))

	## plane direction
	dirname = ["geo_dependence2", "geo_dependence3"]
	filename = ["normal_byfilIncE1_1.00e+14-1.00e+16_t1.5_iem0_iemp0.dat", "normal_bymass_1.00e+14-1.00e+16_lmax10.0_t1.5_iem0_iemp0_nangle100.dat"]
	for dir, file in zip(dirname, filename):
		vec, cc = get_normal_vec(os.path.join(path_fil, dir, path_data, file), nclus)
		ccList += cc
		clus_info = N.column_stack((clus_info, vec))

	## cluster long axis 
	vec, cc = get_clus_vec()
	ccList += cc
	clus_info = N.column_stack((clus_info, vec))

	if nodesonly:
		node = N.loadtxt(os.path.join(path_fil, "massivefil_angle2", path_data, "massive_fil_direc_1.00e+14-1.00e+16_10.0_IncE1.dat"), usecols=[0], dtype='int')
		clus_info = clus_info[node, :]

	return ccList, clus_info
	

##### dot product of vector properties with lines of sight
##### originally in corr_losDotprops_masspc directory
def read_los(dataname = "observed_mass_MGe1.00e+14_fullclus.dat"):
	fullpath = os.path.join(path_fil, "observable", dataname)
	print 'read %s' % fullpath
	los = N.loadtxt(fullpath, usecols=(0, 2, 3, 4))
	## all cluster have same 96 diff los direc  so take clus 0's 96 los direcs
	iuniq = N.nonzero(los[:, 0] == 0)[0] ## nlos
	los = los[iuniq, 1:]
	#los = N.loadtxt(Path + "filaments/observable/" + dataname + '.dat', usecols=(2, 3, 4))
	#los = los[:nlos, :]
	print "nlos = %d" % (len(iuniq))
	return los

def calc_cosang_los_props(los, direc, nclus, takeabs=True):
	nlos = len(los)
	nrow = nlos*nclus
	nprops = direc.shape[1]/3
	print "los shape: {:}".format(los.shape)
	print 'direction shape: {:}'.format(direc.shape)

	for i in xrange(3):
		los[:, i] /= N.sqrt(N.sum(los**2, axis=1)) ## normalize

	cosang = N.empty((nrow, nprops))
	## for each direction property
	for idirec in xrange(nprops):
		a_direc = direc[:, idirec*3:(1+idirec)*3]
		#print "a_direct shape:{:}".format(a_direc.shape)
		norm = N.sqrt(N.sum(a_direc**2, axis=1)) ## nclus x 1 
		for idim in xrange(3):
			a_direc[:, idim] /= norm + 1.e-7## normalize
		## I'll take care of all clusters at once but los by los
		## radier doing cluster by cluster, do los by los
		for ilos in xrange(nlos):
			#print "los: {:}".format(los[ilos, :])
			cosang_ilos = N.abs(N.sum(los[ilos, :]*a_direc, axis=1))
			iwrong = N.where(cosang_ilos > 1.)[0]
			if len(iwrong) > 0:
				print 'iwrong: ', iwrong
				print cosang_ilos[iwrong]
				print norm[iwrong], a_direc[iwrong, :]
				print direc[iwrong, idirec*3:(idirec+1)*3]
				print los[ilos, :]
				i_inerr = N.where(cosang_ilos[iwrong] < 1. + 5.e-4)[0]
				if len(i_inerr) > 0:
					print 'corrected: ', iwrong[i_inerr]
					for ierr in iwrong[i_inerr]:
						if cosang_ilos[ierr] > 1.:
							cosang_ilos[ierr] = 1.
						else:
							cosang_ilos[ierr] = -1.

			cosang[ilos::nlos, idirec] = cosang_ilos

	## take absolute value if takeabs=True
	if takeabs:
		cosang = N.abs(cosang)

	return cosang



#### probability of not detect substructure - los dependent properties
def read_num_of_substruc(filename='substructure22_0.9068_0.txt'):
	fullpath = os.path.join(Path, "filaments/joanne/", filename)
	print 'read %s' % fullpath
	nsubs = N.loadtxt(fullpath)
	nsubs = N.delete(nsubs, N.s_[1:9:], 1)
	return nsubs

def match_clusid(iclus, iclus_org, nsubs):
	igrp = N.where(nsubs[:, 0] == iclus_org)[0]
	return nsubs[igrp, 1]

def get_num_of_substruc(clus_info, nlos):
	"""
	make an array appropriate to stack with cosang arr
	clus_info contains nodes only if nodesonly is true.
	If not, clus_info contains all clusters
	"""
	nsubsArrAll = read_num_of_substruc() ## it is in the order of original id
	nsubsArr = N.zeros(len(clus_info)*nlos)
	## i is required for the calse that iclus is only for nodes
	i = 0
	for iclus, iclus_org in clus_info[:, :2]:
		nsubsArr[i*nlos:(i+1)*nlos] = match_clusid(iclus, iclus_org, nsubsArrAll)
		i += 1
	return nsubsArr

def get_array_los_dot_vprops(nodesonly=True, inc_nsubs=False):
	"""
	to get the array of cosine angle between vector properties and los
	if inc_nsubs = True, include the probability of not detecting substructure along los
	"""
	los = read_los()
	cc_vec, vecprops = get_array_for_all_vecs(nodesonly=nodesonly)
	print "read vecprops: {:}".format(cc_vec)

	nclus = vecprops.shape[0] ## could be the number of nodes if nodesonly==True
	nlos = len(los)

	# calculate the cosine angle
	vpTOsp = calc_cosang_los_props(los, vecprops[:, 3:], nclus)

	if inc_nsubs:
		cc_vec += r"$n_{sub}$"
		nsubsArr = get_num_of_substruc(vecprops[:, :2], nlos)
		vpTOsp = N.column_stack((vpTOsp, nsubsArr))

	return vpTOsp, cc_vec, N.array(vecprops[:, 0], dtype="int")


########## scalar properties

def diff_periodic(arr):
	"""
	before applying the box size
	"""
	for icol in xrange(arr.shape[1]):
		ii = N.array([0])
		while len(ii) != 0:
			ii = N.where(arr[:, icol] < -0.5)[0]
			arr[ii, icol] += 1.
		ii = N.array([0])
		while len(ii) != 0:
			ii = N.where(arr[:, icol] >= 0.5)[0]
			arr[ii, icol] -= 1.
	return arr


def calc_env_mass(iclus, halo_info, inc_center, boxsize=250., rad=10.):
	hshift = halo_info[:, 3:] - halo_info[iclus, 3:]
	## periodic condition
	hshift = diff_periodic(hshift)
	dst = N.sqrt(N.sum(hshift**2, axis=1))
	inear = N.where(dst<=rad/boxsize)[0]
	print "# of halos near the cluster: %d" % len(inear)
	tot_mass = N.sum(halo_info[inear, 2])
	if inc_center == False:
		tot_mass -= halo_info[iclus, 2]
	return tot_mass
	

def get_env_mass(filename=os.path.join(path_fil, "finder", "lcdm250A", "z0.10", "halombig_0.10_10.5.halo"), mass_limit=5.e13, inc_center=False):
	"""
	calculate the total mass of halos (>mass limit)
	"""
	halo_info = N.loadtxt(filename, usecols=N.arange(6))
	print "read %s" % filename
	grpnum = N.where(halo_info[:, 2] >= mass_limit)[0]
	halo_info = halo_info[grpnum, :]
	print "keep halos more massive than %.4e: %d halos" % (mass_limit, len(grpnum))
	clusnum = N.where(halo_info[:, 2] >= 1.e14)[0]
	saveArr = N.empty((len(clusnum), 2))
	for ii, iclus in enumerate(clusnum):
		tot_mass = calc_env_mass(iclus, halo_info, inc_center)
		saveArr[ii, 0] = iclus
		saveArr[ii, 1] = tot_mass
	N.savetxt("env_mass_halos_ge%.2e_incClus%s.dat" % (mass_limit, inc_center), saveArr, fmt="%8d\t%.6e")
	return saveArr[:, 1], [r"$M_{sphere}$"] #[r"$M_{nearhalo(>%.0e)}$" % mass_limit]

def read_env_mass(mass_limit=5.e13, incClus=False):
	"""
	once get_env_mass is done, I can just use get_env_mass function
	"""
	env_mass  = N.loadtxt(os.path.join(path_scat, "corr_scatter_props", "pca_all_scalprops", "env_mass_halos_ge%.2e_incClus%s.dat" % (mass_limit, incClus)), usecols=[1])
	expo = N.int(N.log10(mass_limit))
	fac = mass_limit/(10.**expo)
	return env_mass, [r"$M_{sphere}$"] #, [r"$\Sigma M_{h(>%.0f\times10^{%d})}$" % (fac, expo)]


def calc_ave_pc_i(pc_arr):
	"""
	calcluate the average
	copied from corr_scatter/get_stats_pc/
	pc_arr: (nclus, nobs) arr
	"""
	pc_ave = N.sum(pc_arr, axis=0)/len(pc_arr)
	print 'average: {:}'.format(pc_ave)

	#pc_ave /= N.sqrt(N.dot(pc_ave, pc_ave))

	return pc_ave


def calc_most_likely_pc_i(pc_arr):
	"""
	calculate the most likely pc_i (see the derivation in the note)
	pc_arr is pc_i for all clusters, nclus x nobs array
	copied from corr_scatter/get_stats_pc/
	pc_arr: (nclus, nobs) arr
	"""
	narr = pc_arr.shape[1]
	arr = N.zeros((narr, narr))
	nclus = len(pc_arr)
	for iclus in xrange(nclus):
		arr += N.outer(pc_arr[iclus, :], pc_arr[iclus, :])

	l_multiplier, min_vec = N.linalg.eig(arr)
	isort_l = N.argsort((nclus - l_multiplier))
	pc_likely = min_vec[:, isort_l[0]]

	print 'most likely pc: {:}'.format(pc_likely)

	return pc_likely
	
def calc_dotprod(pc_arr, f_dotprod):
	"""
	pc_arr contains all clusters' pc_i, dimension is pc_arr = nclus x nobs
	f_dotprod is function name: calc_most_likely_pc_i or calc_ave_pc_i
	"""
	pc_rep = f_dotprod(pc_arr)
	dotprod_arr = N.abs(N.sum(pc_rep*pc_arr, axis=1))
	return dotprod_arr
	

def get_eval_pc0(rfilename, nobs):
	"""
	get the sum of absolute eigenvalues and the fractional eigenvalue for pc0
	"""
	fname = os.path.join(path_scat, "corr_scatter", "get_pc", "eigen", "eigenvalvecs_%s.dat" % rfilename)
	evals = N.loadtxt(fname, usecols=[1])
	print "read: %s to get sum of the eigenvalues and fractional eigenvalue of pc0" % fname

	evals = evals.reshape(len(evals)/nobs, nobs)

	arr_all = N.empty((len(evals), 5))

	arr_all[:, 0] = N.sum(evals, axis=1) ## sum of all eigenvalues
	arr_all[:, 1] = N.prod(evals, axis=1) ## product of all eigenvalues
	arr_all[:, 2] = evals[:, 0]/arr_all[:, 0]  ## pc0/pctot
	#arr_all[:, 3] = evals[:, 1]/evals[:, 0] ## pc1/pc0
	arr_all[:, 3] = evals[:, 1]
	arr_all[:, 4] = evals[:, -1]

	#cc = ["totVar", "varVol", "Epc1", "PC2toPC1"]
	#cc = ["totVar", "varVol", "Epc1", "Epc2"]
	#cc = [r"$\Sigma_i PC_{i,m}$", r"$\Pi_i PC_{i, m}$", r"$\frac{PC_{0, m}}{\Sigma_i PC_{i, m}}$", r"$PC_{1, m}$", r"$PC_{4, m}$"]
	#cc = [r"$\Sigma_i PC_{i,m}$", r"$\Pi_i PC_{i, m}$", r"${PC_{0, m}}/{\Sigma_i PC_{i, m}}$", r"$PC_{1, m}$", r"$PC_{4, m}$"]
	cc = [r"$\Sigma \lambda$", r"$\Pi \lambda$", r"${\lambda_{0, M}}/\Sigma {\lambda}$", r"$\lambda_{1, M}$", r"$\lambda_{4, M}$"]
	#cc = [r"$\sum_i PC_{i,mass}$", r"$\prod_i PC_{i, mass}$", r"$\frac{PC_{0, mass}}{\small\sum_i PC_{i, mass}}$", r"$PC_{1, mass}$", r"$PC_{4, mass}$"]

	return arr_all, cc


def get_pc0_dot_pcrep(rfilename, nobs, f_rep = calc_most_likely_pc_i):
	fname = os.path.join(path_scat, "corr_scatter", "get_pc", "eigen", "eigenvalvecs_%s.dat" % rfilename)
	evecs = N.loadtxt(fname)
	evecs = evecs[:, 2:]
	pc0 = evecs[:, 0]
	pc0 = pc0.reshape(len(pc0)/nobs, nobs)
	dotprod_arr = calc_dotprod(pc0, f_rep)
	#cc = [r"${\hat PC_{0,m}}\cdot{\hat PC_{0, m, all}}$"]
	cc = [r"$|\cos\theta_{0, minsq}|$"]	
	return dotprod_arr, cc


def calc_ave_mass_offset(observed, nlos):
	nclus = len(observed)/nlos
	nobs = observed[:, 2:].shape[1]
	m_off = N.zeros((nclus, nobs))
	for iclus in xrange(nclus):
		m_off[iclus, :] = N.average(observed[iclus*nlos:(iclus+1)*nlos, 2:], axis=0)
		m_off[iclus, :] -= observed[iclus*nlos, 1]
	return m_off
		

def get_ang_btw_mass_offset_and_pc(observed, pcfilename, ipc, nlos, normalize=True):
	"""
	dot product between (<M_obs^alpha>_i - M_true) and pc vector
	"""
	nclus = len(observed)/nlos
	nobs = observed[:, 2:].shape[1]

	evals, evecs = rfn.read_principal_component(pcfilename)
	if pcfilename.find("alllos") != -1:
		pc_i = N.zeros((nclus, nobs))
		for i in xrange(nobs):
			pc_i[:, i].fill(evecs[i, ipc])
	else:
		pc_i = evecs[:, ipc].reshape(len(evecs)/nobs, nobs)
	print "to get the angle between the offset, read %s, and used pc%d" % (pcfilename, ipc)
		
	m_off = calc_ave_mass_offset(observed, nlos)
	print m_off.shape
	if normalize:
		m_off_norm = N.sqrt(N.sum(m_off*m_off, axis=1))
		for m_off_each in m_off.T:
			m_off_each /= m_off_norm 
	abscosang = N.abs(N.sum(m_off*pc_i, axis=1))
	
	if normalize:
		print "check cosang > 1"
		ii = N.where(abscosang > 1.)[0]
		print abscosang[ii]
	#cc = [r"${\Delta \vec{M}}\cdot{\hat PC_{%d, m}}$" % ipc]
	cc = [r"$\cos\theta_{\Delta, 0}$"]
	return abscosang, cc


def calc_ave_mass(mest_each):
	"""
	mest_each = nlos x 1 array
	"""
	ivalid = N.where(mest_each > 1.e5)[0]
	mave = N.average(mest_each[ivalid])
	return mave


def get_scat_size_ave(observed, obsname, nobs, nlos):
	"""
	<Mobs> - Mtrue
	"""
	print "get_scat_size_ave"
	print observed.shape
	nclus = len(observed)/nlos
	scatAve_arr = N.empty((nclus, nobs))
	for iclus in xrange(nclus):
		scatter, ikeep = rfn.get_valid_array_altogether(observed[iclus*nlos:(iclus+1)*nlos, 2:], observed[iclus*nlos, 1], "sca")
		scatAve_arr[iclus, :] = N.average(scatter, axis=0)
	cc = []
	for obsname1 in obsname:
		cc.append("$\Delta %s$" % obsname1[1:-1])	
	#cc = [r"$\Delta M_{N_{red}}/M_{true}$", r"$\Delta M_{N_{ph}}/M_{true}$", r"$\Delta M_{SZ}/M_{true}$", r"$\Delta M_{vel}/M_{true}$", r"$\Delta M_{WL}/M_{true}$"]
	return scatAve_arr, cc


def get_scat_size_ave_old(nobs, nlos):
	"""
	<Mobs> - Mtrue
	"""
	observed, obsname = rfn.read_data()
	print "get_scat_size_ave"
	print observed.shape
	print observed[0, :]
	nclus = len(observed)/nlos
	scatAve_arr = N.empty((nclus, nobs))
	for iclus in xrange(nclus):
		scatter, ikeep = rfn.get_valid_array_altogether(observed[iclus*nlos:(iclus+1)*nlos, 2:], observed[iclus*nlos, 1], "sca")
		scatAve_arr[iclus, :] = N.average(scatter, axis=0)
		#for iobs in xrange(nobs):
			#scat_arr[iclus, iobs] = (observed[iclus*nlos, 1] - calc_ave_mass(observed[iclus*nlos:(iclus+1)*nlos, 2+iobs]))/observed[iclus*nlos, 1]
			#scat_arr[iclus, iobs] = (calc_ave_mass(observed[iclus*nlos:(iclus+1)*nlos, 2+iobs]) - observed[iclus*nlos, 1])/observed[iclus*nlos, 1]
	cc = []
	for obsname1 in obsname:
		cc.append("$\Delta %s$" % obsname1[1:-1])	
	#cc = [r"$\Delta M_{N_{red}}$", r"$\Delta M_{N_{ph}}$", r"$\Delta M_{SZ}$", r"$\Delta M_{vel}$", r"$\Delta M_{WL}}$"]
	return scatAve_arr, cc


def get_clus_shape_params():
	"""
	get triaxiality and sphericity
	Hahn et al 2007 & Bett et al 2007
	note that the order of eigenvalues I saved is I3>I2>I1
	the order of corresponding axis length is opposite l1>l2>l3
	l1 = sqrt(5/2M_h(-I1+I2+I3))
	l2 = sqrt(5/2M_h(-I2+I3+I1))
	l3 = sqrt(5/2M_h(-I3+I1+I2))
	tri = (l1^2 - l2^2) / (l1^2 - l3^2), sph = l3/l1 (l1>l2>l3)
	"""
	fname = os.path.join(path_fil, "cluster_momentinertia", path_data, "cluster_eigenvalvec.dat")
	evals = N.loadtxt(fname, usecols=N.arange(1, 4))
	print "read cluster moment of inertia eigenvalues: {:}".format(evals.shape)
	l1sq = evals[:, 1] + evals[:, 2] - evals[:, 0]
	l2sq = evals[:, 2] + evals[:, 0] - evals[:, 1]
	l3sq = evals[:, 0] + evals[:, 1] - evals[:, 2]
	shapes = N.empty((len(evals), 2))

	## triaxiality
	shapes[:, 0] = (l1sq - l2sq)/(l1sq - l3sq)
	## sphericity
	shapes[:, 1] = N.sqrt(l3sq/l1sq)

	#cc = ["tri", "sph"]
	cc = [r"$T$", r"$S$"]

	return shapes, cc
	

def get_masses_in_planes(iplane, nclus):
	"""
	get the all halo mass or all filament mass in different definitions of planes
	"""
	dir = "geo_dependence%d" % iplane
	mass = N.zeros(nclus)

	if iplane == 2:
		print "connected filament plane - all filament mass in 10Mph/h sphere"
		fname = "filmassplane_byfilInc1_1.00e+14-1.00e+16_lmax20.0_dl10.0_t1.5_iem0_iemp0.dat"
		#cc = ["CFP-CFM"]
		cc = [r"${M_{fplane}}$"]

		icol = 3
	elif iplane  == 3:
		print "halo mass plane - all halo mass in plane in 10Mpc/h sphere"
		fname = "massplane_1.00e+14-1.00e+16_lmax10.0_t1.5_iem0_iemp0_nangle100.dat"
		#cc = ["AMP-AHM"]
		cc = [r"${M_{hplane}}$"]
		icol = 3
	elif iplane  == 4:
		print "all filament plane - filament mass in plane in 10Mpc/h sphere"
		fname = "filmassplane_byallfilInc1_1.00e+14-1.00e+16_lmax20.0_dl10.0_t1.5_iem0_iemp0.dat"
		#cc = ["AFP-AFM"]
		cc = [r"${M_{allfilplane}}$"]
		icol = 3
	elif iplane  == 5:
		print "galaxy richness plane - galaxy richness in plane in 10Mpc/h sphere"
		fname = "richnessplane_1.00e+14-1.00e+16_lmax10.0_t1.5_nangle60.dat"
		#cc = ["GRP-AGR"]
		cc = [r"${M_{richplane}}$"]
		icol = 3
	else:
		print "iplane should be [2, 5]"
		raise

	fullname = os.path.join(path_fil, dir, path_data, fname)
	data = N.loadtxt(fullname, usecols=(0, icol))
	inode = N.array(data[:, 0], dtype='int')
	mass[inode] = data[:, 1]
	print "shape: mass={:}".format(mass.shape)

	return mass, cc


def get_frac_in_planes(iplane, nclus):
	"""
	get the *fraction* of all halo mass or all filament mass in different definitions of planes to the mass in the sphere
	"""
	dir = "geo_dependence%d" % iplane
	mass = N.zeros(nclus)

	if iplane == 2:
		print "connected filament plane - all filament mass in 10Mph/h sphere"
		fname = "filmassplane_byfilInc1_1.00e+14-1.00e+16_lmax20.0_dl10.0_t1.5_iem0_iemp0.dat"
		#cc = ["FrCFP-CFM"]
		cc = [r"$f_{M_{fplane}}$"]
		icol = 3
	elif iplane  == 3:
		print "halo mass plane - all halo mass in plane in 10Mpc/h sphere"
		fname = "massplane_1.00e+14-1.00e+16_lmax10.0_t1.5_iem0_iemp0_nangle100.dat"
		#cc = ["FrAMP-AHM"]
		cc = [r"$f_{M_{hplane}}$"]
		icol = 3
	elif iplane  == 4:
		print "all filament plane - filament mass in plane in 10Mpc/h sphere"
		fname = "filmassplane_byallfilInc1_1.00e+14-1.00e+16_lmax20.0_dl10.0_t1.5_iem0_iemp0.dat"
		#cc = ["FrAFP-AFM"]
		cc = [r"$f_{M_{allfilplane}}$"]
		icol = 3
	elif iplane  == 5:
		print "galaxy richness plane - galaxy richness in plane in 10Mpc/h sphere"
		fname = "richnessplane_1.00e+14-1.00e+16_lmax10.0_t1.5_nangle60.dat"
		#cc = ["FrGRP-AGR"]
		cc = [r"$f_{M_{richplane}}$"]
		icol = 3
	else:
		print "iplane should be [2, 5]"
		raise

	fullname = os.path.join(path_fil, dir, path_data, fname)
	data = N.loadtxt(fullname, usecols=(0, icol-1, icol))
	inode = N.array(data[:, 0], dtype='int')
	mass[inode] = data[:, 2]/data[:, 1]
	print "shape: mass={:}".format(mass.shape)

	return mass, cc

	
def get_largest_subgrp_props():
	"""
	calculate the fraction of richness in largest subgroup and the fraction of the distance to largest subgroup form the center to the length of the longest axis of the cluster
	"""
	fname = os.path.join(path_fil, "joanne", "halosubs_11.30_0.9068_13.70.dat")
	halosubs = N.loadtxt(fname)
	print "read %s to get largest subgroup information" % fname

	## 4th column is index of subgroup (lowest index is the biggest)
	ilg_subgrp = N.where(halosubs[:, 4] == 0)[0]
	## n of galaxies in cluster(3rd column) will be the same as n of galaxies in subgroup (6th column) if that row contains the information of the entire halo
	iall_halo = N.where(halosubs[:, 3] == halosubs[:, 6])[0] 

	subgrp_props = N.zeros((len(ilg_subgrp), 2))
	## fraction of richness
	subgrp_props[:, 0] = 1.*halosubs[ilg_subgrp, 6]/halosubs[ilg_subgrp, 3]
	## distance to the subgroup: halosubs[:, 16] is the longest axis of the subgroup
	subgrp_props[:, 1] = N.sqrt(N.sum(halosubs[ilg_subgrp, 21:24]**2, axis=1))/halosubs[iall_halo, 16]

	#cc = ["fR-sb", "fD-sb"]
	cc = [r"$f_{R_{sub}}$", r"$f_{D_{sub}}$"]

	return subgrp_props, cc


def calc_cm(a_clusinfo, gals):
	"""
	a_clusinfo: one cluster info: id by mass, haloid, clus mass, position, velocity, r180, cvir, merger times
	gals: x, y, z, vx, vy, vz, haloid, log infall mass, satellite(1)/center(0), infall time
	calculate center of mass after subtracting center galaxy position to apply periodic condition easily - so later, i don't have to subtract cetner galaxy position
	03/28/2012: use just position: don't use infall mass
	"""
	imems = N.where(gals[:, 6] == a_clusinfo[1])[0]
	##
	shifted_galpos = gals[imems, :3] - a_clusinfo[3:6]
	shifted_galpos = diff_periodic(shifted_galpos)
	cm = N.zeros(shifted_galpos.shape[1])

	for i in xrange(len(cm)):
		cm[i] = N.sum(shifted_galpos[:, i])/len(imems)

	return cm


def calc_cm_old(a_clusinfo, gals, use_logmass):
	"""
	a_clusinfo: one cluster info: id by mass, haloid, clus mass, position, velocity, r180, cvir, merger times
	gals: x, y, z, vx, vy, vz, haloid, log infall mass, satellite(1)/center(0), infall time
	calculate center of mass after subtracting center galaxy position to apply periodic condition easily - so later, i don't have to subtract cetner galaxy position
	"""
	imems = N.where(gals[:, 6] == a_clusinfo[1])[0]
	##
	shifted_galpos = gals[imems, :3] - a_clusinfo[3:6]
	shifted_galpos = diff_periodic(shifted_galpos)
	cm = N.zeros(shifted_galpos.shape[1])
	if use_logmass:
		mass = gals[imems, 7]
		mgaltot = N.sum(mass)
	else:
		mass = 10.**gals[imems, 7]
		mgaltot = N.sum(mass)

	for i in xrange(len(cm)):
		#cm[i] = N.sum(mass*shifted_galpos[:, i])/a_clusinfo[2]
		cm[i] = N.sum(mass*shifted_galpos[:, i])/mgaltot

	return cm


def get_frac_dst_center_to_cm(boxsz=250.):
	fname = os.path.join(path_fil, "joanne", "subh_0.9068.gal")
	gals = N.loadtxt(fname) ## galaxy file
	fname = os.path.join(path_fil, "joanne", "halonumvals_0.9068_m13.70.dat")
	clusinfo = N.loadtxt(fname) ## poitions are the same as the positions fo the center galaxies
	cmArr = N.zeros((len(clusinfo), 3))
	for iclus in xrange(len(clusinfo)):
		print "clus% d" % iclus
		cmArr[iclus, :] = calc_cm(clusinfo[iclus, :], gals)

	## no need division - commented out
	"""
	if use_r180:
		rr = clusinfo[:, 9]
	else:
		rr = clusinfo[:, 2]**(1./3.)
	"""
	print cmArr[:4, :]*boxsz
	fracdst = boxsz*N.sqrt(N.sum(cmArr*cmArr, axis=1))
	#print fracdst[:4] 
	#fracdst /= rr
		
	ff = open(os.path.join(path_scat, "corr_scatter_props", "frac_dst_center_to_cm.dat") , "w")
	ff.write("#number weighted center of mass\n")
	ff.write("#center of mass setting the center galaxy (0, 0, 0)\tfractional distance between the center galaxy and center of mass\n")
	N.savetxt(ff, N.column_stack((cmArr*boxsz, fracdst)))
	ff.close()
	return fracdst, [r"$x_{off}$"] #, [r"$f_{cg-cm}$"]


def read_frac_dst_center_to_cm():
	fracdst = N.loadtxt(os.path.join(path_scat, "corr_scatter_props", "frac_dst_center_to_cm.dat"), usecols=[3])
	return fracdst, [r"$x_{off}$"] #[r"$f_{cg-cm}$"]


def get_concentration():
	fname = os.path.join(path_fil, "joanne", "mass_cvir.dat") 
	clus_cvir = N.loadtxt(fname, usecols=[2])
	print "concentration: read %s" % fname
	print "cvir shape: {:}".format(clus_cvir.shape)

	cc = [r"$c_{vir}$"]

	return clus_cvir, cc


def get_recent_merger():
	"""
	get the time of most recent major merger 1:3 and 1:10
	"""
	fname = os.path.join(path_fil, "joanne", "halonumvals_0.9068_m13.70.dat")
	times = N.loadtxt(fname, usecols=(11, 12))
	print "the time of most recent merger: read %s" % fname
	print "shape: {:}".format(times.shape)

	#cc = ["t03", "t10"]
	cc = [r"$t_{1:3}$", r"$t_{1:10}$"]

	return times, cc


def get_cor_pc_mlikeTOmtrue(pcfilename, mtype_like="mlh", normalize=True, ipc=0):
	path = os.path.join(path_scat, "corr_scatter_props", "corr_mlike_masspc", "data")
	filename = "btw_mlike-%sTOmtrue_pc-%s_norm%s" % (mtype_like, pcfilename, normalize)
	cor = N.loadtxt(os.path.join(path, "pea_%s.dat" % filename))
	#cov = N.loadtxt(os.path.join(path, "cov_%s.dat" % filename))
	#cc = [r"$\langle\frac{M_{like}}{M_{true}}, \frac{a_{%d}^\alpha}{\sqrt{\sum_i(a_i^\alpha)^2}} \rangle$" % ipc]
	cc = [r"$\langle\frac{M_{like}}{M_{true}}, \cos\theta_{%d} \rangle$" % ipc]
	#cc = [r"$\langle \frac{M_{\rm like}}{M_{\rm true}} PC_{%d} \rangle$" % ipc, r"$\langle \frac{M_{\rm like}}{M_{\rm true}} PC_{%d} \rangle_{\rm cov}$" % ipc] 
	#return N.column_stack((cor[:, 1+ipc], cov[:, 1+ipc])), cc
	return cor[:, 1+ipc], cc


def get_node_clus(arr):
	inode = N.where(N.any(arr !=0., axis=1))[0]
	print "only nodes will be used: node=%d" % len(inode)
	return inode


def get_array_for_all_scals(i_rfile, nodesonly=False, mass_limit_env=5.e13, take_abs_scat=False, mtype_like="mo2"):
	"""
	originally from pca_all_scalprops/get_corr_allscals.py
	"""
	nobs = 5
	nlos = 96
	nbox = 10
	filenameArr = rfn.get_filename(nbox)
	rfilename = filenameArr[i_rfile]
	cc = []

	## cluster info
	clus_info = get_clus_info()
	allArr = clus_info.copy()
	observed, obsname = rfn.read_data()

	## scatter size of observables (<Mobs>-Mtrue)/Mtrue
	scatArr, cc_scat = get_scat_size_ave(observed, obsname, nobs, nlos)
	print "scatArr: {:}".format(cc_scat)
	for i in xrange(5):
		print N.min(scatArr[:, i]), N.max(scatArr[:, i])
	if len(cc_scat) != nobs:
		print "scatter array size is wrong"
		raise

	if take_abs_scat:
		scatArr = N.abs(scatArr)

	nclus =  len(scatArr)
	cc += cc_scat
	allArr = N.column_stack((allArr, scatArr))

	## dot prodcut between pc0 and most_likely_pc
	dotprodArr, cc_dotprod = get_pc0_dot_pcrep(rfilename, nobs)
	print "dotprodArr: {:}".format(cc_dotprod)
	print N.min(dotprodArr), N.max(dotprodArr)
	cc += cc_dotprod
	allArr = N.column_stack((allArr, dotprodArr))

	evalprops, cc_eval = get_eval_pc0(rfilename, nobs)
	print "evalArr: {:}".format(cc_eval)
	nevalprops = len(cc_eval)
	cc += cc_eval
	allArr = N.column_stack((allArr, evalprops))

	## cluster shape parameter
	shapeArr, cc_shapes = get_clus_shape_params()
	print "shape: {:}".format(cc_shapes)
	cc += cc_shapes
	allArr = N.column_stack((allArr, shapeArr))

	## get fractional mass in planes
	fm_con, cc_fmcon = get_frac_in_planes(2, nclus)
	fm_ahl, cc_fmahl = get_frac_in_planes(3, nclus)
	print "fractional mass in plane: {:}, {:}".format(cc_fmcon, cc_fmahl)
	print N.min(fm_con), N.max(fm_con)
	print N.min(fm_ahl), N.max(fm_ahl)
	cc += cc_fmcon + cc_fmahl
	allArr = N.column_stack((allArr, fm_con, fm_ahl))

	## environment mass
	envmassArr, cc_env = read_env_mass(mass_limit=mass_limit_env)
	#envmassArr, cc_env = get_env_mass(mass_limit=mass_limit_env)
	print "envmass: {:}".format(cc_env)
	cc += cc_env
	allArr = N.column_stack((allArr, envmassArr))

	## true mass
	mtrue = observed[::nlos, 1]
	cc += [r"$M_{true}$"]
	allArr = N.column_stack((allArr, mtrue))

	## largest subgrp props
	subgrp_props, cc_subgrp = get_largest_subgrp_props()
	print "subgrp props:{:}".format(cc_subgrp)
	cc += cc_subgrp
	allArr = N.column_stack((allArr, subgrp_props))

	## concentraion
	clus_cvir, cc_cvir = get_concentration()
	print "concentration: {:}".format(cc_cvir)
	cc += cc_cvir
	allArr = N.column_stack((allArr, clus_cvir))

	## merger time
	times, cc_times = get_recent_merger()
	print "merger time: {:}".format(cc_times)
	cc += cc_times
	allArr = N.column_stack((allArr, times))

	ncol = len(cc)

	## frac dist
	#fdsttocm, cc_fdsttocm = get_frac_dst_center_to_cm()
	fdsttocm, cc_fdsttocm = read_frac_dst_center_to_cm()
	print "frac dist: {:}".format(fdsttocm.shape)
	cc += cc_fdsttocm
	allArr = N.column_stack((allArr, fdsttocm))

	## cor or cov between projected value on PC0 and mtrue/mlike
	#corcov_mrat, cc_corcov_mrat = get_cor_pc_mlikeTOmtrue(rfilename, mtype_like=mtype_like, normalize=True, ipc=0)
	#print "cor or cov of mlike/mtrue with projected on PC0: {:}".format(corcov_mrat.shape)
	#cc += cc_corcov_mrat
	#allArr = N.column_stack((allArr, corcov_mrat)) 

	## ang between Delta M and PC0
	ang_DM_PC, cc_ang_DM_PC = get_ang_btw_mass_offset_and_pc(observed, rfilename, 0, nlos, normalize=True)
	print "ang between Delta M and PC_0: {:}".format(ang_DM_PC.shape)
	cc += cc_ang_DM_PC
	allArr = N.column_stack((allArr, ang_DM_PC))

	#########

	if nodesonly:
		inode = get_node_clus(N.column_stack((fm_con, fm_ahl)))
		allArr = allArr[inode, :]
	print "array length: %d" % len(allArr)

	return cc, allArr
