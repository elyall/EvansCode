import os
import sys
import numpy as np
import h5py
from scipy.io import loadmat

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn import metrics

from dask.distributed import Client, LocalCluster

from matplotlib import pyplot as plt


def main(filebase, string='_catch', n_X=15, n_iter=500, n_splits=10, max_iter=200, save_out=True):
	
	# Gather data
	data, StimLog, StimID, fn, TrialIndex, NeuronIndex = load(filebase, string) # load in experiment & neural data
	_, N = np.shape(data)
	number_of_neurons = determine_num(N, n_X) # determine # of neurons to subsample at log spacing
	classifiers = return_classifiers(max_iter) # load in classifiers to run

	# Run classifiers
	pred, perc_correct, number_of_neurons, _, classifiers = run(data, 
		StimLog[StimID,:], StimID, number_of_neurons, classifiers=classifiers, 
		n_iter=n_iter, n_splits=n_splits)
	
	# Save output
	save_file = ''
	if save_out:
		save_file = filebase+'.h5'
		save_file = save(pred, perc_correct, number_of_neurons, StimLog, StimID, classifiers, 
			TrialIndex, NeuronIndex, save_file=filebase+'.h5')

	return pred, perc_correct, number_of_neurons, StimLog, StimID, classifiers, save_file


def load(filebase, string):

	# Determine filenames
	expfile = filebase + '.exp'
	if os.path.isfile(filebase + string + '.rois'):
		fn = [filebase + string + '.rois']
	else:
		fn = [filebase + "_depth%01d" % d + string + '.rois' for d in [1,2,3,4]]

	# Load stim info
	h5f = h5py.File(expfile,'r')
	StimID = h5f['/TrialInfo/StimID'][:].astype('int')
	StimLog = h5f['/Experiment/stim/stim'][:].transpose().astype('bool')
	TrialIndex = h5f['TrialIndex'][:].astype('int') - 1
	try:
		NeuronIndex = h5f['NeuronIndex'][:].astype('int') - 1
	except:
		NeuronIndex = []
	h5f.close()

	# Load neural data
	temp = []
	for f in fn:
	    h5f = h5py.File(f,'r')
	    numNeurons = len(h5f['ROIdata/rois/stimMean'])
	    numTrials = len(h5f['ROIdata/DataInfo/StimID'][0])
	    data = np.zeros([numTrials,numNeurons])
	    for n in np.arange(numNeurons):
	        data[:,n] = h5f[h5f['ROIdata/rois/stimMean'][n][0]][:]
	    h5f.close()
	    temp.append(data)
	data = np.concatenate(temp[:],axis=1)

	# keep only desired trials & neurons
	StimID = np.squeeze(StimID[TrialIndex])
	data = np.squeeze(data[TrialIndex,:])
	if NeuronIndex:
		data = np.squeeze(data[:,NeuronIndex])
	# data = data.nan_2_num(data) # convert nan's to 0's
	data = data[:,~np.any(np.isnan(data), axis=0)] # remove neurons that have at least 1 nan

	# # Convert StimLog to vector
	# a = [1,2,4,8,16]
	# StimLog = StimLog.astype('int')
	# for i in np.arange(StimLog.shape[1]):
	# #     StimLog[:, i] *= i+1
	#     StimLog[:, i] *= a[i]
	# StimLog = StimLog.sum(axis=1)

	return data, StimLog, StimID, fn, TrialIndex, NeuronIndex


def save(pred, perc_correct, number_of_neurons, StimLog, StimID, classifiers, TrialIndex, 
	NeuronIndex, save_file='temp.h5'):

	# Write output to hdf5 file
	h5f = h5py.File(save_file, 'w')
	h5f.create_dataset('/perc_correct', data=perc_correct)
	h5f.create_dataset('/pred', data=pred)
	h5f.create_dataset('/number_of_neurons', data=number_of_neurons)
	h5f.create_dataset('/StimLog', data=StimLog)
	h5f.create_dataset('/StimID', data=StimID)
	asciiList = [n.encode("ascii", "ignore") for n in classifiers]
	h5f.create_dataset('/classifiers', (len(asciiList),1), 'S10', asciiList)
	h5f.create_dataset('/TrialIndex', data=TrialIndex)
	h5f.create_dataset('/NeuronIndex', data=NeuronIndex)
	h5f.close()

	return save_file


def determine_num(num_neurons, N=15):
	# n = np.logspace(0,np.log10(num_neurons),N).round().astype('int') # list of number of neurons to sample
	n = np.logspace(0,np.log10(100),N).round().astype('int') # list of number of neurons to sample
	n = n[n<num_neurons]
	n = np.append(n,num_neurons)
	return np.unique(n)


def return_classifiers(max_iter=100):
	from sklearn.linear_model import LogisticRegression, RidgeClassifierCV
	from sklearn.svm import SVC
	from sklearn.multiclass import OneVsRestClassifier
	from sklearn.tree import DecisionTreeClassifier, ExtraTreeClassifier
	from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
	from sklearn.neighbors import KNeighborsClassifier, RadiusNeighborsClassifier
	from sklearn.neural_network import MLPClassifier

	classifiers = {
    # # 'Logistic Regression' : LogisticRegression(multi_class='multinomial',max_iter=max_iter,solver='sag'),
    #'Ridge (CV)': RidgeClassifierCV(), # has formatting issue (pred is Nx1 not Nx5)
    # # 'One vs Rest' : OneVsRestClassifier(SVC(kernel='linear')),
    #'Decision Tree' : DecisionTreeClassifier(),
    #'Extra Tree' : ExtraTreeClassifier(),
    #'Extra Trees' : ExtraTreesClassifier(),
    #'Random Forest': RandomForestClassifier(),
    #'K Neighbors' : KNeighborsClassifier(),
    # 'Radius Neighbors': RadiusNeighborsClassifier(radius=25.0), # very sensitive to radius
    'Multi-layer Perceptron' : MLPClassifier(max_iter=max_iter),
	}
	return classifiers


def run(data, StimLog, StimID, number_of_neurons, classifiers=return_classifiers(),
	n_iter=100, n_splits=10):

	cluster = LocalCluster(n_workers=20)
	e = Client(cluster)

	def classify(n_neurons):

	    n_classifiers = len(classifiers)
	    numTrials, numNeurons = np.shape(data)

	    # initialize outputs
	    # pred = np.zeros([numTrials,n_classifiers,n_iter])
	    pred = np.zeros([numTrials,5,n_classifiers,n_iter])
	    perc_correct = np.zeros([n_classifiers,n_iter])
	    
	    for ind in np.arange(n_iter):

	        # pull out sample of neurons
	        my_sample = np.random.choice(range(numNeurons), n_neurons, replace=False)
	        current = data[:,my_sample]

	        # perform K-fold cross validation (train and test on different pieces of data)
	        skf = StratifiedShuffleSplit(n_splits=n_splits)
	        for train_index, test_index in skf.split(current, StimID): # does skf.split() stratify multilabel data properly?
	            X_train, X_test = current[train_index,:], current[test_index,:]
	            # y_train, y_test = StimLog[train_index], StimLog[test_index]
	            y_train, y_test = StimLog[train_index,:], StimLog[test_index,:]
	                        
	            # fit each classifer
	            for index, (name, classifier) in enumerate(classifiers.items()):
	                classifier.fit(X_train, y_train)
	                # pred[test_index,index,ind] = classifier.predict(X_test)
	                pred[test_index,:,index,ind] = classifier.predict(X_test)
	        
	        # compute performance
	        for index in range(len(classifiers)):
	            # perc_correct[index,ind] = metrics.accuracy_score(StimLog, pred[:,index,ind])
	            perc_correct[index,ind] = metrics.accuracy_score(StimLog, pred[:,:,index,ind])

	    return perc_correct, pred

	# args = ((data, StimLog, StimLog, N, n_iter, n_splits, classifiers) for N in number_of_neurons)
	# Futures = e.map(lambda p: classify(*p), args)
	Futures = e.map(classify, number_of_neurons)
	out = e.gather(Futures)
	e.close()
	
	# combine data
	n_X = len(number_of_neurons)
	n_classifiers = len(classifiers)
	n_trials, _ = np.shape(data)
	# pred = np.zeros([n_trials,n_classifiers,n_iter,n_X])
	pred = np.zeros([n_trials,5,n_classifiers,n_iter,n_X])
	perc_correct = np.zeros([n_classifiers,n_iter,n_X])
	for n in range(len(number_of_neurons)):
	    perc_correct[:,:,n] = out[n][0]
	    # pred[:,:,:,n] = out[n][1]
	    pred[:,:,:,:,n] = out[n][1]

	return pred.astype('bool'), perc_correct, number_of_neurons, StimLog, list(classifiers.keys())


def log_2_ID(pred, StimLog):
	n_T, _, n_classifiers, n_iter, n_X = np.shape(pred)
	n_S, _ = np.shape(StimLog)
	pred = np.moveaxis(pred,0,-1)[np.newaxis,:] # free up first dimension for broadcasting
	pred = np.tile(pred, (n_S,1,1,1,1,1)) # broadcast matrix
	StimLog = np.tile(StimLog[:,:,np.newaxis,np.newaxis,np.newaxis,np.newaxis], (1,1,n_classifiers,n_iter,n_X,n_T)) # broadcast matrix
	test = np.equal(pred, StimLog).all(axis=1) # find which elements match
	pred_id = np.apply_along_axis(np.where, 0, test).squeeze(axis=(0,1)) # find which stim was predicted for each stim
	pred_id = np.moveaxis(pred_id,-1,0) # reorganize to match original
	return pred_id


def load_output(save_file):
	h5f = h5py.File(save_file,'r')
	pred = h5f['/pred'][:]
	perc_correct = h5f['/perc_correct'][:]
	number_of_neurons = h5f['/number_of_neurons'][:]
	StimLog = h5f['/StimLog'][:]
	StimID = h5f['/StimID'][:]
	classifiers = h5f['/classifiers'][:]
	h5f.close()
	return pred, perc_correct, number_of_neurons, StimLog, StimID, classifiers


def load_mult(filenames):
	out = [load_output(f) for f in filenames]
	pred = [out[x][0] for x in range(len(out))]
	perc_correct = [out[x][1] for x in range(len(out))]
	number_of_neurons = [out[x][2] for x in range(len(out))]
	StimLog = [out[x][3] for x in range(len(out))]
	classifiers = [out[x][4] for x in range(len(out))]

	return pred, perc_correct, number_of_neurons, StimLog, classifiers


def plot_performance(perc_correct, number_of_neurons, classifiers, ax=False):
	X = number_of_neurons
	Y = np.mean(perc_correct,axis=1)
	E = 1.96*np.std(perc_correct,axis=1)

	if not ax:
		fig, ax = plt.subplots(1,1,figsize=(10,10))
	else:
		fig = ax.get_figure()

	for ind, (y,e,l) in enumerate(zip(Y,E,classifiers)):
	    ax.plot(X, y, label=l)
	    ax.fill_between(X, y-e, y+e, alpha=0.1, antialiased=True)
	#     ax.errorbar(X, y, yerr=e, label=l)
	ax.set_ylabel('Percent Correct')
	ax.set_xlabel('# of Neurons')
	ax.legend(loc='upper left')
	return fig, ax


def plot_performance_mult(perc_correct, number_of_neurons, classifiers, ax=False):
	if not ax:
		fig, ax = plt.subplots(len(perc_correct),1,figsize=(10,10))
	for a, pc, nn, c in zip(ax.flat,perc_correct,number_of_neurons,classifiers):
		plot_performance(pc, nn, c, a)
	return fig, ax


def plot_confusion_matrices(pred, StimID, number_of_neurons, C=0):
	numTrials, n_classifiers, n_iter, n_X = np.shape(pred)
	# numTrials, _, n_classifiers, n_iter, n_X = np.shape(pred)
	x = np.ceil(np.sqrt(n_X)).astype('int')
	y = np.ceil(n_X/x).astype('int')
	fig, axs = plt.subplots(y,x,sharex=True,sharey=True,figsize=(10,10))
	for ax, s in zip(axs.flat,range(n_X)):
		confusion_matrix = compute_confusion_matrix(pred, StimID, C=C, S=s)
		ax.imshow(confusion_matrix)
		ax.set_title('%d' % number_of_neurons[s])
	fig.text(0.5, 0.04, 'Predicted Label', ha='center')
	fig.text(0.04, 0.5, 'True Label', va='center', rotation='vertical')
	return fig, ax


def compute_confusion_matrix(pred, StimID, C=0, S=-1):
	_, _, n_iter, _ = np.shape(pred)
	confusion_matrix = metrics.confusion_matrix(np.repeat(StimID,n_iter), np.ravel(pred[:,C,:,S]))
	# _, _, _, n_iter, _ = np.shape(pred)
	# confusion_matrix = metrics.confusion_matrix(np.repeat(StimLog,n_iter), np.ravel(pred[:,:,C,:,S]))
	return confusion_matrix


def Savio(N):
	File = [];
	File.append('/global/scratch/elyall/Decoding/7142_220_002')
	File.append('/global/scratch/elyall/Decoding/6994_210_000')
	File.append('/global/scratch/elyall/Decoding/7120_250_003')
	File.append('/global/scratch/elyall/Decoding/7197_160_001')
	File.append('/global/scratch/elyall/Decoding/7734_338_000')
	File.append('/global/scratch/elyall/Decoding/7734_308_001')
	File.append('/global/scratch/elyall/Decoding/7736_300_000')
	File.append('/global/scratch/elyall/Decoding/7736_265_001')
	File.append('/global/scratch/elyall/Decoding/7737_291_000')
	File.append('/global/scratch/elyall/Decoding/7737_326_001')
	File.append('/global/scratch/elyall/Decoding/9445_180_005')
	File.append('/global/scratch/elyall/Decoding/9019_165_000')
	File.append('/global/scratch/elyall/Decoding/9025_180_002')
	main(File[int(N)])

	return File, N


if __name__ == "__main__":
	Savio(sys.argv[1])
	# main(sys.argv[1:])


