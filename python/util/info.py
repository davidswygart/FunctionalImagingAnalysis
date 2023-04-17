import numpy as np
from sklearn.model_selection import KFold
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler


class MI:
    def __init__(self, stim_i, ndims=2, nbins=12):
        self.ndims = ndims
        self.nbins = nbins

        self.axes = tuple(range(1,ndims+1))

        self.stim_i = stim_i
        _,self.stim_counts = np.unique(stim_i, return_counts=True)

        self.stim_prob = np.expand_dims(self.stim_counts / len(stim_i), self.axes)

        self.Py = np.empty([nbins]*ndims, dtype=np.float64)
        self.Py_x = np.empty([len(self.stim_counts)]+[nbins]*ndims, dtype=np.float64)
        self.Py_x_sh = np.empty(self.Py_x.shape, dtype=np.float64)
        # self.Py_x_ind = np.empty(self.Py_x.shape, dtype=np.float64)

        # self.edges = np.empty()    
        self.stim_c = [self.stim_i==c for c in range(len(self.stim_counts))]    

    def __call__(self, resp, shuffle=True):
        self.Py[:], edges = np.histogramdd(resp, self.nbins)
        self.Py[:] = self.Py / self.Py.sum()

        for c in range(len(self.stim_counts)):
            self.Py_x[c],*_ = np.histogramdd(resp[self.stim_c[c],:], bins=edges)

        # print(self.axes, self.Py_x.shape)
        self.Py_x[:] = self.Py_x / self.Py_x.sum(axis=self.axes, keepdims=True)

        HR = -np.nansum(self.Py * np.log2(self.Py))
        HRS = -np.nansum(self.Py_x * np.log2(self.Py_x) * self.stim_prob)

        if not shuffle:
            return HR - HRS

        for c in range(len(self.stim_counts)):
            # resp_sh = (np.random.permutation(resp[stim_c,0]), np.random.permutation(resp[stim_c,1]))
            resp_sh = [np.random.permutation(resp[self.stim_c[c],i]) for i in range(self.ndims)]
            self.Py_x_sh[c],*_ = np.histogramdd(resp_sh, bins=edges)

        self.Py_x_sh[:] = self.Py_x_sh / self.Py_x_sh.sum(axis=self.axes, keepdims=True)

        Py_x_ind = self.Py_x.sum(axis=self.axes[1:], keepdims=True)
        for d in range(1,self.ndims):
            Py_x_ind = Py_x_ind * self.Py_x.sum((*self.axes[:d], *self.axes[d+1:]), keepdims=True)

        HRS_ind = -np.nansum(Py_x_ind * np.log2(Py_x_ind) * self.stim_prob)
        HRS_sh = -np.nansum(self.Py_x_sh * np.log2(self.Py_x_sh) * self.stim_prob)

        return HR - HRS_ind + HRS_sh - HRS


    def QE(self, resp, nreps=20, nsamps=3, samps=None, shuffle=True, return_fit=False):
        if samps is None:
            samps = (resp.shape[0]//np.arange(1,nsamps+1))
        else:
            nsamps = len(samps)

        x = np.tile(1/samps, (nreps,1)).flatten()

        qe = np.empty((nsamps,nreps))

        stim_i = self.stim_i
        stim_c = self.stim_c
        #TODO: should really be a with statement...

        for i in range(nsamps):
            for j in range(nreps):
                # ri = np.random.randint(0,resp.shape[0], samps[i]) #bootstrap
                ri = np.random.permutation(resp.shape[0])[:samps[i]] #permute

                # qe[i,j] = MI_2d(resp[ri], stim_i[ri], shuffle=shuffle)
                self.stim_i = stim_i[ri]
                self.stim_c = [self.stim_i==c for c in range(len(self.stim_counts))]
                qe[i,j] = self.__call__(resp[ri])

        self.stim_i = stim_i
        self.stim_c = stim_c

        b,a,Itrue= np.polyfit(x, qe.T.flatten(), 2)
        if return_fit:
            return Itrue,a,b 
        else:   
            return Itrue

    # fig,ax = plt.subplots(1,3, figsize=(30,5))
    # ax[0].plot(Itrue)
    # ax[1].plot(a)
    # ax[2].plot(b)

def svm_info(data, labels, nlabels, ncomps=10, **kwargs):
    cv = KFold(10, shuffle=True)

    y_hat = np.empty(len(labels))
    for train, test in cv.split(data):
        
        # mn = data[train].min(axis=(0,1), keepdims=True)
        # mx = data[train].
        # max(axis=(0,1))

        # mmd = (data[train] - mn) / (mx-mn)
        mms = MinMaxScaler().fit(data[train])
        mod = SVC(**kwargs)

        if ncomps < data.shape[1] and ncomps < len(train):
            pca = PCA(ncomps, whiten=True)
            comps = pca.fit_transform(mms.transform(data[train]))
            
            mod.fit(comps, labels[train])


            y_hat[test] = mod.predict(pca.transform(mms.transform(data[test])))
        else:
            mod.fit(data[train], labels[train])
            y_hat[test] = mod.predict(mms.transform(data[test]))



    mi = MI(labels, 1, nbins=nlabels) #TODO: this really needs to be categorical...
    return mi.QE(y_hat[:,None], shuffle=False)


# def MI_2d_cont(X,Y, nbins=12):
#     Px, x_e = np.histogramdd(X, nbins)
#     Px = Px / Px.sum()

#     Py, y_e = np.histogramdd(Y, nbins)
#     Py = Py / Py.sum()
