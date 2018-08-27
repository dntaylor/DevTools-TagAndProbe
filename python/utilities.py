import os

latestTrees = {
    'electron' : '2016-06-15_DevTools_TagAndProbe_Electron_80X_v2',
    'muon'     : '2016-06-15_DevTools_TagAndProbe_Muon_80X_v1',
}

def getTagAndProbeDirectory(mode):
    baseDir = '/hdfs/store/user/dntaylor'
    if mod in latestTrees: return os.path.join(baseDir,latestTrees[mode])

def getBinning(obj,var,trig=False,threshold=20):
    '''Return the binning for each object variable'''
    if var=='pt' and trig:
        #bins = [5]
        #for i in range(threshold-5,threshold+10):
        #    if i<=5: continue
        #    bins += [i]
        #bins += [threshold+15, threshold+20, threshold+30, threshold+40]
        #if bins[-1]<50: bins += [50]
        #if bins[-1]<100: bins += [100]
        #if bins[-1]<500: bins += [500]

        # adaptive
        dpt = max([int(threshold/20),1])
        bins = [0]
        for i in range(threshold-5*dpt, threshold+11*dpt, dpt):
            if i: bins += [i]
        for i in range(threshold+12*dpt, threshold+22*dpt, 2*dpt):
            if i: bins += [i]
        for i in range(threshold+25*dpt, threshold+55*dpt, 5*dpt):
            if i: bins += [i]
        if bins[-1]<100: bins += [100]
        if bins[-1]<500: bins += [500]

        return bins
    if obj=='electron':
        if var=='pt':
            return [10, 20, 30, 40, 50, 200]
        if var=='eta':
            return [-2.5, -2.0, -1.479, -0.8, 0., 0.8, 1.479, 2.0, 2.5]
    if obj=='muon':
        if var=='pt':
            return [5, 10, 20, 30, 40, 50, 200]
        if var=='eta':
            return [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]
    if obj=='photon':
        if var=='pt':
            #return [10, 15, 20, 25, 30, 40, 50, 100]
            return [20, 35, 50, 90, 150, 500]
        if var=='eta':
            #return [-2.5, -2.0, -1.479, -0.8, 0., 0.8, 1.479, 2.0, 2.5]
            return [-2.5, -2.0, -1.57, -1.4442, -0.8, 0., 0.8, 1.4442, 1.57, 2.0, 2.5]
    if obj=='charge':
        if var=='pt':
            return [10, 20, 30, 40, 50, 100, 200, 500, 6500]
        if var=='eta':
            return [0., 0.8, 1.479, 2.0, 2.5]
    return []
