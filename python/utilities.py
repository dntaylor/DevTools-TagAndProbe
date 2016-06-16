import os

latestTrees = {
    'electron' : '2016-06-15_DevTools_TagAndProbe_Electron_80X_v2',
    'muon'     : '2016-06-15_DevTools_TagAndProbe_Muon_80X_v1',
}

def getTagAndProbeDirectory(mode):
    baseDir = '/hdfs/store/user/dntaylor'
    if mod in latestTrees: return os.path.join(baseDir,latestTrees[mode])

def getBinning(obj,var):
    '''Return the binning for each object variable'''
    if obj=='electron':
        if var=='pt':
            return [10, 20, 30, 40, 50, 200]
        if var=='eta':
            return [-2.5, -2.0, -1.479, -0.8, 0., 0.8, 1.479, 2.0, 2.5]
    if obj=='muon':
        if var=='pt':
            return [10, 20, 30, 40, 50, 200]
        if var=='eta':
            return [-2.4, -2.1, -1.6, -1.2, -0.9, -0.3, -0.2, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 2.4]
    return []
