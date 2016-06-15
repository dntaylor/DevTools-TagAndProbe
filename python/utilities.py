import os

latestTrees = {
    'electron' : '2016-06-15_DevTools_TagAndProbe_Electron_80X_v2',
    'muon'     : '2016-06-15_DevTools_TagAndProbe_Muon_80X_v1',
}

def getTagAndProbeDirectory(mode):
    baseDir = '/hdfs/store/user/dntaylor'
    if mod in latestTrees: return os.path.join(baseDir,latestTrees[mode])

