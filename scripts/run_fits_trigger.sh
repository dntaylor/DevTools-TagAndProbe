obj="$1"
inputfile="`dirname $0`/inputs_"$obj"_trigger.txt"
while read idConditions; do
    if [ ${idConditions:0:1} == "#" ]; then
        continue
    fi
    args=($idConditions)
    baseDir=/hdfs/store/user/dntaylor/tagAndProbe/$obj
    if [ "${#args[@]}" -gt "1" ]; then
        conds=$(IFS=, ; echo "${args[*]:1}")
        nohup cmsRun DevTools/TagAndProbe/python/oldfitter.py object=$obj inputFileName=$baseDir/data.root conditions=$conds idName=${args[0]} &
    else
        nohup cmsRun DevTools/TagAndProbe/python/oldfitter.py object=$obj inputFileName=$baseDir/data.root idName=${args[0]} &
    fi
done < $inputfile
