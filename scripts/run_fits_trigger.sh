obj="$1"
inputfile=inputs_"$obj"_trigger.txt
while read idConditions; do
    if [ ${idConditions:0:1} == "#" ]; then
        continue
    fi
    args=($idConditions)
    if [ "${#args[@]}" -gt "1" ]; then
        conds=$(IFS=, ; echo "${args[*]:1}")
        nohup cmsRun DevTools/TagAndProbe/python/oldfitter.py object=$obj inputFileName=tagAndProbe/$obj/data.root conditions=$conds idName=${args[0]} &
    else
        nohup cmsRun DevTools/TagAndProbe/python/oldfitter.py object=$obj inputFileName=tagAndProbe/$obj/data.root idName=${args[0]} &
    fi
done < $inputfile
