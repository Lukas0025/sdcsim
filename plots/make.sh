samples=1
iterations=30

timeC=99999
timeRange=1
timeStep=1

strandsC=100
strandsRange=1
strandsStep=1

python3 ../stats.py ./../build/sdcsim ./../sdcasm/bind.sdcasm            --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p bind.png
#python3 ../stats.py ./../build/sdcsim ./../sdcasm/bind_and_unwrap.sdcasm --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p bind_and_unwrap.png
python3 ../stats.py ./../build/sdcsim ./../sdcasm/replace_assist.sdcasm  --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p replace_assist.png
python3 ../stats.py ./../build/sdcsim ./../sdcasm/replace.sdcasm         --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p replace.png
python3 ../stats.py ./../build/sdcsim ./../sdcasm/sticker.sdcasm         --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p sticker.png
python3 ../stats.py ./../build/sdcsim ./../sdcasm/unbind_complete.sdcasm --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p unbind_complete.png
python3 ../stats.py ./../build/sdcsim ./../sdcasm/unbind_partly.sdcasm   --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p unbind_partly.png
python3 ../stats.py ./../build/sdcsim ./../sdcasm/rule110.sdcasm         --strands $strandsC --strands_range $strandsRange --strands_step $strandsStep -t $timeC --time_range $timeRange --time_step $timeStep -T 50 -r 100 -s $samples -i $iterations -p rule110.png