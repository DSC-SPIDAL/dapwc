#!/bin/bash
. init.properties
numNodes=$1
nodesFile=$2
confFileDAPWC=$6
confFileDAMDS=$7


mdstpp=$3
mdsppn=$4
mdsmem=$5
distFile=$8
roundName=round_
roundCount=0
fpwc=pwc
fmds=mds
mdsName=results
mdsFolderName=$fmds"_"$mdsName
doNext=y
trap '' 2
################# Run MDS to get points for visualization ####################
echo "Start Running MDS algorithm on `date`" >> status.txt
#./run.generic.sh $mdstpp $mdsppn $fmds $mdsName $mdsmem g $numNodes $nodesFile $confFileDAMDS
echo "Finised Running MDS algorithm on `date` " >> status.txt
################# End MDS Section #########################

################# start dapwc section #########################
folderName=$fpwc"_"$roundName$roundCount
mkdir -p $folderName
echo "Start Running DAPWC initial round algorithm on `date`" >> status.txt
#./run.dapwc_init.sh $mdstpp $mdsppn $fpwc $roundName$roundCount $mdsmem g $numNodes $nodesFile $confFileDAPWC
echo "Start Running DAPWC initial round algorithm on `date`" >> status.txt
################# End dapwc Section #########################

## create webplotviz files for initial run
mv ./$folderName/cluster* ./$folderName/cluster.txt 
./create_webplotvix.sh $PWD/$folderName/cluster.txt $PWD/$mdsFolderName/points.txt $PWD/$folderName/plot.webplotviz

echo "Created webplotviz file for initial run `date`" >> status.txt
## create default group in WebPlotViz
curl -i -H "Accept: application/json" -X POST https://spidal-gw.dsc.soic.indiana.edu/groups/public/dapwc/create

## upload file to WebPlotViz
desc="Plot for "$runName" initial run"
cp $PWD/$folderName/plot.webplotviz $PWD/$folderName/$runName".initial.webplotviz"
curlFile=$PWD/$folderName/$runName".initial.webplotviz"
curl -F "desc=$desc" -F "file=@$curlFile" https://spidal-gw.dsc.soic.indiana.edu/upload/public/dapwc

echo -e "\e[34mWebPlotViz file \e[92m$PWD/$folderName/plot.webplotviz \e[34mcreated and uploaded. Please check the plot in WebPlotViz DAPWCTempCollection Collection\e[0m"
echo ""
echo -e "\e[93mcheck plot at https://spidal-gw.dsc.soic.indiana.edu/groupdashboard/DAPWCTempCollection \e[0m"
echo ""
echo -e "\e[41mPlease provide clusters that need to be merged in the following format 0+2|4+5|. Join 0 and 2. Also join 4 and 5\e[0m"
echo -e "\e[41mEnter" "|" if no joins are needed"\e[0m"

read pattern

### Check the cluster file and define clusters to merge #################
./join.sh $PWD/$folderName/cluster.txt $PWD/$mdsFolderName/points.txt $PWD/$folderName/cluster_joined.txt $PWD/$folderName/plot_joined.webplotviz $pattern

desc="Plot for "$runName" initial run joined"
cp $PWD/$folderName/plot_joined.webplotviz $PWD/$folderName/$runName".initial.joined.webplotviz"
curlFile=$PWD/$folderName/$runName".initial.joined.webplotviz"
curl -F "desc=$desc" -F "file=@$curlFile" https://spidal-gw.dsc.soic.indiana.edu/upload/public/dapwc
echo -e "\e[34mWebPlotViz file \e[92m$PWD/$folderName/plot_joined.webplotviz \e[34mcreated and uploaded. Please check the plot in WebPlotViz DAPWCTempCollection Collection\e[0m"
echo ""
echo -e "\e[93mcheck plot at https://spidal-gw.dsc.soic.indiana.edu/groupdashboard/DAPWCTempCollection \e[0m"
echo ""

echo -e "\e[41mDo you wish to continue y or n\e[0m"
read doNext

resultsFolderOld=results_$roundCount
resultsFolderNew=""
pwcFolderOld=$fpwc"_"$roundName$roundCount
pwcFolderNew=""
while [ "$doNext" == "y" ] 
do
echo -e "\e[31m Starting DAPWC round $roundCount \e[0m"
echo  -e "\e[32mEnter clusters that need to be further clustered\e[0m"
echo -e "\e[41mEnter comma seperated list of clusters to run DAPWC on \e[0m"
read commaclusters
echo -e "\e[41mEnter comma seperated list to specifiy the # of sub clusters for each cluster\e[0m"
read commaSubclusters

echo running cluster extractor 
############# run cluster extractor ##################
folder=results_$roundCount
mkdir -p $resultsFolderOld
((roundCount++))
resultsFolderNew=results_$roundCount
pwcFolderNew=$fpwc"_"$roundName$roundCount

./extract.sh $distFile $PWD/$pwcFolderOld/cluster_joined.txt $commaclusters $commaSubclusters $PWD/$resultsFolderOld

echo running DAPWC for each cluster
./run_all_in_folder.sh $mdstpp $mdsppn $mdsmem $numNodes $nodesFile $PWD/$resultsFolderOld/ $roundCount

echo ran DAPWC for all specified clusters results are available in $PWD/$resultsFolderOld/
echo combining results to single file. plot files will be available for each spereated clusters as well. 
newFolderName=$fpwc"_"$roundName$roundCount
. $PWD/$resultsFolderOld/conf_shared.prop
./collate.sh $PWD/$pwcFolderOld/plot_joined.webplotviz $PWD/$pwcFolderNew/$pwcFolderNew"_cluster_{0}" cluster-M"{0}"-C"{1}"txt $PWD/$pwcFolderNew/ $numClusters $subClusters
echo -e "\e[34mcollated cluster file \e[92m$PWD/$pwcFolderNew/cluster.txt\e[0m"
echo -e "\e[34mcollated plot file \e[92m$PWD/$pwcFolderNew/plot.webplotviz\e[0m"
## upload files
desc="Plot for "$runName" round "$roundName$roundCount" run"
cp $PWD/$pwcFolderNew/plot.webplotviz $PWD/$pwcFolderNew/$runName"."$roundName$roundCount".webplotviz"
curlFile=$PWD/$pwcFolderNew/$runName"."$roundName$roundCount".webplotviz"
curl -F "desc=$desc" -F "file=@$curlFile" https://spidal-gw.dsc.soic.indiana.edu/upload/public/dapwc

echo -e "\e[34mWebPlotViz file \e[92m$PWD/$pwcFolderNew/plot.webplotviz \e[34mcreated and uploaded. Please check the plot in WebPlotViz DAPWCTempCollection Collection\e[0m"
echo ""
echo -e "\e[93mcheck plot at https://spidal-gw.dsc.soic.indiana.edu/groupdashboard/DAPWCTempCollection \e[0m"
echo ""

echo -e "\e[41mPlease provide clusters that need to be merged in the following format 0+2|4+5|. Join 0 and 2. Also join 4 and 5\e[0m"
echo -e "\e[41m Enter | if no joins are needed \e[0m"

read pattern
./join.sh $PWD/$pwcFolderNew/cluster.txt $PWD/$pwcFolderNew/plot.webplotviz $PWD/$pwcFolderNew/cluster_joined.txt $PWD/$pwcFolderNew/plot_joined.webplotviz $pattern
resultsFolderOld=$resultsFolderNew
pwcFolderOld=$pwcFolderNew
echo done joining clusters
echo -e "\e[34mJoined cluster file \e[92m$PWD/$pwcFolderNew/cluster_joined.txt\e[34m and plot file \e[92m$PWD/$pwcFolderNew/plot_joined.webplotviz\e[34m created Check plot to identify clusters that need to be seperated. Execute following command on local machine to download plot file\e[0m"
## upload files
desc="Plot for "$runName" round "$roundName$roundCount" joined run"
cp $PWD/$pwcFolderNew/plot_joined.webplotviz $PWD/$pwcFolderNew/$runName"."$roundName$roundCount".joined.webplotviz"
curlFile=$PWD/$pwcFolderNew/$runName"."$roundName$roundCount".joined.webplotviz"
curl -F "desc=$desc" -F "file=@$curlFile" https://spidal-gw.dsc.soic.indiana.edu/upload/public/dapwc

echo -e "\e[34mWebPlotViz file \e[92m$PWD/$pwcFolderNew/plot.webplotviz \e[34mcreated and uploaded. Please check the plot in WebPlotViz DAPWCTempCollection Collection\e[0m"
echo ""
echo -e "\e[93mcheck plot at https://spidal-gw.dsc.soic.indiana.edu/groupdashboard/DAPWCTempCollection \e[0m"
echo ""

echo -e "\e[41mDo you wish to continue with another round ? y or n\e[0m"
read doNext
done
trap 2
echo Done processing
