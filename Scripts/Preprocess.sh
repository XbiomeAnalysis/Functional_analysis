#! /bin/bash
# This shell script is used for generating KO profiles and Pathway profiles from KEGG annotation result of HUMANN2.

helpFunction()
{
   echo ""
   echo "Usage: $0 -i inDir -o outDir -p prefix -m method (optional)"
   echo -e "\t-i Input directory where all genefamily profile files are stored"
   echo -e "\t-o Output directory"
   echo -e "\t-p Prefix of output files"
   echo -e "\t-m Method for calculation of pathways' abundance, select from 1:median; 2:harmonic mean; 3: average of top 50% abundance KOs (default: 1)"
   echo ""
   echo "Example: sh Preprocess.sh -i /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Scripts/testIn/ \
    -o /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Scripts/testOut/ -p testme"
   echo "Note: You can either excute this script step by step in \
   /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Analysis_template/Functional_Analysis_Kegg_Pathway_preprocess.Rmd"
   echo ""
   exit 1 # Exit script after printing help
}

while getopts "i:o:p:m:" opt
do
   case "$opt" in
      i ) inDir="$OPTARG" ;;
      o ) outDir="$OPTARG" ;;
      p ) prefix="$OPTARG" ;;
      m ) method="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$inDir" ] || [ -z "$outDir" ] || [ -z "$prefix" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


## Main scripts

### Merge all genefamily profiles of HUMANN2 result
docker run -i --rm -u $(id -u):$(id -g) -v $inDir:/in -v $outDir:/out harbor.xbiome.com/xbiome/environments/humann2:2.8.1-2b8c5c3 bash    \
-c "merge_metaphlan_tables.py /in/*genefamilies.tsv > /out/merged_kegg_pathway_profile.tsv"

### Remove suffix in sample name
sed -i 's#_genefamilies##g' $outDir/merged_kegg_pathway_profile.tsv
sed -i 's$^ID$# Gene Family$g' $outDir/merged_kegg_pathway_profile.tsv

### Transform all species code with species name in merged profile file & map kegg gene ID to KO entries
python /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Scripts/rename-kegg-gene.py    \
-in $outDir/merged_kegg_pathway_profile.tsv -db /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Scripts/db/    \
-out $outDir/$prefix.kegg -species 1

### Map kegg gene entries to KO entries and generates profile table of all KOs.
python /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Scripts/reminpath.py    \
-i $outDir/$prefix.kegg -k $outDir/$prefix\_knumber.txt -m $outDir/$prefix.input

### Use minpath to extract the mapping information between Knumers and KEGG pathways
python /share/work/runtime/softwares/MinPath/MinPath.py -ko $outDir/$prefix.input    \
-report $outDir/$prefix.report -details $outDir/$prefix.details

### Calculate the abundance of Kegg pathways
python /share/projects/Analytics/analytics/Function_Analysis/Tongbangzhuo/Phase1/Kegg/Scripts/kegg_abundance.py    \
-path $outDir/$prefix.report -detail $outDir/$prefix.details    \
-in  $outDir/$prefix.kegg -out $outDir/$prefix\_pathway.txt -method ${method-1}

## Inherit UMAPPED entries' abundance to pathway profile
grep "UNMAPPED" $outDir/merged_kegg_pathway_profile.tsv >> $outDir/$prefix\_pathway.txt

### Substitute Pathway ID with Pathway names
awk -F '\t' '{OFS = "\t"} NR==FNR{a[$1] = $2} NR>FNR{if (a[$1]){$1 = a[$1]; print $0} else{print$0}}'    \
/share/work/runtime/softwares/MinPath/data/KEGG-pathway.txt $outDir/$prefix\_pathway.txt > $outDir/$prefix\_pathway_profile.txt
