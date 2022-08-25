#!/usr/bin/env python                                                                                                         
import sys
import gzip
import os
import re
import pandas as pd
import numpy as np
import statistics 

def loadMinpath(inputFile,pathwayType):
    inputF=open(inputFile,"r")
    reportPath=list()
    for inline in inputF:
        infor=inline.strip().split(" ")
        for num, name in enumerate(infor, start=0):
            if name==pathwayType:
                if infor[num+1]=="1":
                    reportPath.append(infor[1])
    inputF.close()
    return (reportPath)

def FromDetailsToTable(input_details):
    res_dict = dict()
    pathwayTotalDict=dict()
    f = open(input_details,'r')
    for line in f.readlines():
        line = line.strip()
        infor=line.split(" ")
        if  infor[0]== "path":
            key=infor[1]
            pathwayTotalDict[key]=int(infor[3])
            res_dict[key] = list()
        else:
            res_dict[key].append(infor[0])
    f.close()
    return (res_dict,pathwayTotalDict)

def checkLine(inputLine):    
    if inputLine.startswith("#") or inputLine.startswith("UNMAPPED") or (not inputLine.startswith("K")):
        return False
    if "|" in inputLine or "," in inputLine:
        return False
    return True


def calListAbun(inputList,method):
    if method=="harmonic_mean":
        #print (inputList)
        inputList=list(filter(lambda a: a !=0, inputList))  #harmonic value will be 0 if any 0 present in the list
        #print (inputList)
        if len(inputList)==0:
            return (0)
        else:
            return (statistics.harmonic_mean(inputList))
    elif method=="median":
        return (statistics.median(inputList))
    elif method=="majority":
        inputList.sort(reverse=True)
        lst=inputList[:(int(len(inputList)/2))]
        #print (inputList)
        #print (lst)
        return(float(sum(lst) / len(lst)))
    else:
        return (0)
    
def calPathAbundance(KODict,pathwayGeneSet,method,pathwayTotalDict):
    pathwayAbuDict=dict()
    for (pathwayID,geneSet) in pathwayGeneSet.items():
        GeneSetAbun=list(map(lambda x:KODict[x] if x in KODict.keys() else 0,geneSet))
        GeneSetAbun=GeneSetAbun+[0]*(pathwayTotalDict[pathwayID]-len(GeneSetAbun))
        #print (pathwayID)
        pathwayAbu=calListAbun(GeneSetAbun,method)
        #print (pathwayID+":"+str(pathwayAbu))
        pathwayAbuDict[pathwayID]=pathwayAbu
    return (pathwayAbuDict)

        
def loadKOsample(inputFile,pathwayGeneSet,reportPath,method,outputFile,pathwayTotalDict):
    inputF=open(inputFile,"r")
    outputF=open(inputFile+".tmp","w")
    reportPathDict=dict()
    commentLine=""
    for inline in inputF:
        if inline.startswith("#"):
            commentLine+=inline
        if checkLine(inline):
            print (inline.strip(),file=outputF)
    inputF.close()
    outputF.close()
    knumber = pd.read_csv(inputFile+".tmp",sep="\t",header=None)
    outR=open(outputFile,"w")
    for (columnName, columnData) in knumber.iteritems():
        if columnName>0:
            mydict = dict(zip(knumber[0], columnData.values))
            pathwayAbundanceDict=calPathAbundance(mydict,pathwayGeneSet,method,pathwayTotalDict)
            for eachPath in reportPath:
                if eachPath not in reportPathDict:
                    reportPathDict[eachPath]=[]
                if eachPath in pathwayAbundanceDict:
                    reportPathDict[eachPath].append(str(pathwayAbundanceDict[eachPath]))
                else:
                    reportPathDict[eachPath].append("0")
    print (commentLine.strip(),file=outR)
    for (pathway,abu) in reportPathDict.items():
        print (pathway+"\t"+"\t".join(abu),file=outR)
    os.remove(inputFile+".tmp")
    outR.close()
    
# 00901006_k_number.txt文件为用 全姐的脚本处理 humann2 kegg 结果后得到的文件，需要将knumber对应的行提取出来
# 因为新技术报告是两个样本，所以只有3列，如果有多个样本，names 需要按照样本数修改。
#k_number = pd.read_table("./00901006_k_number.txt",sep="\t",names=["Knumber","before","after"])

# 将details文件的格式整理为如下格式：
# KO    Knumber
# 00010 k000001
# 00010 k000002
# 00010 k000003
# .............
#lst = []
#for k,v in input_details.items():
#    for i in v:
#        lst.append((k,i))
#details_data = pd.DataFrame(lst,columns=["KO","Knumber"])

# 最后将基因丰度文件与KO-Knumber文件merge
#merge_data = pd.merge(k_number,details_data,left_on="Knumber",right_on="Knumber")

# 接下来 groupby KO 就可以得到每个KO 的相对丰度
#data_groupby_KO = merge_data.groupby("KO").sum().apply(lambda x: x*100/sum(x))

def main():
    usage="python3 kegg_abundance.py -path <minpath report from minPath> -detail <minpath detail from minPath> -in <KO table from processed humann2> -out <output file> -report <1:minpath reported 2: naive pathways> -method <1:median; 2:harmonic mean; 3: average of top 50% abundance KOs >\n"
    usage+="parameters:\n"
    usage+='%-10s' % ("-path:")+"file .ko.minpath generated by Minpath \n"
    usage+='%-10s' % ("-detail:")+"file .ko.minapth.detail generated by Minpath \n"
    usage+='%-10s' % ("-in:")+" input file that is generated by processed humann2 \n"
    usage+='%-10s' % ("-out:")+" output file (default: output.tsv) \n"
    usage+='%-10s' % ("-report:")+"1: the pathway reported by minpath; 2: naive pathways (default: 1)\n"
    usage+='%-10s' % ("-method:")+"The method to calculate the pathway abundance: 1:median; 2:harmonic mean; 3: average of top 50% abundance KOs (default: 1)\n"
    methodDict={1:"median",2:"harmonic_mean",3:"majority"}
    reportDict={1:"minpath",2:"naive"}  #naive 1  minpath 
    
    inputFile=None
    minpath=None
    minpathDetail=None
    report=1
    method=1
    output="output.tsv"

    for idx in range(len(sys.argv)):
        if (sys.argv[idx] == "-in") and (len(sys.argv) > idx + 1):
            inputFile=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-h") or (sys.argv[idx] == "-help"):
            print (usage)
            sys.exit()
        elif (sys.argv[idx] == "-path") and (len(sys.argv) > idx + 1):
            minpath=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-out") and (len(sys.argv) > idx + 1):
            output=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-detail") and (len(sys.argv) > idx + 1):
            minpathDetail=sys.argv[idx + 1]
        elif (sys.argv[idx] == "-report") and (len(sys.argv) > idx + 1):
            report=int(sys.argv[idx + 1])
        elif (sys.argv[idx] == "-method") and (len(sys.argv) > idx + 1):
            method=int(sys.argv[idx + 1])
            
    if report in reportDict and method in methodDict:
        report=reportDict[report]
        method=methodDict[method]
    else:
        print (usage)
        sys.exit()
    if inputFile and minpath and minpathDetail:
        pathwayResult=loadMinpath(minpath,report)
        (KOinPathway,pathwayTotalDict)=FromDetailsToTable(minpathDetail)
        loadKOsample(inputFile,KOinPathway,pathwayResult,method,output,pathwayTotalDict)
    else:
        print (usage)
        sys.exit()


if __name__ == '__main__':
    main()