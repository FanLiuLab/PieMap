##################################################drawPiemap##################################################
#'@date 2019-08-13
#'@author Li Yi
#'@title draw the map of worldwide allele-frequency distributions for SNPs based on 1000 Genome Project
#'@description This function can calculate the allele-frequency for SNPs using data of 1000 Genome Project and generate the csv file that is prepared for Mapviewer software. And can draw the map plot of worldwide allele-frequency distributions for SNPs using R.
#'@environment R;r packages ("rworldmap");this script applies to the sever:192.168.118.92. 
#'@param workpath (character string)the absolute URL of the output files.
#'@param effect (character string)the absolute URL of text file including SNP information. This file did not including header. The first column of this file was rsid of SNPs. The second column of this file was the effect allele.
#'@param effectAllele (logical)if TRUE, change the allele of SNPs in 1000 Genome Project based the effect allele you supply, default F.
#'@param piecolor (vector of character string)choose the two colors for pie plot, default is blue and yellow.
#'@param landcolor (character string)choose the color for land of map plot, default is grey79.
#'@return Generate the csv file of the allele-frequency for SNPs and the map plot of this data.
#'@export The header of export csv file. pop: population; freq: the frequency of A1; num: the number of alleles; A1 ; A2.
#'@example 
#'workpath="/thinker/storage/org/liufanGroup/public_software/share/mapPiePDF-output"
#'effect="/thinker/storage/org/liufanGroup/public_software/share/mapPiePDF-example.txt"
#'source("/thinker/storage/org/liufanGroup/public_software/share/mapPiePDF-ly.R")
#'drawPiemap(workpath,effect,effectAllele=FALSE)
drawPiemap<-function(workpath,effect,piecolor=c("blue","yellow"),landcolor="grey79",effectAllele=FALSE){
    library("rworldmap")
    setwd(workpath)
    effectinfo<-read.table(effect,header=F,stringsAsFactors=F)
    if(dim(effectinfo)[2]==2){
        colnames(effectinfo)<-c("SNP","EFA")
    }
    coorinfo<-read.csv("/thinker/storage/org/liufanGroup/liufan/data/1000G/sample_coorINFO.csv",header=T,stringsAsFactors=F)
    mysample<-read.table("/thinker/storage/org/liufanGroup/liufan/data/1000G/sample.txt",header=T)
    snplist<-effectinfo[,1]
    write.table(snplist,file=paste(workpath,"/snplist.txt",sep=""),row.names=F,quote=F,col.names=F)
    system(paste("plink --bfile /thinker/storage/org/liufanGroup/liufan/data/1000G/1000G --extract ",workpath,"/snplist.txt"," --make-bed --recodeA --out ",workpath,"/extract_result",sep=""))
    mydata1<-read.table("extract_result.bim",header=F,stringsAsFactor=F)
    mydata2<-read.table("extract_result.raw",header=T,,stringsAsFactor=F)
    mydata2$FID<-as.character(mydata2$FID)
    mydata2$IID<-as.character(mydata2$IID)
    for(i in 7:dim(mydata2)[2])
        {
        mydata2[,i][is.na(mydata2[,i])]=mean(mydata2[,i],na.rm=TRUE)
        }
    pop<-c("GBR","FIN","CHS","PUR","CDX","CLM","IBS","PEL","PJL","KHV","ACB","GWD","ESN","BEB","MSL","STU","ITU","CEU","YRI","CHB","JPT","LWK","ASW","MXL","TSI","GIH")
    mydata3<-mydata2[,-(2:6)]
    mydata4<-merge(mysample,mydata3,by.x="sample",by.y="FID")
    num<-numeric(length(pop))
    freq<-numeric(length(pop))
    efa<-numeric(length(pop))
    otha<-numeric(length(pop))
    mydata4$sample<-as.character(mydata4$sample)
    mydata4$pop<-as.character(mydata4$pop)
    mydata6<-mydata1
    if(effectAllele==FALSE)
    {
        for(j in 3:dim(mydata4)[2])
        {
            for(i in 1:length(pop))
            {
                num[i]<-length(which(mydata4$pop==pop[i]))
                med<-c(which(mydata4$pop==pop[i]))
                freq[i]<-(sum(mydata4[med,j]))/2/num[i]
            }
            efa<-freq*num
            otha<-num-efa
            mydata5<-data.frame(pop,freq,num,efa,otha)
            mydata5$pop<-as.character(mydata5$pop)
            mydata7<-mydata6
            info<-c(mydata7[j-2,2],paste("ref_allele ",mydata7[j-2,5],sep=""),paste("other_allele ",mydata7[j-2,6],sep=""))
            colnames(mydata5)<-c("pop","freq","num",info[2],info[3])
            write.csv(mydata5,file=paste("Ref_map_",info[1],"_",info[2],".csv",sep=""),quote=F,row.names=F)
            drawdata<-merge(mydata5,coorinfo,by.x="pop",by.y="pop")
            pdf(paste(workpath,"/Ref_",info[1],"_mapPies.pdf",sep=""),bg="white",width=14.7,height=7)
            mapPies(drawdata,nameX = "Longitude",nameY="Latitude",zColours = piecolor,nameZs = c(info[2],info[3]),addCatLegend = F,mapRegion="world",barOrient='vert',oceanCol="white",landCol=landcolor,borderCol="white")
            points(x=c(-172,-172),y=c(-32,-41),pch=19,col=piecolor,cex=2.9)
            legend(x=c(-177,-177), y = c(-23,-38), c(info[2],info[3]),cex=1.15,bty="n")
            text(x=-152,y=-20,labels=info[1],cex = 1.3)
            dev.off()
        }
    }else
    {
        sort<-c(1:dim(mydata6)[1])
        mydata6<-cbind(mydata6,sort)
        mydata7<-merge(mydata6,effectinfo,by.x="V2",by.y="SNP")
        mydata7<-mydata7[order(mydata7$sort),]
        index<-which(mydata7$V5!=mydata7$EFA)
        for(m in 1:length(index))
        {
            mydata4[,index[m]+2]<-abs(mydata4[,index[m]+2]-2)
        }
        othA<-mydata7$V6
        othA[index]<-mydata7$V5[index]
        mydata7<-cbind(mydata7,othA)
        for(j in 3:dim(mydata4)[2])
        {
            for(i in 1:length(pop))
            {
                num[i]<-length(which(mydata4$pop==pop[i]))
                med<-c(which(mydata4$pop==pop[i]))
                freq[i]<-(sum(mydata4[med,j]))/2/num[i]
            }
            efa<-freq*num
            otha<-num-efa
            mydata5<-data.frame(pop,freq,num,efa,otha)
            mydata5$pop<-as.character(mydata5$pop)
            info<-c(mydata7[j-2,1],paste("Effect Allele ",mydata7[j-2,8],sep=""),paste("Non-Effect Allele ",mydata7[j-2,9],sep=""))
            colnames(mydata5)<-c("pop","freq","num",info[2],info[3])
            write.csv(mydata5,file=paste("map_",info[1],"_",info[2],".csv",sep=""),quote=F,row.names=F)
            drawdata<-merge(mydata5,coorinfo,by.x="pop",by.y="pop")
            pdf(paste(workpath,"/Efa_",info[1],"_mapPies.pdf",sep=""),bg="white",width=14.7,height=7)
            mapPies(drawdata,nameX = "Longitude",nameY="Latitude",zColours = piecolor,nameZs = c(info[2],info[3]),addCatLegend = F,mapRegion="world",barOrient='vert',oceanCol="white",landCol=landcolor,borderCol="white")
            points(x=c(-172,-172),y=c(-32,-41),pch=19,col=piecolor,cex=2.9)
            legend(x=c(-177,-177), y = c(-23,-38), c(info[2],info[3]),cex=1.15,bty="n")
            text(x=-152,y=-20,labels=info[1],cex = 1.3)
            dev.off()
        }
    }
    rmfiles <- list.files(workpath,pattern="extract_result")
    for(fil in rmfiles)
    {
        system(paste0("rm -f ",workpath,"/",fil))
    }
    system(paste0("rm -f ",workpath,"/snplist.txt"))
}
workpath="/thinker/storage/org/liufanGroup/public_software/share/mapPiePDF-output"
effect="/thinker/storage/org/liufanGroup/public_software/share/mapPiePDF-example.txt"
drawPiemap(workpath,effect,effectAllele=FALSE)