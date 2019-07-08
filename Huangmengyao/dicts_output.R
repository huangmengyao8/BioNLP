# 需要加载的包
library(foreign)
library(dplyr)
library(data.table)
library(stringr) 
library(bindrcpp)
library(tidyr)

#设置工作路径
#setwd("D:\\工作\\平安实习\\Bionlp\\BioNLP\\tax")
#path <- "D:\\工作\\平安实习\\Bionlp\\BioNLP\\tax\\"

#========================taxonomy
#####################(1)存储格式1
namesdata <- read.csv("names.dmp",header=F,sep = "|",quote="")
namesdata <- namesdata[,-5]                #去除空列
colnames(namesdata) <- c("tax_id","name_txt","unique_name","name_class")
namesdata <- namesdata[,-3]                #仅保留id name,name_class
namesdata$tax_id <- as.character(namesdata$tax_id) #变成字符型
namesdata$tax_id <- str_trim(namesdata$tax_id,side="both") #去除字符前后的空格
namesdata$name_txt <- as.character(namesdata$name_txt)
namesdata$name_txt <- str_trim(namesdata$name_txt,side="both")
namesdata$name_class <- as.character(namesdata$name_class)
namesdata$name_class <- str_trim(namesdata$name_class,side="both")

#写入文件，以|分隔字段
#write.table(namesdata,paste(path,"TAX1.txt",sep=""),sep="|",quote=F,row.names = F)

#taxonomy_trim
namesdata1 <- namesdata
normids <- read.table("BioNLP-OST-2019_BB-norm_Microorganism-ids.txt",header=F,sep = "",quote="")
#normids <- data.frame(normids[600001:400000,])
colnames(normids) <- "tax_id"
namesdata1 <- merge(namesdata1,normids,by="tax_id")
#write.table(namesdata,paste(path,"tax1_trim.txt",sep=""),sep="|",quote=F,row.names = F)

#####################(2)存储格式2
###合并name,name_class
temp <- namesdata1[,1:2]
temp$tax_id <- as.numeric(temp$tax_id)
temp1 <- namesdata1[,c(1,3)]
tax_id <- temp$tax_id
name_txt <- temp$name_txt
name_class<- temp$name_class
n <- n_distinct(tax_id) 

newdata1 <- as.data.frame(matrix(NA,ncol=2,nrow = n))
colnames(newdata1) <- c("tax_id","name_txt")
newdata2 <- as.data.frame(matrix(NA,ncol=2,nrow = n))
colnames(newdata2) <- c("tax_id","name_class")

system.time(
  for(i in 1:n){
    data <- filter(temp,tax_id %in% unique(tax_id)[i])
    newdata <- as.data.frame(t(c(data[1,1],data[,2])))
    colnames(newdata)[1] <- "tax_id"
    newdata <- unite(newdata,"name_txt",-tax_id)
    
    newdata1[i,1] <- unique(tax_id)[i]
    newdata1[i,2] <- newdata$name_txt
    
    data1 <- filter(temp1,tax_id %in% unique(tax_id)[i])
    newdata3 <- as.data.frame(t(c(data1[1,1],data1[,2])))
    colnames(newdata3)[1] <- "tax_id"
    newdata3 <- unite(newdata3,"name_class",-tax_id)
    
    newdata2[i,1] <- unique(tax_id)[i]
    newdata2[i,2] <- newdata3$name_class
    print(i)
  }
)

TAX2 <- na.omit(cbind(newdata1[,1:2],newdata2[,2]))
colnames(TAX2) <- c("tax_id","name_txt","name_class")
write.table(TAX2,paste(path,"TAX2_trim747643_1081798.txt",sep=""),sep="|",quote=F,row.names = F)

#========================OntoBiotope
OBTdata <- read.table("OntoBiotope_BioNLP-OST-2019.obo",header=F,sep = "\t",quote="",skip = 3)
OBTdata$V1 <- as.character(OBTdata$V1)
id <- OBTdata[grep(pattern="id:",OBTdata[,1]),] #取出id:开头的行
id <- substr(id,5,1000)                         #仅保留id码
name <- OBTdata[grep(pattern="name:",OBTdata[,1]),] #取出name:开头的行
name <- substr(name ,7,1000)                        #仅保留name名
id <- as.character(id)
name <- as.character(name)
OBTdata1 <- data.frame(id,name)

#写入文件，以|分隔字段
#write.table(OBTdata1,paste(path,"OBT.txt",sep=""),sep="|",quote=F,row.names = F)


