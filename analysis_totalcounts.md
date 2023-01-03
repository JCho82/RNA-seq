# Create a file to analyse your data in R workspace
### This work is refered to the standford RNA seq analysis website. (https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html)
In order to make a Rdata file, you need to download and install R (https://www.r-project.org/) and Rstudio (https://posit.co/download/rstudio-desktop/).
Once you install them, start with R.


Go to packages and then click load packages. Then, select Rcmdr (if this is the first time to load Rcmdr, you will see a pop-up saying that you need to install the package. Just click one of servers listed (just select the closest server to your location), then automatically you will see it is automatically installed).
Then, you will see a pop-up box labeled with R commander. 


To load your the total counted data (in text file), you need to alter some parts.
Open excel, and then load your totalcounts.txt (the excel will suggest you how to divide your data. Just follow the suggestions (delimited -> Tab -> General). It should be OK).


Remove the unnecessay infomation like 1) the first row showing basic info 2) the second column = chr 3) the third column = start 4) the fouth column = end 5) the fifth column = strand 6) the sixth column = length. 


Modify column information like aaaaa.bam (to S1)


Since the first column is geneid, I suggest it blank (I mean the A1 position).


The final excel file should look like this.



![123](https://user-images.githubusercontent.com/105310312/210405420-99a88765-1d46-4daa-96de-9ec1c8eba797.png)

Save it as (sample.xlsx).


Now, you can work with R commander.


Go to data -> import data -> from excel file.

You will see a pop-up box.


1) enter the data name (like Mydata)


2) Mark all three 1) variable names in first row 2) Row names in first column (being set as default, it wouldn't be marked. YOU NEED TO MARK IT) 3) convert character data to factors.

3) Leave the missing data indicator as default <empty cell>  

  
  
 ![1234](https://user-images.githubusercontent.com/105310312/210407051-ac6d6d00-e4a7-45ed-8f5b-1b95530034e1.png)

  
Click View data set. Then you should be able to see the first column and row gray-colored.
  
  
Go to file -> save R workspace as -> save as type -> R data files -> type your data name (like Mydata) 
  
  
Now, you have a file to work on R studio.
  

# EdgeR
  
We are going to use EdgeR to analyze the total-counts.
  
First, you need to install edgeR in R (https://bioconductor.org/packages/release/bioc/html/edgeR.html).

  Then, file -> new project -> create project (choose either new director which you should move your Rdata file or existing directory which you have saved your Rdata file)
  
  Then type below in the terminal.
``` 
library(edgeR)
```
  This will load edgeR package for your analysis.
  
  
Then, let's start analyzing your data.
  

  ```
  load("Mydata.Rdata")
  ```
  
  ```
  head(Mydata)
  ```
  
  ```
  DataGroups <- c("WT", "WT", "WU", "WU")
  ```
                
```
d <- DGEList(counts=Mydata,group=factor(DataGroups))
```
  
```
d
```
  
```
dim(d)
```

```
head(d$counts)
```
  
```
head(cpm(d))
```
  
```
apply(d$counts, 2, sum)
```
This shows the sum of the counts for each sample.
```
keep <- rowSums(cpm(d)>100) >= 2
```
We must have at least 100 counts per million (calculated with cpm() in R) on any particular gene that we want to keep. In this example, we're only keeping a gene if it has a cpm of 100 or greater for at least two samples.
 
```
d <- d[keep,]
```
      
```
dim(d)
```
      
```
d$samples$lib.size <- colSums(d$counts)
```
  
```
d$samples
```
  
```
d <- calcNormFactors(d)
```
       
```
d
```
       
```
d1 <- estimateCommonDisp(d, verbose=T)
```
  
```
names(d1)
```
  
```
d1 <- estimateTagwiseDisp(d1)
```
        
```
names(d1)
```
      
```
plotBCV(d1)
```
      
```
et12 <- exactTest(d1, pair=c(1,2))
```
  
```
topTags(et12, n=10)
```
  
```
out1 <- topTags(et12, n=885)
```
        
```
write.csv(out1, file="final.csv")
```
          
  
  
  

 
  
  

  
  



