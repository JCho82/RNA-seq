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
Load your file, which have created in R. The folder you chose should contain the file.
  ```
  head(Mydata)
  ```
Shows the head data (a several rows) in your file.  
  ```
  DataGroups <- c("WT", "WT", "WU", "WU")
  ```
Assign your samples into a group. Since S1 and S2 (from my samples) are biological replicates, the S1 and S2 are assigned as WT.
S3 and S4 are another biological replicates so that they are grouped into WU. If you have more samples, you can assign like this pattern like ("WT", "WT", "WU", "WU", "WZ", "WZ").
                
```
d <- DGEList(counts=Mydata,group=factor(DataGroups))
```
This commands put the counts of each sample into each group.
```
d
```
Shows how they were assigned.
```
dim(d)
```
shows total featured genes and the number of samples (dimension info).
```
head(d$counts)
```
shows the number of counts again. You can compare this to the next step. 
```
head(cpm(d))
```
This shows the head counts of cpm which is calculated by counts per million. 
```
apply(d$counts, 2, sum)
```
This shows the sum of the counts for each sample.
```
keep <- rowSums(cpm(d)>100) >= 2
```
We must have at least 100 counts per million on any particular gene that we want to keep. This command allow us to keep a gene if it has a cpm of 100 or greater for at least two samples.
 
```
d <- d[keep,]
```
d is reassigned with the genes with more than 100 counts per million (from at least two samples)      
```
dim(d)
```
shows how many genes have passed through the cpm limitation.    
```
d$samples$lib.size <- colSums(d$counts)
```
This puts the counts of the filtered genes into lib.size (you can compare this to the previous).   
```
d$samples
```
shows the results.  
```
d <- calcNormFactors(d)
```
This commands put calculated normalization factors (to scale the raw library sizes) to the counts in "d."
Please look at https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors.
     
```
d
```
You can see the normalization factors assigned to each samples.
     
```
d1 <- estimateCommonDisp(d, verbose=T)
```
Assign d1 as factors to estimate a common dispersion value across all genes.
Please take a look at https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/estimateCommonDisp.  
```
names(d1)
```
shows what types of values were assigned.  
```
d1 <- estimateTagwiseDisp(d1)
```
Add another factors to estimate the tagwise negative binomial dispersions.
Please refer to https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/estimateTagwiseDisp.      
```
names(d1)
```
shows the added factors.      
```
plotBCV(d1)
```
This commands plots the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).      
https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/plotBCV      
```
et12 <- exactTest(d1, pair=c(1,2))
```
Using d1 (various statistic factor), compare the second group (2, WU) to the first group (1, WT).  
```
topTags(et12, n=10)
```
From the analyzed data above, show top 10 (n=10) genes.  
```
out1 <- topTags(et12, n=885)
```
Put all the gene which passed through the cpm limit (885, take look at values obtained from cpm command) to out1         
```
write.csv(out1, file="final.csv")
```
Save the results into csv.file.
After making this file, it would be better to save it as xlsx file in excel which is quite compatible to any other softwares.

This final data would look like below;
    
![1234](https://user-images.githubusercontent.com/105310312/210437203-41db0b5a-5a1e-4b38-b346-c4a521d4a33b.png)



        
logFC = log2 value (that is log fold change)
        
        
logCPM are the log counts per million, which can be understood as measuring expression level.
 
 
