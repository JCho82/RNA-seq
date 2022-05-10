# Environment settings

First, you should be able to use all codes in Conda environment.
If you don't have it, please contact your server administrator since intalling the Conda/anaconda would be prohibited by yourself.

### When you have access to conda

You may have to write code like this or your server administrator may set up a path for you to use conda.

```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

In order not to run this command each time you start a new session on discovery, it would be better to set it in a path.

nano is an editor to modify a file.
```
nano .bashrc
```
In the editor, add the code for source in .bashrc file.

Then, you need to make folders to store your own conda environments

```
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```
-p represents -path or -parents.

If you have an environment file with yml, you can directly install your environment.

```
conda env create -f /dartfs-hpc/scratch/rnaseq1/environment.yml
```

Otherwise, you should make your environment by yourself.

```
conda create --name JC82env
```
Then activate your environment.
```
conda activate JC82env
```
In your environment, install your programs of interest such as python, cutadapter and etc.
If you want to make your own env.yml to share others, type the code below in the activated env.

```
conda env export > /home/JC/JC82env.yml
```
### Something useful to work in command line.

If your files are too big but in your server, please use "link" command. To unlink, type unlink with targets.
./ indicates your current position.
```
ln -s /dartfs-hpc/rc/lab/G/GSR_Active/Labs/Cheung/RNA-seq_3-7-22/*fastq.gz ./
```
To unlink,
```
unlink *fastq.gz
```

if you feel sick of typing same things again and again, use alias and add it to .bashrc.
```
alias myfolder='cd /dartfs-hpc/qqqbbsdsdf/JC/'
```

If you add the code above in .bashrc, when you type myfolder, it performs 'cd /dartfs-hpc/qqqbbsdsdf/JC/'
keep in mind that you should reload .bashrc if you change or upload the file.
```bash
source ~/.bashrc
```





