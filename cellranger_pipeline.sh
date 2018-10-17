#1 cellranger download&installation

tar -xzvf cellranger-2.2.0.tar.gz    #Download and unpack the Cell Ranger file in any location
tar -xzvf refdata-cellranger-GRCh38-1.2.0.tar.gz    #Download and unpack any of the reference data files(GRCh38 in this example) in a convenient location
export PATH=cellranger-2.2.0path:$PATH    #Prepend the Cell Ranger directory to your $PATH

#2 verify installation

cellranger testrun --id=tiny    #This test can take up to 60 minutes on a sixteen-core workstation;Whether the test pipestance succeeds or fails, you will then see:Saving diagnostics to tiny/tiny.mri.tgz

#2.5 cellranger mkfastq可将Illumina sequencer's BCLs -> FASTQ files && 也可使用其他的FASTQ来源
#cellranger mkfastq --id=tiny-bcl \
#                   --run=/path/to/tiny_bcl \
#                   --samplesheet=cellranger-tiny-bcl-samplesheet-1.2.0.csv


#3 FASTQ文件命名及位置

#i) how to name FASTQs
#[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz

#Where Read Type is one of:

#I1: Sample index read (optional)
#R1: Read 1
#R2: Read 2


#ii) like this example
#Limb
#├── Limb-0710_S1_L001_R1_001.fastq.gz
#├── Limb-0710_S1_L001_R2_001.fastq.gz
#├── Limb-0710_S1_L002_R1_001.fastq.gz
#├── Limb-0710_S1_L002_R2_001.fastq.gz
#├── Limb-0710_S1_L003_R1_001.fastq.gz
#├── Limb-0710_S1_L003_R2_001.fastq.gz
#├── Limb-0710_S1_L004_R1_001.fastq.gz
#├── Limb-0710_S1_L004_R2_001.fastq.gz
#├── Limb-0710_S1_L005_R1_001.fastq.gz
#├── Limb-0710_S1_L005_R2_001.fastq.gz
#├── Limb-0710_S1_L006_R1_001.fastq.gz
#├── Limb-0710_S1_L006_R2_001.fastq.gz
#├── Limb-0710_S1_L007_R1_001.fastq.gz
#├── Limb-0710_S1_L007_R2_001.fastq.gz
#├── Limb-0710_S1_L008_R1_001.fastq.gz
#├── Limb-0710_S1_L008_R2_001.fastq.gz
#└── md5.md5

#iii) run_args
#All samples	--fastqs=/PATH/TO/PROJECT_FOLDER
#Process SAMPLENAME from all lanes	--fastqs=/PATH/TO/PROJECT_FOLDER --sample=SAMPLENAME
#Process SAMPLENAME from lane 1 only	--sample=SAMPLENAME --fastqs=/PATH/TO/PROJECT_FOLDER --lanes=1


#4 generate single-cell gene counts for a single library

#i) args
#--id	A unique run ID string: e.g. sample345
#--fastqs	Any folder containing fastq files;If you have multiple libraries for the sample, you will need to run cellranger count on them individually, and then combine them with cellranger aggr
#--sample	
#--transcriptome    Path to the Cell Ranger compatible transcriptome reference e.g.    For a human-only sample, use /opt/refdata-cellranger-GRCh38-1.2.0    For a human and mouse mixture sample, use /opt/refdata-cellranger-hg19-and-mm10-1.2.0
#--expect-cells    Expected number of recovered cells. Default: 3,000 cells.
#--r1-length    Hard-trim the input R1 sequence to this length
#--r2-length    Hard-trim the input R2 sequence to this length.
#--lanes    Lanes associated with this sample
#--localcores    Restricts cellranger to use specified number of cores to execute pipeline stages. By default, cellranger will use all of the cores available on your system.
#--localmem    Restricts cellranger to use specified amount of memory (in GB) to execute pipeline stages. By default, cellranger will use 90% of the memory available on your system. Please note that cellranger requires at least 16 GB of memory to run all pipeline stages.
#--indices    Sample indices associated with this sample. 

#tree
#GRCh38/(ref)     Limb/(fastq.gz)       Limb1016/(output)   nohup.out(run_infomation)

cellranger count --id=Limb1016 --transcriptome=./GRCh38 --fastqs=./Limb --sample=Limb-0710 --expect-cells=10000 --localcores=20

#The pipeline will create a new folder named with the sample ID you specified (e.g. /home/jdoe/runs/sample345) for its output. If this folder already exists, cellranger will assume it is an existing pipestance and attempt to resume running it.注意不能和输入参数中的数据文件夹重名，否则会报错"is not a pipestance directory"


#5 Output Files(outs/)

#web_summary.html	Run summary metrics and charts in HTML format
#metrics_summary.csv	Run summary metrics in CSV format
#possorted_genome_bam.bam	Reads aligned to the genome and transcriptome annotated with barcode information
#possorted_genome_bam.bam.bai	Index for possorted_genome_bam.bam
#filtered_gene_bc_matrices	Filtered gene-barcode matrices containing only cellular barcodes in MEX format
#filtered_gene_bc_matrices_h5.h5	Filtered gene-barcode matrices containing only cellular barcodes in HDF5 format
#raw_gene_bc_matrices	Unfiltered gene-barcode matrices containing all barcodes in MEX format
#raw_gene_bc_matrices_h5.h5	Unfiltered gene-barcode matrices containing all barcodes in HDF5 format
#analysis	Secondary analysis data including dimensionality reduction, cell clustering, and differential expression
#molecule_info.h5	Molecule-level information used by cellranger aggr to aggregate samples into larger datasets.
#cloupe.cloupe	Loupe Cell Browser visualization and analysis file
