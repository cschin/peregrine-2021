# Peregrine-2021: A faster and minimum genome assembler for long-reads with good enough accuracy with the Rust language

## System requirement: 

A modern Linux workstation or compute node with enough disk, CPUs, and RAM. It is better to have a good number of CPUs (my testing system has 20 cores) and a good amount of RAM (~ total 1.5x of the reads data set). For example, for 100G sequences, it is probably good to have at least 150G RAM. A smaller amount, e.g., 32G, works, but you will need some manual setup for effective computation.
## Some Ballpark Performance Summary 

With a proper hardware (e.g. ~1Tb RAM), Peregrine-2021 had successful assembled a total 30G diploid genome (2n = 30G) with a contig N50 = 55.2Mb for a large diploid genome. (For who might want to know more details, I ran it as a unpaid service so I don't have much other infomation. [Link is the the graph of the assembly](https://twitter.com/infoecho/status/1330617986185457669?s=20&t=FlHjuWCHslvjxVdyZpU1gQ).)

For a typical human-size assembly, a much cheaper compute instance with from 128G to 512G RAM can work well. (see this blog [Accelerating genome assembly with AWS Graviton2](https://aws.amazon.com/blogs/publicsector/accelerating-genome-assembly-aws-graviton2/)), and it takes only 2 to 3 hours wall o'clock time to get an assembly. (We also provide a "fast" mode eliminating one error correct stage for perfect reads as input.)

On the accuracy side, here, we only have some rough numbers from the earlier version compared to other assemblers. Ironically, it could take more effort and resources to do a comprehensive benchmark for a publication than writing the assembler itself. Unfortunately, benchmarking work is a luxury that I currently do not have. However, suppose if anyone is interested in taking a shot running the code fine-tuning the parameters to evaluate the results correctly, in that case, I might help a bit from time to time. 

## Usage:

1.	Put the paths to the sequence read files (fasta / fasta.gz / fastq / fastq.gz, compressed with the standard gzip but not bgzip) in a file, e.g. `reads.lst`, so the Peregrin-2021 assembler can find the read data. For example, this shows the content of a `reads.lst` file:
```
$ cat seq.lst 
/wd/CHM13_20k/m64062_190803_042216.fastq.gz
/wd/CHM13_20k/m64062_190804_172951.fastq.gz
/wd/CHM13_20k/m64062_190806_063919.fastq.gz
/wd/CHM13_20k/m64062_190807_194840.fastq.gz
```

2.	Make sure your have enough disk (preferably SSD storage or high performance network filesystem) for a working directory.  Let’s call the working directory `asm_out`.

3.	Execute: `pg_asm reads.lst asm_out` from command line / shell, some potentially useful intermediate files and the assembled contigs will be in the directory `asm_out/`

```
$ pg_asm seq.lst asm_out >& log &
```

4.	There are a number of options that you can try to tune for optimizing the assembly results. Here is the full usage information of `pg_asm`.

```
❯ pg_asm --help
pg_asm peregrine-2021 0.4.9 (arm_config:58e666e+, release build, linux [x86_64] [rustc 1.58.0 (02072b482 2022-01-11)])
Jason Chin <jason@omnibio.ai>

Peregrine-2021 genome assembler,
LICENSE: http://creativecommons.org/licenses/by-nc-sa/4.0/

USAGE:
    pg_asm [FLAGS] [OPTIONS] <input_reads> <work_dir> [ARGS]

FLAGS:
        --fast          run the assembler in the fast mode
    -h, --help          Prints help information
        --keep          keep intermediate files
        --no_resolve    disable resolving repeats / dups at the end
    -V, --version       Prints version information

OPTIONS:
    -b, --bestn <bestn>              number of best overlaps for initial graph [default: 6]
    -k <k>                           Kmer size [default: 56]
    -l <layout_method>               layout version [default: 2]
        --log <log>                  log level: DBBUG or INFO (default)
    -c, --min_ec_cov <min_ec_cov>    Minimum error coverage [default: 1]
    -r <r>                           Reduction factor [default: 6]
    -t, --tol <tol>                  Alignment tolerance [default: 0.01]
    -w <w>                           Window size [default: 80]

ARGS:
    <input_reads>    Path to a file that contains the list of reads in .fa .fa.gz .fastq or fastq.gz formats
    <work_dir>       The path to a work directory for intermediate files and the results
    <NTHREADS>       Number of threads
    <NCHUNKS>        Number of partition
```

You can reduce the `Reduction factor` and `Window Size` to increase the sensitivity to 
detect overlaps. I found `-r 4 -w 64` may be better for human assembly than the default
parameters. 

## Example Results: (`v0.2.0, main:6d91294, release build, linux [x86_64]`)

### Input data set: Human CHM13 dataset

- 5,567,153 reads, average length = 18,028 bp, about 30x human genome size

- Systme: 20 core Intel(R) Xeon(R) CPU E5-2630 v4 @ 2.20GHz, 512G ram


### Run Peregrine-2021 in Fast-mode (without read level error correction)
 

#### CPU time for the fast-mode:
```
usr: 152025s = 42.5 cpu hours, 

sys: 900s = 0.25 cpu hours, 

wall clock time = 2:45:39
``` 
(In contrast, it takes [HiCanu](https://www.biorxiv.org/content/10.1101/2020.03.14.992248v3) > 4000 CPU hour to get an assembly)

#### Assembly Summary Statistics with the fast-mode:

```
total size: 3,034,243,471
max size: 201,005,110
N50 size: 81,361,265
N90 size: 17,301,175
Number of Contigs: 237
Number of Contigs > 100kb: 151
```

#### CHM13 BAC evaluation result with the fast-mode setting:

```
******************* BAC SUMMARY ******************
 TOTAL    : 341
 BP       : 51532183
************** Statistics for: _asm_p_ctg-2.fa ****************
BACs closed: 327 (95.8944); BACs attempted: 338 %good = 96.7456; BASES 49454907 (95.969)
Median:          99.9844
MedianQV:        38.06875
Mean:            99.93617
MeanQV:          31.94969
***** STATS IGNORING INDELS ********************
Median:          100
MedianQV:        Inf
Mean:            99.99177
MeanQV:          40.84457
**********************************************
```
(This [HiCanu paper preprint](https://www.biorxiv.org/content/10.1101/2020.03.14.992248v3) reported resolving only 326 out of the 341 BAC BACs.)

#### CHM13 BAC evaluation result with the fast-mode setting but without final contig level consensus:
(This pretty much just reflects the quality of the input reads.)

```
******************* BAC SUMMARY ******************
 TOTAL    : 341
 BP       : 51532183
************** Statistics for: asm-p_cns-nocns_fast.fa ****************
BACs closed: 327 (95.8944); BACs attempted: 339 %good = 96.4602; BASES 49454907 (95.969)
Median:          99.7163 
MedianQV:        25.47141 
Mean:            99.67459 
MeanQV:          24.87566 
***** STATS IGNORING INDELS ********************
Median:          99.9782 
MedianQV:        36.61544 
Mean:            99.97005 
MeanQV:          35.23557 
**********************************************
```

### Run Peregrine-2021 in Standard-mode (one round of the read level error correction)

#### CPU time for the standard setting which involves one round of read level error correction: 

```
usr: 355257s = 98.7 cpu hours,  

sys: 2052s = 0.6 cpu hours, 
  
wall clock time = 5:47:23
```
(In contrast, it takes [HiCanu](https://www.biorxiv.org/content/10.1101/2020.03.14.992248v3) > 4000 CPU hour to get an assembly)

#### Assembly Summary Statistics of the standard Setting:

```
total size: 3,039,592,838
max size: 142,204,433
N50 size: 83,143,579
N90 size: 16,250,509
Number of Contigs: 227
Number of Contigs > 100kb: 157
```

#### CHM13 BAC evaluation result:

```
******************* BAC SUMMARY ******************
 TOTAL    : 341
 BP       : 51532183
************** Statistics for: _asm_p_ctg.fa ****************
BACs closed: 330 (96.7742); BACs attempted: 341 %good = 96.7742; BASES 49874385 (96.783)
Median:          99.9911
MedianQV:        40.5061
Mean:            99.94317
MeanQV:          32.45434
***** STATS IGNORING INDELS ********************
Median:          100
MedianQV:        Inf
Mean:            99.99405
MeanQV:          42.25218
**********************************************
```
#### CHM13 BAC evaluation result without final contig level consensus:

```
******************* BAC SUMMARY ******************
 TOTAL    : 341
 BP       : 51532183
************** Statistics for: asm-p_cns-nocns.fa ****************
BACs closed: 330 (96.7742); BACs attempted: 341 %good = 96.7742; BASES 49874385 (96.783)
Median:          99.9832
MedianQV:        37.74691
Mean:            99.93411
MeanQV:          31.81188
***** STATS IGNORING INDELS ********************
Median:          99.9994
MedianQV:        52.21849
Mean:            99.9927
MeanQV:          41.36713
**********************************************
```

  
## Some FAQs: 

* Q: Why do you write a new assembler? There are already many others (FALCON, HifiAsm, Flye, Shasta, HiCanu, etc.).

    A: We demonstrated that it was possible to use Spare Hierarchical MiniMiER (SHIMMER) to assemble a human-size genome with 100 wall clock minutes with the Peregrine genome assembler. A standard approach filters out high duplicated sequences to get the significant part of a gnome assembled. And, most assemblers adapt some repeat filtering schemes to make the run time for genome assembly acceptable while keeping the most helpful information of a genome. 

    However, while it was reasonable to do a repeat-suppressing assembly, the narrative about what an assembler should do is changing. Genomic researchers may want to get more about the repeat even if it needs additional compute power/energy ( > 20x for a human genome compared to a Peregrine run) to get the repetitive sequences in a genome. I think it is worth showing the same technique that we used in the original Peregrine-2021 can also get the high repeat content assembled with only moderate increases of computation cost/energy.

    Another reason for this Peregrine-2021 assembler is that we were not happy with the C / Python hybrid approach used in the original Peregrine assembler. While C / Python combination is very efficient for rapid development, it has too many caveats. It could be interesting to learn something new as well. Following Richard Feynman's wisdom, "What I cannot create, I don't understand." To better understand how to handle repeats and understand the Rust programming language, I created Peregrine-2021 from late-2020 to mid-2021. I want to push Peregrine-2021 to the next stage to apply it to genomics research work; unfortunately, I realized it is too demanding to take this as a hobby or a side-gig. Given that I won't be able to push it too far by myself, it might still be helpful for others who are interested in using this. I finally decide to release the code for non-commercial use.

* Q: I can't run it, can you help?

    A: it depends. If it is something that is straight forward, I am happy to help. If it is more involved, then I simply can not do it as I only have limited resources. 

* Q: The results are bad, do I do something wrong?
  
    A: Possiblely. Like all other assemblers, Peregrine-2021 is design with some sepcifications in mind and test accordingly. Depending on the input data and the parameters used to run the assemblers, one might be able to improve the results. Unfortunately, there is no universal simple answer for how to investigate and improve now. It can be some trivial mistakes or very intensive investigation that related to inital sequence methods or even the genome biology itself. 

* Q: I use it and generate some good results for publication, how do I cite it?

    A: While Peregrine-2021 is mostly based on the original ideas described in the preprint [Human Genome Assembly in 100 Minutes](https://www.biorxiv.org/content/10.1101/705616v1), it is a different codebase. I currently have no plan to write up another manuscript specifically for peregrine-2021. Yes, please cite the preprint if you find it is useful.

* Q: Is any published work using Peregrine-2021?

    A: Check [this](https://scholar.google.com/scholar?cites=12093836825307648052). Most citations are for benchmarking purposes (of the older version), and some are for the related ideas. There are a couple of papers using it for generating results. As this is mostly a "hobby" project which lacks resources for promoting its usage, I am grateful that the Chief Scientific Officier of Medicinal Genomics, Keven McKernan, provides data sets for me to test it and push the results into papers.

    * [A draft reference assembly of the Psilocybe cubensis genome](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC8220353/) 

    * [Cannabis Genome: Jamaican Lion strain ](https://www.medicinalgenomics.com/jamaican-lion-data-release/), Note Peregrine-2021 assembled thie genome in an old Mac Pro with only 64G RAM.

* Q: I don't have a large memory machine, how do I run Peregrine-2021?

    A: For efficiency, it will be great to put all sequence data and the smaller index data in the RAM. Currently, my suggestion is that the RAM should be about 1.5x of the total sequence data. However, the data is accessed through the memory-mapped file (MMAP) mechanism, and the chunking machinery for parallelization can help if one does not have a big memory machine. The code can access the data from the disk through the MMAP file. In such a case, high efficient NVME SSD will help. I had successfully assembled a human genome using a 32G RAM machine. However, I won't recommend that is the right way to go, given that a medium-size RAM machine is relatively cheap to rent now.     

* Q: Can you write better Rust code?

    A: I guess I could, but I was literally learning Rust and developing new algorithms at the same time. Rust is a big language, there is still a lot to learn

* Q: Why do you choose CC BY-NC-SA 4.0 license to release the code?

    A: oh, that will be a fire side or beer hour conversation.
## How does Peregrine-2021 work?
	
Peregrine-2021 is still just another OLC assembler. 

1. It uses the SHIMMER index for finding overlap candidates as described in the Peregrine preprint but it is more aggressive for overlapping repetitive reads rather than filtering them out.

2. The overlaper performs analysis to identify read overlaps within the same haplotype. We did not use a de-Bruijn graph approach for this. We think our method is likely more computational efficient than using a de-Bruijn graph to separate haplotypes.

3. It adapts the techniques using partial homopolymer compression for separating the reads from different haplotype.

## Build For X86_64

1. Check [Rust Installation](https://www.rust-lang.org/tools/install)

2. Run [`cargo install`](https://doc.rust-lang.org/cargo/commands/cargo-install.html) or [`cargo build`](https://doc.rust-lang.org/cargo/commands/cargo-build.html). Make sure you set up the environment variable `PATH` to the directory of the built binaries or you can run the excutable `pg_asm` with full path. If you use `cargo build`, make sure you compile it with the `--release` option for optimization.

3. `pg_asm` will run the assembly pipeline end-to-end. If it fails, it does not re-use the existing data when one runs `pg_asm` again. The assemblers is much faster than other assemblers, so it is less important to re-use intermidate data. That has been said, the built will contains executables (e.g. `pg_build_idx`, `pg_ovlp`, etc.) for each assembly steps which one can chain them together with their favorite workflow engine for re-using and re-starting an assembly pipeline.

4. A Dockerfile is provided for creating a Docker image. It also provide information to build the assembler from a clean environment.

5. To compile in aarch64, it will need some configuration changes to get the best performance. The memory alloctor package needs to be patched for aarch64. See [https://github.com/cschin/mimalloc_rust/tree/aarch64_build] (https://github.com/cschin/mimalloc_rust/tree/aarch64_build).

## Other utility command line tools

```
pg_build_sdb  # convert fasta/fastq/fasta.gz/fastq.gz read data into a simple binary data base for the assembler to fetch the reads. 

pg_build_idx  # build the SHIMMER index from the reads for overlapping

pg_build_sdb  # build the sequence database

pg_dedup      # perform all contigs to all contigs alignment to remove duplicates 

pg_dp_graph   # take overlap data file as input to generate the layout file using a polyploid aware layout algorithm

pg_getreads   # generate fasta file for a subset of reads from the sequence database

pg_graph      # (obsoleted) convert the overlap information between the reads into an assembly group

pg_layout     # convert the assembly graph to paths and generate the contig fasta file

pg_ovlp_ec    # perform error correction from the haplotype specific overlaps

pg_ovlp       # generate haplotype specific overlaps between the reads

pg_resolve    # this tool aligns all contigs to themselve to identify haplotype-related contigs
```

--
Jason Chin (twitter: @infoecho)

first version: Nov. 16, 2020

current version: Fab. 5, 2022


