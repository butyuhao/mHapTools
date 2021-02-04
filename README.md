### Build example

```bash
cd mHapTools
cd htslib-1.10.2
./configure --prefix=`pwd`
make
make install
cd ..
g++ -o mhaptools  haptk.cpp convert.cpp mhap.cpp merge.cpp beta.cpp summary.cpp utils.cpp -I ./htslib-1.10.2/htslib -I ./include  -L ./htslib-1.10.2/ -lhts -std=c++11
export LD_LIBRARY_PATH=`pwd`/htslib-1.10.2/lib
```

### Commands

* **convert** 

  Convert SAM/BAM format file to mHap format file. It takes an indexed Bisulfite-seq BAM and CpGs position files as inputs to extract DNA methylation haplotypes. 

  

* **merge**

  Merge multiple sorted mHap files, produce a single sorted mHap file.

  

* **beta**

  Output summary of CpG site-level methylation from mHap files. It is similar to Bismark DNA methylation caller but uses mHap as inputs.

  

* **summary**

  Computes the total number of reads, methylated CpG sites, total CpG sites, DNA methylation discordant reads for given genomic regions or genome wide. 

  

### Details

#### convert

- **-i** input file, SAM/BAM format, should be sorted by samtools.
- **-n** non-directional, do not group results by the direction of reads.
- **-b** bed file, one query region per line.
- **-c** CpG file, gz format.
- **-r** region. **chr1:2000-200000**
- **-o** output filename. (default: out.mhap.gz)

#### merge

* **-i** input file, multiple .mhap.gz files to merge.
* **-c** CpG file, gz format.
* **-o** output filename. (default: merge.mhap.gz)

#### beta

* **-i** input file, .mhap.gz format.
* **-c** CpG file, gz format.
* **-o** output filename. (default: beta.txt)
* **-s** group results by the direction of mHap reads.
* **-b** bed file, one query region per line.

#### summary

* **-i** input file, mhap.gz format.

```c++
//Generate index for .mhap.gz file
tabix -b 2 -e 3 -p bed file.mhap.gz
```

* **-c** CpG file, gz format.
* **-b** bed file of query regions.
* **-r** query region, e.g. chr1:2000-20000.
* **-o** output fiename. (summary.txt | summary_genome_wide.txt)
* **-s** group results by the direction of mHap reads.
* **-g** get genome-wide result.
