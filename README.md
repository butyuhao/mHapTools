### Command

* convert 

Convert SAM/BAM file to mhap file.

* merge

Merge two mHap files.

### Options

#### convert

- **-i** input file, SAM/BAM format, should be sorted by samtools.
- **-n** non-directional
- **-b** bed file of query regions.
- **-c** CpG file, gz format.
- **-r** region. **chr1:2000-200000**
- **-o** output path. (default: out.mhap)

#### merge

* **-i** input file, mhap format.
* **-c** CpG file, gz format.
* **-o** output path. (default: out.mhap)

#### beta

* **-i** input file, mhap format
* **-c** CpG file, gz format.
* **-o** output path. (default: beta.txt)
* **-s** if specified, the results are grouped by the direction of mhap reads.
* **-b** bed file of query regions.

#### summary

* **-i** input file, mhap.gz format or mhap format (if opt -g is specified)
* **-c** CpG file, gz format.
* **-b** bed file of query regions.
* **-r** query region, e.g. chr1:2000-20000.
* **-o** output path. (summary.txt | summary_genome_wide.txt)
* **-s** if specified, the results are grouped by the direction of mhap reads.
* **-g** genome-wide result.

### Build example

```bash
cd haptools
cd htslib-1.10.2
./configure --prefix=`pwd`
make
make install
cd ..
g++ -o haptools  convert.cpp mhap.cpp haptk.cpp merge.cpp beta.cpp summary.cpp -I ./htslib-1.10.2/htslib -I ./include  -L ./htslib-1.10.2/ -lhts -std=c++11
export LD_LIBRARY_PATH=`pwd`/htslib-1.10.2/lib
```

### Index File

Index files are needed when you use `-r`  to specify a region or when you use `-b` to specify a bed file. Suppose you have a bam file named `XXX.bam` and a cpg file named `YYY.gz`, the names of corresponding index files should be `XXX.bam.bai` and `YYY.gz.tbi`. Index files should be placed in the same directory as indexed files.

