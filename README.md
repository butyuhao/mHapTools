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

* convert 

Convert SAM/BAM file to mHap file.

* merge

Merge two mHap files.

* beta

Get the percentage of methylated bases, number of methylated bases and number of unmethylated bases across reads on each CpG position.

* summary

Get number of reads, number of methylated bases, number of total bases, K4plus and nDR on each CpG position or within region(s).

### Details

#### convert

- **-i** input file, SAM/BAM format, should be sorted by samtools.
- **-n** non-directional
- **-b** bed file of query regions.
- **-c** CpG file, gz format.
- **-r** region. **chr1:2000-200000**
- **-o** output path. (default: out.mhap.gz)

#### merge

* **-i** input file, .mhap.gz format.
* **-c** CpG file, gz format.
* **-o** output path. (default: merge.mhap.gz)

#### beta

* **-i** input file, .mhap.gz format.
* **-c** CpG file, gz format.
* **-o** output path. (default: beta.txt)
* **-s** if specified, the results are grouped by the direction of mHap reads.
* **-b** bed file of query regions.

#### summary

* **-i** input file, mhap.gz format.

```c++
//Generate index for .mhap.gz file
tabix -b 2 -e 3 -p bed file.mhap.gz
```

* **-c** CpG file, gz format.
* **-b** bed file of query regions.
* **-r** query region, e.g. chr1:2000-20000.
* **-o** output path. (summary.txt | summary_genome_wide.txt)
* **-s** if specified, the results are grouped by the direction of mHap reads.
* **-g** genome-wide result.
