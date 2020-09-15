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

```c++
//Convert from mhap file to bgzip file
cat file.mhap | sort -k1,1 -k2,2n | bgzip > file.mhap.gz
//Convert from bgzip file to tabix index file
tabix -b 2 -e 3 -p bed file.mhap.gz
```

* **-c** CpG file, gz format.
* **-b** bed file of query regions.
* **-r** query region, e.g. chr1:2000-20000.
* **-o** output path. (summary.txt | summary_genome_wide.txt)
* **-s** if specified, the results are grouped by the direction of mhap reads.
* **-g** genome-wide result.

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

