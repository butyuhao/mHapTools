### Command

* convert 

Convert SAM/BAM file to hap file.

* merge

Merge two hap files.

### Options

#### convert

- **-i** input file, SAM/BAM format, should be sorted by samtools.
- **-a** aligner. [**BISMARK** | **BSMAP** | **UNKNOWN**(default)]
- **-b** bed file of query regions.
- **-c** CpG file, gz format.
- **-r** region. **chr1:2000-200000**
- **-o** output path. (default: out.hap)

#### merge

* **-i** input file, hap format.
* **-c** CpG file, gz format.
* **-o** output path. (default: out.hap)

#### beta

* -i input file, hap format
* -c CpG file, gz format.
* -o output path. (default: beta.txt)
* -s if specified, the results are grouped by the direction of hap reads.

### Build example

```bash
mkdir build
cd build
cmake ..
make
./haptools convert -i /XXX/XXX.bam -c /XXX/XXX.gz -r X:2000-20000 -a BISMARK -o out.hap 
```

### Index File

Index files are needed when you use `-r`  to specify a region or when you use `-b` to specify a bed file. Suppose you have a bam file named `XXX.bam` and a cpg file named `YYY.gz`, the names of corresponding index files should be `XXX.bam.bai` and `YYY.gz.tbi`. Index files should be placed in the same directory as indexed files.

