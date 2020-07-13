### Options
- **-i** path of the sam/bam file
- **-a** aligner [**BISMARK** | **BSMAP** | **UNKNOWN**]
- **-b** path of the bed file
- **-c** path of the cpg file
- **-r** region **chr1:2000-200000**
- **-o** output path

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

