###Options
- **-i** path of the bam file
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

