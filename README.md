The homepage: https://jiantaoshi.github.io/mHap/index.html

The PDF file of the paper: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab458/6305824

If you find mHapTools is helpful, please cite:

> ```latex
> @article{zhang2021dna,
>   title={The DNA methylation haplotype (mHap) format and mHapTools},
>   author={Zhang, Zhiqiang and Dan, Yuhao and Xu, Yaochen and Zhang, Jiarui and Zheng, Xiaoqi and Shi, Jiantao},
>   journal={Bioinformatics},
>   year={2021}
> }
> ```

### Build example

```bash
cd mHapTools
cd htslib-1.10.2
./configure --prefix=`pwd`
make
make install
cd ..
gcc -std=gnu99 -O2 -c convert.c bsread.c -I ./include -I ./htslib-1.10.2
g++ -O2 -o mhaptools haptk.cpp mhap.cpp merge.cpp beta.cpp summary.cpp utils.cpp convert.o bsread.o -I ./htslib-1.10.2/htslib -I ./include -L ./htslib-1.10.2/ -lhts -std=c++11
export LD_LIBRARY_PATH=`pwd`/htslib-1.10.2/lib
bash test/run_test.sh    # regression test of convert
```

If `configure` stops because the bzip2 or lzma headers are missing, use
``./configure --prefix=`pwd` --disable-bz2 --disable-lzma`` (only CRAM files
compressed with these codecs then cannot be read).

`convert` is written in C (`convert.c`, `bsread.c`) since version 0.11; the
other commands are C++. Any htslib >= 1.10 can replace the bundled copy.

### Commands

* **convert**

Convert SAM/BAM format file to mHap format file. It takes a coordinate-sorted Bisulfite-seq (or EM-seq, TAPS) BAM and CpGs position files as inputs to extract DNA methylation haplotypes. 

* **merge**

Merge multiple sorted mHap files, produce a single sorted mHap file.

* **beta**

Output summary of CpG site-level methylation from mHap files. It is similar to Bismark DNA methylation caller but uses mHap as inputs.

* **summary**

Computes the total number of reads, methylated CpG sites, total CpG sites, DNA methylation discordant reads,  methylated reads for given genomic regions or genome wide. 

### Details

#### convert

- **-i** input file, SAM/BAM/CRAM format, sorted by coordinate (`-` for stdin; `-r` and `-b` need an index).
- **-c** CpG file (chr, 1-based C position, ...), gz format; read per contig if tabix-indexed.
- **-r** region. **chr1:2000-200000**
- **-b** bed file, one query region per line (overlapping regions are merged, every read is used once).
- **-n** non-directional, do not group results by the direction of reads (strand `*`).
- **-m** sequencing mode. ( **TAPS** | **BS** (default)  )
- **-o** output filename. (default: out.mhap.gz; bgzipped and tabix-indexed (`tabix -b 2 -e 3`) if it ends in .gz)
- **-q** minimum mapping quality. (default: 10)
- **-F** skip reads with any of these flags. (default: 0xF04, unmapped, secondary, QC-fail, duplicate, supplementary)
- **-B** minimum base quality of a CpG call. (default: 0)
- **-L** maximum distance between the mates of a fragment, in bp. (default: 1000)
- **--max-unconv** drop incompletely converted reads: at least 3 C outside CpGs, making up more than FLOAT of the C and T bases outside CpGs (G and A on the bottom strand). (default: 0.2; -1 = no filter; not applied with `-m TAPS`)
- **--max-ch** drop reads with more than INT unconverted cytosines outside CpGs. Bases filled in at fragment ends during library preparation and genuine non-CpG methylation are counted too, so this removes many good reads (37% on a BISCUIT sample with `--max-ch 2`). (default: -1, no filter)
- **--split** cut a fragment at a sequenced CpG without a call instead of dropping the fragment.
- **--qname** one line per record with the read name in column 7 (no collapsing).
- **--no-index** do not write the tabix index.
- **-T** reference FASTA (CRAM input).
- **-@** extra threads for decompression and compression.

How reads are converted (since 0.11):

1. Reads are laid out on the reference with their CIGAR: soft clips and insertions are skipped, deletions leave a gap.
2. The bisulfite strand comes from the aligner's tag when present (Bismark `XG`, bwa-meth / BISCUIT `YD`, BSMAP `ZS`), otherwise from the flags (bottom = read 1 reverse or read 2 forward). Non-directional libraries (CTOT/CTOB reads) need the tags.
3. Reads are filtered by flags (`-F`), mapping quality (`-q`) and conversion (`--max-unconv`). The conversion filter needs no reference: outside CpGs, C makes up almost none of the C and T bases of a converted top-strand read and about 40% of those of an unconverted one (G and A on the bottom strand).
4. Mates are merged into one fragment. A CpG is called only if its C and G both lie on the read and the read shows the CpG context (top strand: C/T followed by G; bottom strand: C followed by G/A); a CpG on which the mates disagree has no call. These are the rules of [wgbs_tools](https://github.com/nloyfer/wgbs_tools) `patter`.
5. An mHap haplotype covers consecutive CpGs, so the calls are cut where a CpG has no call. CpGs between two mates that do not overlap were not sequenced, and the fragment gives two records. A sequenced CpG without a call drops the fragment (with `--split`, the fragment is cut there instead; cutting turns one molecule into several records and lowers read-level statistics such as the M-score).
6. Records are sorted, identical records collapsed into the count, and written bgzipped with a tabix index.

Changes from 0.10: soft clips and indels were ignored (bases were read by offset from the start of the sequence), only the flags gave the strand and improper pairs were dropped, mates were merged only when their CpG spans overlapped, any non-C/T (non-G/A) base at a CpG dropped the read, only Bismark reads were checked for conversion (any `X`/`H`/`U` in `XM`), there was no MAPQ filter, a CpG file whose contig names differed from the BAM (`chr1` / `1`) gave an empty result without error, and no index was written. On simulated data with 0.5% sequencing errors, 0.11 calls 99.8% of CpGs correctly on reads with indels or soft clips (0.10: 97.7% paired-end, 92.7% single-end) and on non-directional CTOT/CTOB reads with strand tags (0.10: 40.8%), and 99.8% when 5% of the fragments are unconverted (0.10, on a BAM without Bismark's `XM` tag: 96.4%). On a BISCUIT paired-end sample (2.0M reads, hg38), 0.11 takes 13.9 s and 63 MB (0.10: 30.9 s and 293 MB). The conversion code and its validation against wgbs_tools are shared with [mHapASM](https://github.com/JiantaoShi/mHapASM).

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
tabix -b 2 -e 3 file.mhap.gz
```

* **-c** CpG file, gz format.
* **-b** bed file of query regions.
* **-r** query region, e.g. chr1:2000-20000.
* **-o** output fiename. (summary.txt | summary_genome_wide.txt)
* **-s** group results by the direction of mHap reads.
* **-g** get genome-wide result.



