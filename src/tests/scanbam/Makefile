FILES := \
	ex-pe-hard-leading.sorted.bam \
	ex-pe-hard-leading.sorted.bam.bai \
	ex-pe-hard-trailing.sorted.bam \
	ex-pe-hard-trailing.sorted.bam.bai \
	ex-pe-soft-neg.sorted.bam \
	ex-pe-soft-neg.sorted.bam.bai \
	ex-pe-soft-pos.sorted.bam \
	ex-pe-soft-pos.sorted.bam.bai \
	ex-pe-tlen.sorted.bam \
	ex-pe-tlen.sorted.bam.bai \
	ex-pe-tid.sorted.bam \
	ex-pe-tid.sorted.bam.bai \
	ex-pe-orient.sorted.bam \
	ex-pe-orient.sorted.bam.bai

.PHONY: all
all: $(FILES)

.PHONY: clean
clean:
	rm -f *.sorted.*

%.sorted.sam: %.sam
	samtools sort -O SAM -o $@ $<

%.bam: %.sam
	samtools view -O BAM -o $@ $<

%.bam.bai: %.bam
	samtools index $<