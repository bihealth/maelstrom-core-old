FILES := \
	ex-duplicates.bam \
	ex-duplicates2.bam

.PHONY: default
default: $(FILES)

%.bam: %.sam
	samtools view --no-PG -O BAM -o $@ $<
