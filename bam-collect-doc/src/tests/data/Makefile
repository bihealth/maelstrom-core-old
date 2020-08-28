FILES := \
	ex.sorted.bam \
	ex.sorted.bam.bai

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
