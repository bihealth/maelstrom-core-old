FILES := \
	ex-delly-filter.vcf.gz \
	ex-delly-filter.vcf.gz.tbi \
	ex-delly-nofilter.vcf.gz \
	ex-delly-nofilter.vcf.gz.tbi

.PHONY: default
default: $(FILES)

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.tbi: %.vcf.gz
	tabix -f $<
