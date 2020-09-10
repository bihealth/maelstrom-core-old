FILES := \
	ex-cnmops-filter.vcf.gz \
	ex-cnmops-filter.vcf.gz.tbi \
	ex-cnmops-nofilter.vcf.gz \
	ex-cnmops-nofilter.vcf.gz.tbi \
	ex-delly-filter.vcf.gz \
	ex-delly-filter.vcf.gz.tbi \
	ex-delly-nofilter.vcf.gz \
	ex-delly-nofilter.vcf.gz.tbi \
	ex-manta-filter.vcf.gz \
	ex-manta-filter.vcf.gz.tbi \
	ex-manta-nofilter.vcf.gz \
	ex-manta-nofilter.vcf.gz.tbi

.PHONY: default
default: $(FILES)

.PHONY: clean
clean:
	rm -f $(FILES)

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.tbi: %.vcf.gz
	tabix -p vcf -f $<
