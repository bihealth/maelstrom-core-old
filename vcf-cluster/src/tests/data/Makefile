FILES := \
	ex-delly-1.vcf.gz \
	ex-delly-1.vcf.gz.tbi \
	ex-delly-2.vcf.gz \
	ex-delly-2.vcf.gz.tbi \
	ex-cnmops-1.vcf.gz \
	ex-cnmops-1.vcf.gz.tbi \
	ex-cnmops-2.vcf.gz \
	ex-cnmops-2.vcf.gz.tbi

.PHONY: default
default: $(FILES)

.PHONY: clean
clean:
	rm -f $(FILES)

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.tbi: %.vcf.gz
	tabix -f $<
