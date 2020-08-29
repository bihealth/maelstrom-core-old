LANG := C
SHELL := bash

#	ex-delly-snvs.vcf.gz
#	ex-delly-snvs.vcf.gz.tbi
FILES := \
	ex-delly-doc.vcf.gz \
	ex-delly-doc.vcf.gz.tbi \
	ex-delly-pesr.tsv.gz \
	ex-delly-pesr.tsv.gz.tbi \
	ex-delly-svs.vcf.gz \
	ex-delly-svs.vcf.gz.tbi

.PHONY: default
default: $(FILES)

%.tsv.gz: %.tsv
	( \
		grep '^#' $<; \
		grep -v '^#' $< | sort -k1,1 -k2,2n -k3,3n; \
	) \
	| bgzip -c >$@

%.tsv.gz.tbi: %.tsv.gz
	tabix -p bed -f $<

%.vcf.gz: %.vcf
	( \
		grep '^#' $<; \
		grep -v '^#' $< | sort -k1,1 -k2,2n -k3,3n; \
	) \
	| bgzip -c >$@

%.vcf.gz.tbi: %.vcf.gz
	tabix -f $<
