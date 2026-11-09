# Changelog

## 1.0.0 (2026-11-09)


### Features

* **config:** provide portable migration defaults and Biowulf profile ([ff546ad](https://www.github.com/shukwong/TriosCompass_v2/commit/ff546add1cbd519fd916b0a29e6f28dda458675d)), closes [#6](https://www.github.com/shukwong/TriosCompass_v2/issues/6)
* **qc:** add Telomerecat and filtered BAM exports ([02620b3](https://www.github.com/shukwong/TriosCompass_v2/commit/02620b326b340a8432d542683d9ed9f6d04ffe6f)), closes [#2](https://www.github.com/shukwong/TriosCompass_v2/issues/2)
* **README.md:** Add new section about the improvements to call DNMs with Strelka ([b266390](https://www.github.com/shukwong/TriosCompass_v2/commit/b26639016e98b47567b5cc5df8a488e56bc5dc9f))
* **reports:** add offline interactive trio SV reports ([a842aea](https://www.github.com/shukwong/TriosCompass_v2/commit/a842aea9f3bad67368659ffd92105790a9cbb90c)), closes [#3](https://www.github.com/shukwong/TriosCompass_v2/issues/3)
* **rufus:** add independent trio calling and record summaries ([d0c76de](https://www.github.com/shukwong/TriosCompass_v2/commit/d0c76deda9854302dd5079fb5f73d0f141873052)), closes [#1](https://www.github.com/shukwong/TriosCompass_v2/issues/1)
* **scripts/prepare_yml_from_vcf.pl:** add new perl script to generate yaml file for the use of igv_snapshot_maker ([4002e82](https://www.github.com/shukwong/TriosCompass_v2/commit/4002e82ed26bf2acc88200f275329adeae10de09))
* **Snakefile:** new improvement for Strelka DNM calls ([36a6726](https://www.github.com/shukwong/TriosCompass_v2/commit/36a6726943af9f0473edfc3150e1121bc2c994ee))
* **validation:** add portable full-trio output audit ([05bc19f](https://www.github.com/shukwong/TriosCompass_v2/commit/05bc19fab5a6659f6724369b3f9d2d14cfd3e59b))


### Bug Fixes

* add memory requirements to GATK CombineGVCFs step ([c24b2fe](https://www.github.com/shukwong/TriosCompass_v2/commit/c24b2fea46ee670a9952f79adff2e3c1f22d7360))
* **ci:** grant release automation scoped write permissions ([1480d0c](https://www.github.com/shukwong/TriosCompass_v2/commit/1480d0c44238c1f31f0e3b279e852e5b3c91acca))
* **env:** retain pkg_resources for Snakemake 7 ([9dbd8ec](https://www.github.com/shukwong/TriosCompass_v2/commit/9dbd8ec1acdf042bc0ce4ff939de5b5528759e1d)), closes [#6](https://www.github.com/shukwong/TriosCompass_v2/issues/6)
* **workflow/profiles/slurm/config.yaml:** change memory for glnexus_dv and aggregate_phase ([55605b1](https://www.github.com/shukwong/TriosCompass_v2/commit/55605b1325f6716618e359d23df4dbba46068c50))
* **workflow/profiles/slurm/config.yaml:** for the slurm setting used in CCAD (some of them might not be compatible at Biowulf) ([b6733d9](https://www.github.com/shukwong/TriosCompass_v2/commit/b6733d99f86aa2ec938ffc542bda033950472322))
* **workflow/report/dnm.rst:** fix the typo in rst ([957f593](https://www.github.com/shukwong/TriosCompass_v2/commit/957f593de2078334a088d51314b946f489db4956))
* **workflow/rules/dnSTR.smk:** dnSTR summary for each trio (not family) ([8bc1426](https://www.github.com/shukwong/TriosCompass_v2/commit/8bc1426e3d43bb12172aec7f767291296c39ff45))
* **workflow/rules/dnSV_manta.smk:** replace manta docker image to fix: taskExitCode -11 ([0513d22](https://www.github.com/shukwong/TriosCompass_v2/commit/0513d22322b1ec93d46708dbb6a326c42c330bd9))
* **workflow/rules/hipstr.smk:** exit due to zgrep and also fix one typo ([2d1746e](https://www.github.com/shukwong/TriosCompass_v2/commit/2d1746e427acf6de32472e754b65293358d3e5e8))
* **workflow/rules/multiqc.smk:** use docker://multiqc/multiqc:v1.33 to fix issue in the use of multiqc conda or wrapper ([6b5659d](https://www.github.com/shukwong/TriosCompass_v2/commit/6b5659d6c3bba521afa5c45e323934bdce4706da))
* **workflow/rules/pedigree.smk:** sorted pedigree file in the order: father,mother,child ([b92463c](https://www.github.com/shukwong/TriosCompass_v2/commit/b92463c09dd9f9007ac520e0db9e3dec6ae735af))
* **workflow:** address runtime and output-audit review findings ([4684108](https://www.github.com/shukwong/TriosCompass_v2/commit/4684108c639115e122ac5e8f7897dc9e35817906))
* **workflow:** containerize tools and repair GIAB recovery failures ([f5013ff](https://www.github.com/shukwong/TriosCompass_v2/commit/f5013ffaa47a31abf3728909f3e058d1eba6049f))
* **workflow:** repair GATK memory and HipSTR output staging ([364b9ca](https://www.github.com/shukwong/TriosCompass_v2/commit/364b9cac20edd003d841bebef2ab95199de38515))
