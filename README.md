# snakemake-pgx

  <p align="center">
    A Snakemake workflow for extracting PGx information from WES data
    <br />
    <a href="https://github.com/jucosie/snakemake-pgx"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/jucosie/snakemake-pgx/issues">Report Bug</a>
    ·
    <a href="https://github.com/jucosie/snakemake-pgx/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#obtain a copy of this workflow">Obtain a copy of this workflow</a></li>
        <li><a href="#configure the workflow">Configure the workflow</a></li>
        <li><a href="#install snakemake">Install Snakemake</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#example">Example</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This tool was developed with the aim of providing a unified solution to extract relevant PGx information from a set of specific genes in WES data. It brings together tools for genotyping HLAs (``Optitype``) and haplotyping of pharmacogenes (``Aldy``), together with position-specific variant calling (``GATK4``) and coverage data (``Mosdepth``). 

![Pipeline](https://github.com/jucosie/snakemake-pgx/blob/main/pipeline.png){width=70%,height=70%}
### Built With

* Python


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

We have tried to include all the tools necessary for the pipeline to work as wrappers, or via Conda envirotnments. Even so, it is necessary to install locally: 
* [GATK4 (2.1.0)](https://github.com/broadinstitute/gatk). GATK4 requires Java 8 to run. The path to the GATK wrapper script must be in the user's ``PATH`` variable.
* [Picard](https://github.com/broadinstitute/picard/releases). The path to the ``.jar`` file has to be specified in the configuration file. 

In addition, the user needs to download the following reference files and indicate their path in the configuration file: 
* BWA index of the human genome (GRCh38).
* Fasta sequence of the human genome (GRCh38), accompanied by the ``.dict`` and``.fai`` files.
* ``1000G_omni2.5.hg38.vcf.gz``
* dbSNP in VCF.
* HLA reference fasta file.


It is recommended that these reference files are downloaded from the Broad Institute's [resource bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) to avoid incompatibilities. The latter file can be obtained from the [Optitype repository](https://github.com/FRED-2/OptiType/tree/master/data). 

### Obtain a copy of this workflow

Clone the repository:
```sh
git clone https://github.com/jucosie/snakemake-pgx.git
```
Or download the ZIP.

### Configure the workflow
Configure the workflow according to your needs via editing the files in the ``config/`` folder. Adjust ``config.yaml`` with the absolute paths to the reference materials, and ``samples.tsv`` to specify your sample setup. 

* ``positions_coverage.bed`` contains the positions of interest from which the sequencing depth information will be extracted.
* ``variants_sorted_0.bed`` is an auxiliary file for the script that generates the report.
* ``VC_positions_chr_sorted_0.bed`` contains the intervals over which Variant Calling shall be applied.

The resources (CPUs and memory) required for each rule can be adjusted according to the available capacities. 

### Install Snakemake
Install Snakemake using Conda:
```sh
conda create -c bioconda -c conda-forge -n snakemake snakemake
```
For installation details, see the instructions in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

<!-- USAGE EXAMPLES -->
## Usage
This workflow was designed modularly, so that rules can be executed one at a time or all at the same time (by executing the final rule, ``all_report``). The names of the general rules can be found in the ``Snakefile``.

First, we activate the conda environment:
```sh
conda activate snakemake
```
If you want to run the workflow locally:
```sh
snakemake --use-conda -c <num_cores> all_<rule>
```
Snakemake also offers the possibility to run in different computing environments (Slurm, SGE, etc.). To do so, it is necessary to generate specific configurations following the instructions available [here](https://github.com/Snakemake-Profiles/doc).

NOTE: Be careful, if you execute the Variant Calling rule more than once, you must delete the directory that creates the GenomicsDB command before each execution.

## Example
We provide as an example of the expected results of the pipeline the final reports of two publicly available whole exome samples:
* [NA11829](https://www.internationalgenome.org/data-portal/sample/NA11829) (SRR710128).
* [NA07000](https://www.internationalgenome.org/data-portal/sample/NA07000) (SRR766039).


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Julia Corell - jucosie@alumni.uv.es

Project Link: [https://github.com/jucosie/snakemake-pgx](https://github.com/jucosie/snakemake-pgx)







<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username
