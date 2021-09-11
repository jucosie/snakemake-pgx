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
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

This tool was developed with the aim of providing a unified solution to extract relevant PGx information from a set of specific genes in WES data. It brings together tools for genotyping HLAs (``Optitype``) and haplotyping of pharmacogenes (``Aldy``), together with position-specific variant calling (``GATK4``) and coverage data (``Mosdepth``). 


### Built With

* Python


<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

We have tried to include all the tools necessary for the pipeline to work as wrappers, or via Conda envirotnments. Even so, it is necessary to install locally: 
* [GATK4 (2.1.0)](https://github.com/broadinstitute/gatk). GATK4 requires Java 8 to run. The path to the gatk wrapper script must be in the user's ``PATH`` variable.
* [Picard](https://github.com/broadinstitute/picard/releases). The path to the ``.jar`` file has to be specified in the configuration file. 

In addition, the user needs to download the following reference files and indicate their path in the configuration file.: 
* BWA index of the human genome (GRCh38).
* Fasta sequence of the human genome (GRCh38), accompanied by the ``.dict`` and``.fai`` files.
* ``1000G_omni2.5.hg38.vcf.gz``
* dbSNP in VCF.
* HLA reference fasta file.
It is recommended that these reference files are downloaded from the Broad Institute's [resource bundle](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) to avoid incompatibilities. The latter file can be obtained from the [Optitype repository](https://github.com/FRED-2/OptiType/tree/master/data). 
### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/jucosie/snakemake-pgx.git
   ```
Or download the ZIP.



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



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
