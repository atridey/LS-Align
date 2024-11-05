## LS-Align Project
An implementation of the LS-Align small molecule alignment algorithm (DOI: [10.1093/bioinformatics/bty081](https://doi.org/10.1093/bioinformatics/bty081)) using openbabel to add input/output functionality and build a comprehensive tool.

### Setup
Requires scipy, numpy, and openbabel in a **conda environment**—the pip version of openbabel rarely works without causing various issues.
Query and template molecule files of an [openbabel supported format](https://open-babel.readthedocs.io/en/latest/FileFormats/Overview.html) should be included in the working directory.

### Usage
Run `main.py` and enter the query and template filenames as prompted in the console. The aligned template will be exported to `aligned_query.pdb` in the working directory.

### Evaluation
`rmsd.py` is a testing tool of mine that evaluates how well two molecules have been superimposed by calculating their root mean square deviation, which I have included here. To use, input the filenames when prompted—RMSD will be printed to console.
