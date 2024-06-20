## RFDiffusion Jobs
### Setup
- Clone the RFDiffusion Repository and the path of the cloned RFDiffusion directory in the `.env` file.
- Setup a conda environment `SE3nv` following the installation instructions as mentioned in the RFDiffusion repo.

To start RF Diffusion Jobs, please run python run_rfdiffusion.py with the required flags. 

<!-- TODO: ADD INFO ABOUT THE REQUIRED FLAGS -->

The job command for six different job types are mentioned in `rfdiffusion_jobs.sh`.

The job name can be broken down as follows:
- tr: trimmed receptor/binder
- r366_381: Diffused over residues between 366 and 381
- fixed: represents that a fixed length of residues are diffused for all variants.
- flexible: variable length of residues are diffused in different variants.
