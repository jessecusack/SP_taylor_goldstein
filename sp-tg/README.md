Installation instructions
-------------------------

1. Install miniconda3 (https://docs.conda.io/en/latest/miniconda.html) to manage python.
1. Download this git repository. 
1. Contact me for the raw data, which should be placed in the raw_data directory within the repository.
1. Create and install the conda environment as follows:

First run the following within the base directory of this repository,

```
conda env create -f environment.yml
```

which creates a new conda environment called `sp-tg`.

Then install the kernel for the environment so that it can be used with jupyter lab or jupyter notebooks.

```
conda activate sp-tg
python -m ipykernel install --user --name sp-tg --display-name sp-tg
```

With the kernel activated and the data downloaded, you can then run the reprocessing scripts contained the `reprocessing` directory, e.g.

```
cd reprocessing
python reprocess_MP_data.py
python reprocess_TOWYO_data.py
``` 
     
which will create a lot of files in the `proc_data` directory. 

After following the above instructions, the analysis code (contained in `analysis`) should then run successfully. Just don't forget to either activate the conda environment first (`conda activate sp-tg`) or set up notebooks/lab that uses that environment. 