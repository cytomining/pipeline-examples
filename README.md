# Perform an example image-base profiling pipeline

We provide an example image-based profiling pipeline using `pycytominer` and `cytominer-eval`.
We use publicly-available Cell Painting datasets to demonstrate how to use tools in the cytominer ecosystem.

## Step 0: Setup

* Step 0.0 - Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).
  * After installing miniconda, restart your terminal. You should now see the (base) prefix on your command line.
  * Optionally, also install [mamba](https://github.com/mamba-org/mamba) with `conda install mamba -c conda-forge`
* Step 0.1 - Install [aws cli](https://docs.aws.amazon.com/cli/latest/userguide/install-macos.html)
* Step 0.2 - Clone this github repository (forking is optional):

```bash
# Make sure you are navigated to the directory of your choice
git clone git@github.com:cytomining/pipeline-examples.git
```

* Step 0.3 - Create the conda environment, which includes pycytominer and cytominer-eval packages.

```bash
# Make sure you navigated into the repository folder after cloning
# cd pipeline-examples
conda env create --force --file environment.yml

# Or, if you installed mamba
mamba env create --force --file environment.yml
```

* Step 0.4 - Activate the conda environment

```bash
conda activate cytominer-example
```

* Step 0.5 - Alternatively, these two packages can be installed via `pip`

```bash
pip install git+https://github.com/cytomining/pycytominer@c1aa34b641b4e07eb5cbd424166f31355abdbd4d
pip install git+https://github.com/cytomining/cytominer-eval@6f9d350badd0a18b6c1a76171813aaf9a52f8d9f
```

## Step 1: Download one plate of single cell Cell Painting data (2GB .SQLite file)

Data are not included in this repository.
You must run the code in `0.download.sh`, which requires the AWS command line interface.

```bash
# Download one example plate from AWS
./download.sh
```

## Step 2: Perform an image-based profiling pipeline

* Step 2.0 - Run the command `jupyter notebook` in your terminal, in the top level directory.

```bash
# Make sure the cytominer-example environment is activated
jupyter notebook
```

* Step 2.1 - Navigate to `1.profile.ipynb` and follow along!

## Step 3: Evaluate profile quality

* Step 3.0 - Navigate to `2.evaluate.ipynb` and follow along!
