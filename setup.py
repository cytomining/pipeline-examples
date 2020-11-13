import glob
import os.path

import setuptools

NOTEBOOKS = glob.glob(os.path.join("notebooks", "*.ipynb"))

setuptools.setup(
    install_requires=[
        "notebook>=6.1.5",
        "numpy>=1.19.4",
        "pandas>=1.1.4",
        "plotnine>=0.7.1",
        "pycytominer@https://github.com/cytomining/pycytominer/archive/master.zip",
        "requests>=2.25.0"
    ],
    name="cytominer-pipeline-examples",
    package_data={
        "cytominer_pipeline_examples": NOTEBOOKS
    },
    package_dir={
        "": "src"
    },
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.8, <4",
    version="1.0.0"
)
