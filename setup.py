from setuptools import find_packages, setup

setup(
    name="egfr_microsim",
    version="1.0.0",
    description="eGFR based CKD progression microsimulation model",
    url="https://github.com/StanfordHPDS/data_transformation",
    packages=find_packages(),
    install_requires=[
        "GitPython==3.1.43",
        "matplotlib==3.9.0",
        "pandas==2.2.2",
        "scipy==1.14.0",
        "seaborn==0.13.2",
        "setuptools==70.0.0",
        "jinja2==3.1.4",
        "pyarrow==14.0.2",
        "progressbar2==4.2.0",
        "tqdm==4.66.4",
        "scikit-learn",
        "numpy==1.26.4",
    ],
)
