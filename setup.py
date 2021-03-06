from setuptools import setup, find_namespace_packages
import glob

with open("README.md") as f:
    long_description = f.read()

setup(
    name="amplicon_benchmark",
    version="0.0.2",
    author="John Sundh",
    url="https://github.com/johnne/amplicon_benchmark/",
    description="Python package to test taxonomic assignment tools on "
                "amplicon databases",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",

    python_requires=">=3.8",
    install_requires=[
        "biopython",
        "tqdm",
        "pandas",
        "importlib_resources",
        "snakemake",
    ],
    package_dir={"": "src"},
    packages=find_namespace_packages("src"),
    scripts=glob.glob('src/amplicon_benchmark/scripts/*.py'),
    package_data={
        "amplicon_benchmark": [
            "Snakefile",
            "config.schema.yaml",
            
        ]
    },
    entry_points={
        "console_scripts": [
            "amplicon_benchmark = amplicon_benchmark.__main__:main"
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ]
)