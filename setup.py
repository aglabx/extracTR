import re

from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

version = "0.2.0"

setup(
    name="extracTR",
    version=version,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={"": ["README.md"]},
    python_requires=">=3.6",
    include_package_data=True,
    scripts=[],
    license="BSD",
    url="https://github.com/aglabx/extracTR",
    author="Aleksey Komissarov",
    author_email="ad3002@gmail.com",
    description="Extract and analyze satellite DNA from raw sequences.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "argparse",
        "networkx",
        "tqdm",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    entry_points={
        'console_scripts': [
            'extracTR = extracTR.extractr:run_it',
        ],
    },
)
