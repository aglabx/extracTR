import re

from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Читаем requirements.txt
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

version = "0.2.3"

setup(
    name="extractr",
    version=version,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={"": ["README.md"]},
    python_requires=">=3.6",
    include_package_data=True,
    scripts=[],
    license="BSD",
    url="https://github.com/aglabx/extractr",
    author="Aleksey Komissarov",
    author_email="ad3002@gmail.com",
    description="Extract and analyze satellite DNA from raw sequences.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,  # Используем прочитанные зависимости
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
            'extractr = extractr.extractr:run_it',
        ],
    },
)
