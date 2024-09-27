import re
from setuptools import find_packages, setup
import os

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Попытка прочитать requirements.txt, если он существует
if os.path.exists('requirements.txt'):
    with open('requirements.txt') as f:
        requirements = f.read().splitlines()
else:
    # Если файл не найден, используем минимальный набор зависимостей
    requirements = [
        "tqdm",
        "intervaltree",
        "aindex2",
    ]

version = "0.2.11"

setup(
    name="extracTR",
    version=version,
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={
        "": ["README.md"],
        "extractr": ["requirements.txt"],  # Включаем requirements.txt в пакет
    },
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
    install_requires=requirements,
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
            'extracTR = extractr.extractr:run_it',
        ],
    },
)
