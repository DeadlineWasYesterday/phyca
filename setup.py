#!/usr/bin/python3

import ast, os, sys

try:
    from setuptools import __version__ as setuptools_version
    from setuptools import Command
    from setuptools import Extension
    from setuptools import setup
except ImportError:
    sys.exit(
        "We need the Python library setuptools to be installed. "
        "Try running: python -m ensurepip"
    )

setuptools_version_tuple = tuple(int(x) for x in setuptools_version.split("."))
if setuptools_version_tuple < (70, 1) and "bdist_wheel" in sys.argv:
    # Check for presence of wheel in setuptools < 70.1
    # Before setuptools 70.1, wheel is needed to make a bdist_wheel.
    # Since 70.1 was released including
    # https://github.com/pypa/setuptools/pull/4369,
    # it is not needed.
    try:
        import wheel  # noqa: F401
    except ImportError:
        sys.exit(
            "We need both setuptools AND wheel packages installed "
            "for bdist_wheel to work. Try running: pip install wheel"
        )


# Make sure we have the right Python version.
MIN_PY_VER = (3, 8)
if sys.version_info[:2] < MIN_PY_VER:
    sys.stderr.write(
        ("ERROR: phyca requires Python %i.%i or later. " % MIN_PY_VER)
        + ("Python %d.%d detected.\n" % sys.version_info[:2])
    )
    sys.exit(1)
# elif sys.version_info[:2] == (3, 9):
#     sys.stderr.write("WARNING: phyca support for Python 3.9 is now deprecated.\n")

REQUIRES = ["numpy", "pandas>=2", "matplotlib", "seaborn", "biopython", "BioNick", 'scipy']

PACKAGES = ["phyca"]

EXTENSIONS = []

setup(
    name="phyca",
    version='0.0.2',
    author="Md Nafis Ul Alam",
    author_email="deadlinewasyesterday@gmail.com",
    url="https://github.com/DeadlineWasYesterday/phyca",
    description="Universal ortholog based Phylogenomic toolkit.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: Freely Distributable",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    entry_points={
            'console_scripts': [
                'phyca=phyca:main',
                ],
            },
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    #package_data={'phyca' :['CUS.tsv']},
    include_package_data=True,  # done via MANIFEST.in under setuptools
    install_requires=REQUIRES,
    python_requires=">=%i.%i" % MIN_PY_VER,
)