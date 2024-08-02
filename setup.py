import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Cthulhu", # Replace with your own username
    version="0.9.5",
    author="Ryan MacDonald, Arnav Agrawal",
    author_email="ryanjmac@umich.edu, aa687@cornell.edu",
    description="A python package to calculate atomic and molecular cross sections for substellar atmospheres.",

    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MartianColonist/Cthulhu",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7.6',
    install_requires = ['numpy',
                        'scipy',
                        'matplotlib',
                        'h5py',
                        'numba>=0.56',
                        'requests',
                        'bs4',
                        'tqdm',
                        'pandas',
                        'lxml',
                        'hitran-api'],
    zip_safe=False
)
