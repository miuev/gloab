import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gloab",
    version="0.0.1",
    author="Evan Miu",
    author_email="evm@pitt.edu",
    description="generate and analyze graph laplacians from atoms and their bonds",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/miuev/gloab",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    install_requires=['numpy>=1.17.2',
                      'ase>=3.22.1']
)
