import ez_setup

ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open("crazydoc/version.py").read())  # loads __version__

setup(
    name="crazydoc",
    version=__version__,
    author="Zulko",
    url="https://github.com/Edinburgh-Genome-Foundry/crazydoc",
    description="Read genetic sequences from stylized docx files",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="dna-sequence bioinformatics systems-biology docx",
    packages=find_packages(exclude="docs"),
    install_requires=["Biopython", "python-docx"],
)
