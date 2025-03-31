from setuptools import setup, find_packages

version = {}
with open("crazydoc/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="crazydoc",
    version=version["__version__"],
    author="Zulko",
    url="https://github.com/Edinburgh-Genome-Foundry/crazydoc",
    description="Read genetic sequences from stylized docx files",
    long_description=open("pypi-readme.rst").read(),
    license="MIT",
    keywords="dna-sequence bioinformatics systems-biology docx",
    packages=find_packages(exclude="docs"),
    install_requires=["Biopython", "python-docx"],
)
