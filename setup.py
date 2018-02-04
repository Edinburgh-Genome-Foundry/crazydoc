import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, find_packages

exec(open('crazydoc/version.py').read()) # loads __version__

setup(
    name='crazydoc',
    version=__version__,
    author='Zulko',
    description='Read genetic sequences from stylized docx files',
    long_description=open('README.rst').read(),
    license='see LICENSE.txt',
    keywords="dna-sequence bioinformatics systems-biology docx",
    packages=find_packages(exclude='docs'),
    install_requires=['Biopython', 'python-docx'])
