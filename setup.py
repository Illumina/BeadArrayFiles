from distutils.core import setup
setup(
    name='IlluminaBeadArrayFiles',
    description='Library to read data format related to Illumina genotyping bead arrays',
    author='Illumina',
    author_email='ryankelley@illumina.com',
    packages=['IlluminaBeadArrayFiles'],
    package_dir={'IlluminaBeadArrayFiles' : 'module'},
    install_requires=[
        'future',
        'numpy',
    ],
    version='1.3.3'
)
