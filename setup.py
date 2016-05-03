from distutils.core import setup
setup(
    name='IlluminaBeadArrayFiles',
    description='Library to read data format related to Illumina genotyping bead arrays',
    author = 'Illumina',
    author_email = 'ryankelley@illumina.com',
    packages = ['IlluminaBeadArrayFiles'],
    package_dir = {'IlluminaBeadArrayFiles' : 'module'},
    version = '0.9.0'
)
