from setuptools import setup
import versioneer

requirements = [
    "Cython>=0.29.24",
    "pandas>=1.3.3",
    "scipy>=1.9.0",
    "tqdm>=4.64.1",
    "pysam>=0.17.0",
    "numpy>1.24",
    "pyrle>=0.0.31",
    "sorted-nearest>=0.0.38",
    "pyranges>=0.0.118",
    "kneed>=0.8.1",
    "matplotlib>=3.6.0",
    "taxopy>=0.12.0",
    "python-datatable>=1.1.3",
    "networkx>=3.2.1",
    "psutil>=5.9.8",
    "numba>=0.60.0",
]

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        "setuptools>=18.0",
        "Cython>=0.29.24",
    ],
    name="bam-filter",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple tool to filter references from a BAM file using different filter types",
    license="GNUv3",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/bam-filter",
    packages=["bam_filter"],
    entry_points={"console_scripts": ["filterBAM=bam_filter.__main__:main"]},
    install_requires=requirements,
    keywords="bam-filter",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
