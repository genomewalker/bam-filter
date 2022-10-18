from setuptools import setup
import versioneer

requirements = [
    "Cython>=0.29.24",
    "pandas>=1.3.3",
    "scipy>=1.9.0",
    "tqdm==4.50.0",
    "pysam>=0.17.0",
    "numpy>=1.21.2",
    "pyrle>=0.0.31",
    "sorted-nearest<=0.0.33",
    "pyranges>=0.0.112",
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
    description="A simple tool to filter references from a BAM  file using different filter types",
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
