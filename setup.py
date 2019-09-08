import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="vcfproj",
    version="0.0.3",
    author="Hirak Sarkar",
    author_email="hsarkar@umd.edu",
    description="Project vcf level coordinates from chromosomes to transcript",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hiraksarkar/vcfproj",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=[
        'pandas',
        'numpy',
        'ncls'
    ]
)
