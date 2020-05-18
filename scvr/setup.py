import setuptools
import scvr

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scvr", # Replace with your own username
    version=scvr.__version__,
    author="Huidong Chen",
    author_email="huidong.chen@mgh.harvard.edu",
    description="single cell VR preprocess",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pinellolab/singlecellvr",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
          'pandas>=0.21',
          'numpy>=1.14.0',
          'networkx>=2.1',
          'anndata>=0.7',
          'loompy>=2.0',
          'matplotlib>=3.0',
          'scipy>=1.3'
      ],
    entry_points = {'console_scripts': ['scvr=scvr.command_line:main']}
)