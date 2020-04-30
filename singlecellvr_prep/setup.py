import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scvr_prep", # Replace with your own username
    version="1.0",
    author="Huidong Chen",
    author_email="huidong.chen AT mgh DOT harvard DOT edu",
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
          'networkx==2.1'
      ],
    entry_points = {'console_scripts': ['scvr_prep=scvr_prep.command_line:main']}
)