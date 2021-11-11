# probably create a package
# TODO find out how this is done best

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="umidedup",
    version="0.5.0",
    author="Nico Borgsmüller, Fabian Crespo, Michael Schneider",
    author_email="anonymous@example.com",
    description="Package to deduplicate UMIs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/umi-clean",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
    ],
)
