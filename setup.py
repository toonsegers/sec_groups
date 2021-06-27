import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="secgroups",
    version="0.0.1",
    description="Implements the Secure Group scheme.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/toonsegers/sec_groups/",
    author="Toon Segers",
    author_email="a.j.m.segers@tue.nl",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
    ],
    packages=["sec_groups", "sec_groups.tools"],
    include_package_data=True,
    # install_requires=["mpyc"],
    python_requires='>=3.6',
)
