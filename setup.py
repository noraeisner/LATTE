import pathlib
from setuptools import setup, find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "readme.md").read_text()

# This call to setup() does all the work
setup(
    name="LATTE",
    version="0.1.0",
    description="Lightcurve Analysis Tool for Transiting Exoplanets",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/noraeisner/LATTE",
    author="Nora L. Eisner",
    author_email="nora.eisner@new.ox.ac.uk",
    license="GNU Lesser General Public LIcense v3.0",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(exclude = ("tests",)),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "LATTE=__main__:main",
        ]
    },
)