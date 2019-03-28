import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biotrees",
    version="1.0.0",
    author="Tomás Martínez Coronado, Gabriel Riera",
    author_email="t.martinez@uib.eu, gabriel.riera@uib.es",
    description="A set of tools to work with phylogenetic trees",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/biocom-uib/biotrees",
    packages=["biotrees"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
    ],
)

