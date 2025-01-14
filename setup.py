from setuptools import setup, find_packages

setup(
    name="oviz",
    version="0.1.0",
    description="A Python package for visualizing the 3D orbits of star clusters and associations.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Cameren Swiggum",
    author_email="cameren.swiggum@univie.ac.at",
    url="https://github.com/CSwigg/oviz/",
    packages=find_packages(),
    include_package_data=True,  # Includes files specified in MANIFEST.in
    package_data={
        "oviz.themes" : ["*.yaml"],
    },
    install_requires=[
        "numpy",
        "galpy",
        "astropy",
        "pandas",
        "plotly",
        "dash"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)