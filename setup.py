from setuptools import setup, find_packages

setup(
    name="oviz",
    version="0.1.0",
    description="Interactive 3D, Sky, and time visualization for Gaia and Galactic ISM data.",
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
    ],
    extras_require={
        "docs": ["sphinx>=7", "sphinx-rtd-theme>=2"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
