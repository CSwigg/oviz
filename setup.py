from pathlib import Path

from setuptools import find_namespace_packages, setup

setup(
    name="oviz",
    version="0.1.0",
    description="Interactive 3D, Sky, and time visualization for Gaia and Galactic ISM data.",
    long_description=Path("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Cameren Swiggum",
    author_email="cameren.swiggum@univie.ac.at",
    url="https://github.com/CSwigg/oviz/",
    license_expression="MIT",
    # Namespace discovery also lists the data-only directories (themes and the
    # viewer's web sources), so their files ship as package data.
    packages=find_namespace_packages(include=["oviz", "oviz.*"]),
    include_package_data=True,  # Includes files specified in MANIFEST.in
    package_data={
        "oviz.themes": ["*.yaml"],
        # The viewer runtime is bundled from these sources at write time.
        "oviz.viewer": ["web/template.html", "web/styles/*.css", "web/src/**/*.js", "web/assets/**/*.jpg"],
    },
    install_requires=[
        "numpy",
        "pandas",
        "astropy",
        "galpy",
        "scipy",
        "matplotlib",  # volume colormaps
        "pillow",  # PNG volume atlases and Paper images
        "pyyaml",
        "webcolors",
    ],
    extras_require={
        "docs": ["sphinx>=7", "sphinx-rtd-theme>=2"],
        "test": ["pytest"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    # importlib.resources reads the theme directory as a namespace package,
    # which needs Python 3.10.
    python_requires=">=3.10",
)
