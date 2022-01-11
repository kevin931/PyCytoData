import setuptools

VERSION = "0.0.0"

setuptools.setup(
    name = "PyCytoData",
    version = VERSION,
    description = "A Python Interface to HDCytoData",
    packages=["PyCytoData"],
    python_requires=">=3.9",
    install_requires=["fcsparser", "pandas", "numpy"],
    test_requires=["pytest",
                   "pytest-cov",
                   "pytest-mock",
                   "coverage"],
    classifiers = [
        "Programming Language :: Python :: 3 :: Only",
        "Natural Language :: English"
    ]
)