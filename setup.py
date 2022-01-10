import setuptools

VERSION = "0.0.1"

setuptools.setup(
    name = "PythonTemplate",
    version = VERSION,
    description = "A Template for your Python Package",
    packages=["PythonTemplate"],
    python_requires=">=3.9",
    install_requires=[],
    test_requires=["pytest",
                   "pytest-cov",
                   "pytest-mock",
                   "coverage"],
    classifiers = [
        "Programming Language :: Python :: 3 :: Only",
        "Natural Language :: English"
    ]
)