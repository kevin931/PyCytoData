import setuptools
import os
import sys
import shutil
import distutils.cmd

from typing import List

VERSION = "0.1.3"

class PypiCommand(distutils.cmd.Command):
    
    description = "Build and upload for PyPi."
    user_options = []
    
    def initialize_options(self):
        pass
    
    
    def finalize_options(self):
        pass
    
    
    def run(self):
        shutil.rmtree("dist/")
        
        wheel_file = "PyCytoData-{}-py3-none-any.whl".format(VERSION)
        tar_file = "PyCytoData-{}.tar.gz".format(VERSION)
        
        os.system("{} setup.py sdist bdist_wheel".format(sys.executable))
        os.system("twine upload dist/{} dist/{}".format(wheel_file, tar_file))
    
    
class CondaCommand(distutils.cmd.Command):
    
    description = "Build and upload for conda."
    user_options = []
    
    
    @staticmethod
    def move_assets(origin: str, destination: str, exclude: List[str], new_destination_dir: bool) -> None:
        
        if origin[-1] != "\\" and origin[-1] != "/":
            origin += "/"
            
        if destination[-1] != "\\" and destination[-1] != "/":
            destination += "/"
        
        if new_destination_dir:    
            if os.path.isdir(destination):
                raise ValueError("Destination directory already exists.")
            else:
                os.mkdir(destination)
            
        all_files = os.listdir(origin)
        
        for files in all_files:
            if files in exclude:
                pass
            else:
                origin_path = origin + files
                destination_path = destination + files
                shutil.move(origin_path, destination_path)
        
    
    def initialize_options(self):
        pass
    
    
    def finalize_options(self):
        pass
    
    
    def run(self):
        self.move_assets("./PyCytoData/data/", "../temp_assets/", [], True)
        shutil.rmtree("./PyCytoData/data/")
        try:
            shutil.rmtree("dist_conda/")
        except FileNotFoundError:
            pass
        os.system("conda build . --output-folder dist_conda/ -c bioconda")
        os.system("anaconda upload ./dist_conda/noarch/pycytodata-{}-py_0.tar.bz2".format(VERSION))
        
        self.move_assets("../temp_assets/", "./PyCytoData/data/", [], True)
        shutil.rmtree("../temp_assets/")


setuptools.setup(
    name = "PyCytoData",
    version = VERSION,
    description = "An Elegant Data Analysis Tool for CyTOF",
    author="PyCytoData Developers",
    url="https://github.com/kevin931/PyCytoData",
    long_description_content_type = "text/markdown",
    long_description = open("README.md").read(),
    license="MIT",
    packages=["PyCytoData"],
    python_requires=">=3.7",
    install_requires=["fcsparser", "pandas", "numpy>=1.20"],
    test_requires=["pytest",
                   "pytest-cov",
                   "pytest-mock",
                   "coverage"],
    classifiers = [
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Natural Language :: English"
    ],
    cmdclass = {"pypi": PypiCommand,
                "conda": CondaCommand
                }
)