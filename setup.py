from distutils.core import setup
from setuptools import find_packages
import alphascreen

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='alphascreen',
      version=alphascreen.__version__,
      description='Streamlined preparation and analysis of AlphaFold protein-protein interactions.',
      author='Sami Chaaban',
      author_email='chaaban@mrc-lmb.cam.ac.uk',
      url='http://pypi.python.org/pypi/alphascreen/',
      #long_description=long_description,
      #long_description_content_type='text/markdown',
      packages=find_packages(),
      entry_points={
          "console_scripts": [
            "alphascreen = alphascreen.__main__:main",
            ],
      },
      install_requires=["pandas","unipressed", "matplotlib", "numpy", "opencv-python", "openpyxl", "biopython"],
      python_requires='>=3.8'
     )
