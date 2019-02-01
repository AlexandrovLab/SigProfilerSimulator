from setuptools import setup, find_packages
import os

def readme():
	with open('README.rst') as f:
		return(f.read())


setup(name='SigProfilerSimulator',
		version='0.1.3',
		description='SigProfiler simulator tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=find_packages(),
		install_requires =[
			"SigProfilerMatrixGenerator",
			"sigProfilerPlotting"],
		include_package_data=True,
		zip_safe=False)
