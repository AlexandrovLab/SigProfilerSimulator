from setuptools import setup, find_packages
import os
import shutil


#remove the dist folder first if exists
if os.path.exists("dist"):
	shutil.rmtree("dist")

def readme():
	with open('README.rst') as f:
		return(f.read())

VERSION = '1.0.2'

def write_version_py(filename='SigProfilerSimulator/version.py'):
	# Copied from numpy setup.py
	cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILERSIMULATOR SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
	
	"""
	fh = open(filename, 'w')
	fh.write(cnt % {'version': VERSION,})
	fh.close()

write_version_py()

setup(name='SigProfilerSimulator',
		version=VERSION,
		description='SigProfiler simulator tool',
		url='',
		author='Erik Bergstrom',
		author_email='ebergstr@eng.ucsd.edu',
		license='UCSD',
		packages=find_packages(),
		install_requires =[
			"SigProfilerMatrixGenerator>=1.0.14",
			"sigProfilerPlotting>=1.0.3",
			"fastrand>=1.2"],
		include_package_data=True,
		zip_safe=False)
