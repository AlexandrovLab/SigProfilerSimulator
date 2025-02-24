from setuptools import setup, find_packages
import os
import shutil


# remove the dist folder first if exists
if os.path.exists("dist"):
    shutil.rmtree("dist")


def readme():
    this_directory = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(this_directory, "README.md"), encoding="latin-1") as f:
        long_description = f.read()
        return long_description


VERSION = "1.2.0"


def write_version_py(filename="SigProfilerSimulator/version.py"):
    # Copied from numpy setup.py
    cnt = """
# THIS FILE IS GENERATED FROM SIGPROFILERSIMULATOR SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
Update = 'v1.2.0: Require Pandas and Numpy >= 2.0.0 and Python >= 3.9'
	
	"""
    fh = open(filename, "w")
    fh.write(
        cnt
        % {
            "version": VERSION,
        }
    )
    fh.close()


write_version_py()

setup(
    name="SigProfilerSimulator",
    version=VERSION,
    description="SigProfiler simulator tool",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="",
    author="Erik Bergstrom",
    author_email="ebergstr@eng.ucsd.edu",
    license="UCSD",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "SigProfilerMatrixGenerator>=1.3.0",
        "sigProfilerPlotting>=1.4.0",
        "fastrand>=1.2",
    ],
    include_package_data=True,
    zip_safe=False,
)
