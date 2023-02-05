import os

from setuptools import setup
from setuptools.command.build_py import build_py

with open("README.md", "rt") as fh:
    long_description = fh.read()


class PublishCommand(build_py):
    """Publish package to PyPI"""
    def run(self):
        """Publish"""
        os.system("rm -rf dist")
        os.system("python3 setup.py sdist"
                  "&& python3 setup.py bdist_wheel"
                  "&& twine upload dist/*whl dist/*gz")


with open("requirements.txt", "rt") as f:
    requirements = [r.strip() for r in f.readlines()]


setup(
    name='sma_finder',
    version="1.2",
    description="A tool for diagnosing spinal muscular atrophy (SMA) using exome or genome sequencing data",
    install_requires=requirements,
    cmdclass={
        'publish': PublishCommand,
    },
    entry_points = {
        'console_scripts': [
            'sma_finder = sma_finder:main',
        ],
    },
    packages=[],
    long_description_content_type="text/markdown",
    long_description=long_description,
    python_requires=">=3.7",
    license="MIT",
    scripts=["sma_finder.py"],
    keywords=["Spinal Muscular Atrophy", "SMA", "SMN", "SMN1", "SMN2"],
    test_suite="setup.test_suite",
    url='https://github.com/broadinstitute/sma_finder',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
)
