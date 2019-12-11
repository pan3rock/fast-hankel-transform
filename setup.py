#!/usr/bin/env python
from setuptools import setup

setup(name='fht',
      version='0.1',
      description='fast hankel transform',
      author='Lei Pan',
      author_email='panlei7@gmail.com',
      license='MIT',
      packages=['fht'],
      install_requires=[
          'scipy',
          'numpy',
          'matplotlib',
      ],
      zip_safe=False)
