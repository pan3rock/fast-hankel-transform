#!/usr/bin/env python
from setuptools import setup
import os
package_name = 'fht'
files_so = [t for t in os.listdir(package_name) if t.endswith('.so')]
print(files_so)

setup(name='fht',
      version='0.1',
      description='fast hankel transform',
      author='Lei Pan',
      author_email='panlei7@gmail.com',
      license='MIT',
      packages=['fht'],
      package_dir={'': '.'},
      package_data={'': files_so},
      install_requires=[
          'scipy',
          'numpy',
          'matplotlib',
      ],
      zip_safe=False)
