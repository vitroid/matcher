#!/usr/bin/env python

# from distutils.core import setup, Extension
from setuptools import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs
import os
import codecs
import re

#Copied from wheel package
here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(os.path.dirname(__file__), 'matcher', '__init__.py'),
                 encoding='utf8') as version_file:
    metadata = dict(re.findall(r"""__([a-z]+)__ = "([^"]+)""", version_file.read()))



long_desc = "".join(open("README.md").readlines())

setup(ext_modules=[Extension("matcher.csmatcher", ["C/csmatcher.c",
                                                   "C/smatcher.c",
                                                   "C/bst.c",
                                                   "C/common.c",
                                                   "C/pairlist.c"],
                                                   extra_compile_args = ["-std=c99",],
                                                   include_dirs=get_numpy_include_dirs()),
                   Extension("matcher.cmatcher2", ["C/cmatcher2.c",
                                                   "C/matcher2.c",
                                                   "C/bst.c",
                                                   "C/common.c",
                                                   "C/pairlist.c",
                                                   "C/svd.c",
                                                   "C/neighborlist.c"],
                                                   extra_compile_args = ["-std=c99",],
                                                   include_dirs=get_numpy_include_dirs())],
      headers=["C/pairlist.h",
               "C/common.h",
               "C/bst.h",
               "C/smatcher.h",
               "C/svd.h",
               "C/neighborlist.h",
               "C/matcher2.h"],
      name='genice_matcher',
      version=metadata['version'],
      zip_safe=False,
      description='Match atomic environments.',
      long_description=long_desc,
      long_description_content_type="text/markdown",
      classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.5",
        ],
      author='Masakazu Matsumoto',
      author_email='vitroid@gmail.com',
      url='https://github.com/vitroid/matcher/',
      keywords=['matcher',],
      license='MIT',
      packages=['matcher',],
      install_requires=['numpy', 'yaplotlib'],
      entry_points = {
          'genice_format_hook0': [
              'smatcher = matcher.smatcher:hook0',
              'matcher2 = matcher.matcher2:hook0',
              ],
          'genice_format_hook1': [
              'smatcher = matcher.smatcher:hook1',
              'matcher2 = matcher.matcher2:hook1',
          ],
          'console_scripts': [
              'matcher2 = matcher.matcher2:main',
              'matcher-visualize = matcher.matcher2yap:main',
              'smatcher-visualize = matcher.smatcher2yap:main',
          ],

          }
)
