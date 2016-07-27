from setuptools import setup, find_packages, Extension
from distutils.errors import DistutilsError, CCompilerError
import warnings

extensions = [Extension(
        "_fastdtw",
        ["_fastdtw.pyx"],
        language="c++",
        include_dirs=[],
        libraries=["stdc++"]
    )]

classifiers = [
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 3',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Topic :: Scientific/Engineering'
]

kwargs = {
    'name': 'fastdtw',
    'version': '0.2.2',
    'author': 'Kazuaki Tanida',
    'url': 'https://github.com/slaypni/fastdtw',
    'description': 'Dynamic Time Warping (DTW) algorithm with an O(N) time and memory complexity.',
    'license': 'MIT',
    'keywords': ['dtw'],
    'install_requires': ['six'],
    'packages': find_packages(),
    'ext_modules':  extensions,
    'py_modules': ['fastdtw'],
    'classifiers': classifiers
}

try:
    setup(**kwargs)
except SystemExit:
    del kwargs['ext_modules']
    warnings.warn('compilation failed. Installing pure python package')
    setup(**kwargs)
