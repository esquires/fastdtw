from setuptools import setup, find_packages, Extension

extensions = [Extension(
        "fastdtw",
        ["fastdtw.pyx"],
        language="c++",
        include_dirs=[],
        libraries=["stdc++"],
        extra_link_args=['-std=c++11'],
        extra_compile_args=['-std=c++11'],
    )]

for e in extensions:
    e.cython_directives = {'boundscheck': False,
                           'cdivision': True,
                           'wraparound': False}
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
    'ext_modules' :  extensions,
    'classifiers': classifiers
}

try:
    setup(**kwargs)
except: 
    del kwargs['ext_modules']
    from pprint import pprint
    pprint(kwargs)
    kwargs['py_modules'] = ['fastdtw']
    print('cython not available on this computer. '
          'Falling back on pure python package')
    setup(**kwargs)
