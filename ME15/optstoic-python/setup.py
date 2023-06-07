"""Credit: https://github.com/pypa/sampleproject"""
from setuptools import setup, find_packages
from codecs import open
import os

current_dir = os.path.dirname(os.path.abspath(__file__))

with open(os.path.join(current_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

install_requires = [
    'pandas>=0.18.0',
    'xlrd',
    'scipy>=0.17.0',
    'numpy>=1.11.1',
    'sympy',
    'graphviz>=0.4.8',
    'PuLP>=1.6.1',
    'future'
]

test_requires = [
    'nose'
]

setup(
    name='optstoicpy',
    version='0.5.0',
    description='optStoic python package',
    long_description=long_description,
    url='http://www.maranasgroup.com/software.htm',
    author='Chiam Yu Ng',
    author_email='ngchiamyu@gmail.com',
    license='GNU GPLv3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.8',
    ],
    packages=find_packages(exclude=['build',
                                     'data',
                                     'docs',
                                     'examples']),
    install_requires=install_requires+test_requires,
    test_suite='nose.collector',
    tests_require=test_requires,
    extras_require={
        "Jupyter": ['notebook', 'ipykernel']
    },
    package_dir={'optstoicpy': 'optstoicpy'},
    package_data={
        'optstoicpy': [ 'data/*.csv',
                        'data/*.json',
                        'data/optstoic_db_v3/*.txt',
                        'data/optstoic_db_v3/*.json',
                        'data/optstoic_db_v3/*.pkl'],
    },
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('data', ['data/cofactors.csv',
    #                       'kegg_compound.json'])],

)