# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open("README.md", "r") as f:
    readme = f.read()

classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11"
        ]

entry_points = {
    'console_scripts': [
        'resoncalc = resoncalc.input:process_command'
    ]
}

setup(
    name='resoncalc',
    version='0.1.0',
    description='Calculate resonances eigenstates for radial symmetric potential',
    long_description=readme,
    author='Petr Rasek',
    author_email='petr.rasek.w@gmail.com',
    url='https://github.com/hydratk/resoncalc',
    license='MIT',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    classifiers=classifiers,
    zip_safe=False,
    entry_points=entry_points,
    keywords='resoncalc',
    install_requires=[
    'numpy>=1.25.1',
    'scipy>=1.11.1',
    'sympy>=1.12',
    'matplotlib>=3.7.2'
],
)
