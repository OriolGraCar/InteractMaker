#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

requirements = ['biopython',
                'pymol',
                'numpy',
                'TK']

readme = open('README.md').read()
setup(
    name='BioMacromplex',
    version='1.0.0',
    author='G.C. Oriol & S.M. Alvaro',
    author_email='oriol.gracar@gmail.com & alvaro.serrano.morras@gmail.com',
    url='https://github.com/OriolGraCar/InteractMaker',
    description='',
    long_description=readme,
    requires=requirements
)