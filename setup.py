from setuptools import setup, find_packages

setup(
    name='PAModelpy',
    version='0.0.5.1',
    packages=find_packages(),
    install_requires=[
        'PAModelpy',
        'matplotlib==0.1.6',
        'scipy',
        'resource',
        'jupyter',
        'gurobipy==9.5.2'
    ]
)
