from setuptools import setup, find_packages

setup(
    name='PAModelpy',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'PAModelpy',
        'matplotlib==0.1.6',
        'scipy',
        'time',
        'resource',
        'PIL',
        'jupyter'
    ]
)
