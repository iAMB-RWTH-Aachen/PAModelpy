from setuptools import setup, find_packages

setup(
    name='PAModelpy',
    version='0.0.4.8',
    packages=find_packages(),
    install_requires=[
        'PAModelpy',
        'matplotlib',
        'scipy',
        'resource',
        'jupyter'
    ]
)
