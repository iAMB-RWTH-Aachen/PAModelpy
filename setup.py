from setuptools import setup, find_packages

setup(
    name='PAModelpy',
    version='0.0.4.7',
    packages=find_packages(),
    install_requires=[
        'PAModelpy',
        'matplotlib',
        'scipy',
        'resource',
        'jupyter'
    ]
)
