from setuptools import setup, find_packages

setup(
    name='PAModelpy',
    version='0.0.4.7',
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
