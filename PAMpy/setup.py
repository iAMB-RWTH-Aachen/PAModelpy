from setuptools import setup, find_packages

setup(
    name='PAMpy',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'PAModelpy',
        'plotly==5.13.0',
        'matplotlib==0.1.6',
        'scipy',
        'time',
        'resource',
        'PIL',
        'jupyter'
    ]
)
