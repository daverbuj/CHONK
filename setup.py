from setuptools import setup, find_packages

setup(
    name='chonk',
    version='0.1.0',
    description='Detect and genotype germline and somatic structural variants from short-read WGS data',
    url='https://github.com/dantaki/chonk',
    author='Danny Antaki, Dan Averbuj',
    author_email='dantaki@ucsd.edu',
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        'pysam',
        'numpy',
        'scikit-learn',
    ],
    entry_points={
        'console_scripts': [
            'chonk=chonk.cli:main',
        ],
    },
)
