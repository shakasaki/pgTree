from setuptools import setup, find_packages
version = '0.0.1'

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pgTree',
    version=version,
    packages=find_packages(exclude=('tests', 'docs')),
    include_package_data=True,
    install_requires=[

    ],
    url='',
    license='LGPL v3',
    author='',
    author_email='',
    description='',
    keywords=['']
)
