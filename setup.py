from setuptools import setup, find_packages

def calculate_version(inputfile):
    version_list =  [x.split('\'')[1] for x in open(inputfile, 'r')
                     if x.startswith('__version__')]
    if version_list:
        return version_list[0]
    else:
        return '1.0'

package_version = calculate_version('./RefProtDB/_version.py')

setup(
    name='RefProtDB',
    version=package_version,
    author='Alex Saltzman',
    packages=find_packages(),
    install_requires=[
        'Click',
    ],
    entry_points="""
    [console_scripts]
    RefProtDB=RefProtDB.cli:cli
    """,
    package_data={
    }
)
