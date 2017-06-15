from setuptools import setup
from setuptools.command.install import install
from os.path import expanduser
from shutil import copyfile

setup(
    name='Clusperin',
    version='0.1',
    author='Franca AF (Alexander da Franca Fernandes)',
    author_email='alexander@francafernandes.com.br',
    license='BSD',
    description='Tool to cluster proteins.',
    long_description='Tool to cluster proteins using Galpering method.',
    scripts=['bin/clusperin'],
    packages=[ 'clusperin' ],
    platforms='Linux',
    url='http://bioinfoteam.fiocruz.br/clusperin',
    install_requires=[
            'configparser',
            'datetime ',
            'pprint',
            ],
)


