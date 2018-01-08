from setuptools import setup, find_packages

# A list of authors and their email addresses.
authors=("Christopher Woods <chryswoods@gmail.com>, "
         "Lester Hedges <lester.hedges@gmail.com, "
         "Antonia Mey <antonia.mey@gmail.com")

setup(name='BioSimSpace',
      version='1.0',
      description='BioSimSpace: Making biomolecular simulation a breeze.',
      author=authors,
      url='https://github.com/michellab/BioSimSpace',
      license='GPLv2',
      packages=find_packages()
      )
