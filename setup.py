from setuptools import setup,find_packages
from shigeifinder import __version__

def readme():
    with open('README.md') as f:
        return f.read()


setup(name='shigeifinder',
      version=__version__,
      description='In silico clustering and serotyping of Shigella and Enteroinvasive E. coli',
      long_description=readme(),
      classifiers=[
          'License :: OSI Approved :: GPLv3',
          'Programming Language :: Python :: 3.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.',
          'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics Shigella EIEC serotyping',
      url='https://github.com/LanLab/ShigEiFinder',
      author='Michael Payne',
      author_email='michael.payne@unsw.edu.au',
      license='GPLv3',
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['shigeifinder=shigeifinder.shigeifinder:main'],
      },
      zip_safe=False)
