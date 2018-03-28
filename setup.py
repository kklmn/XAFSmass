from distutils.core import setup

long_description = r"""
A program for calculating the mass of XAFS [X-Ray Absorption Fine Structure]
samples. The chemical formula parser understands parentheses and weight
percentage, also in nested form. XAFSmass reports the quantity (weight,
thickness or pressure) together with the expected height of the absorption
edge. The GUI is provided by Qt. The documentation is included.

Dependencies: numpy, pyparsing and matplotlib are required. Qt must be provided
by either PyQt4, PyQt5 or PySide.
"""

setup(name='XAFSmass',
      version='1.3.8',
      description='A program for calculating the mass of XAFS samples. '
                  'For synchrotron users.',
      long_description=long_description,
      author='Konstantin Klementiev, Roman Chernikov',
      author_email='konstantin.klementiev@gmail.com, '
                   'rchernikov@gmail.com',
      download_url='https://pypi.org/project/XAFSmass',
      url='http://xafsmass.readthedocs.io',
      platforms='OS Independent',
      license='MIT License',
      packages=['XAFSmass'],
      package_data={'XAFSmass': ['data/*.*', 'help/*.*', 'help/_images/*.*',
                                 'help/_images/math/*.*', 'help/_sources/*.*',
                                 'help/_static/*.*']},
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 3',
                   'License :: OSI Approved :: MIT License',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Scientific/Engineering :: Visualization'],
      )
