from setuptools import setup
from setuptools import find_packages

# tensorflow-gpu==1.13.1
install_requires = [
    'keras==2.2.4',
    'docopt==0.6.2',
    'schema==0.7.0',
    'pandas==0.23.4',
    'numpy==1.15.4',
    'scipy==1.1.0',
    'matplotlib==3.0.2',
    'sklearn',
    'cooler==0.8.3',
    'hic2cool==0.5.1',
    'm2r==0.2.1'
]

setup(name='DeepHicIntegrator',
      version='0.1',
      description='Integration of chromatin conformation and epigenetic features using Convolutional Autoencoder',
      author='Hélène Kabbech',
      author_email='helene.kabbech@gmail.com',
      url='https://github.com/kabhel/DeepHicIntegrator',
      license='GNU',
      install_requires=install_requires,
packages=find_packages())