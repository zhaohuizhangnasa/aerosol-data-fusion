from setuptools import setup, find_packages
import codecs
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Setting up
setup(
    name = 'pyroscopegriddingcpu',         # Package name
    packages = ['pyroscopegriddingcpu'],   
    version = '1.5.1.0',      
    license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
    description = 'Data fusion package for transforming L2 satellite to L3 spatial-temporal gridded data',  
    long_description = long_description,
    long_description_content_type='text/markdown',
    author = 'Sally Zhao, Neil Gutkin',                 
    author_email = 'zhaosally0@gmail.com',     
    url = 'https://github.com/jwei-openscapes/aerosol-data-fusion',   # github repository  
    keywords = ['data fusion', 'satellite', 'L2', 'L3'],   # Keywords
    #packages = find_packages(),
    entry_points ={
        'console_scripts': [
            'pyroscopegriddingcpu = pyroscopegriddingcpu.gtools:main',
        ],
    },
    install_requires=[    
            'pyhdf',        
            'numpy ==1.23.5',
            'joblib ==1.2.0',
            'netCDF4 ==1.6.3',
            'pandas ==2.0.0',
            'pyyaml ==6.0',
            'argparse',
            
        ],
    classifiers=[
    'Development Status :: 3 - Alpha',     
    'Intended Audience :: Developers',     
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',   
    'Programming Language :: Python :: 3',      #supported versions
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8'
    ],
)
