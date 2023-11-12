import sys
if sys.version_info < (3,):
    sys.exit('NvTK requires Python >= 3.7')
from pathlib import Path
from setuptools import setup, find_packages

try:
    from NvTK import __author__, __email__, __version__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = 'Jiaqili@zju.edu.cn'

setup(
    name='NvTK',
    version=__version__,
    # cmdclass=versioneer.get_cmdclass(),
    description='NvTK is a deep learning package for mapping sequence to single-cell data.',
    long_description=Path('README.rst').read_text('utf-8'),
    url='https://github.com/JiaqiLiZju/NvTK',
    author=__author__,
    author_email=__email__,
    maintainer='hanyuwu@zju.edu.cn',
    maintainer_email='hanyuwu@zju.edu.cn',
    license='BSD',
    python_requires='>=3.7',
    install_requires=[
        'numpy',
        'h5py=2.10.0',
        'scikit-learn',
        'tqdm',
        'torch',
        'networkx',
        'captum=0.5.0',
        'tensorboard',
        'pillow',
    ],
    extras_require=dict(
        tune=['ray[tune]=1.10.0'],
        doc=['sphinx', 'sphinx_rtd_theme'],
    ),
    packages=find_packages(),
    include_package_data=True,
    # entry_points=dict(
    #     console_scripts=['nvtk=nvtk.cli:console_main'],
    # ),
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Deeplearning',
    ],
)
