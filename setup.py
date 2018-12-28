from setuptools import setup

setup(
    version="0.1",
    name='GeneSetDistance',
    packages=['gsd'],
    include_package_data=True,
    install_requires=[
        'scipy',
        'pandas',
        'biomart',
        'anytree',
        'jsonpickle',
        'numpy',
        'sklearn',
        'rpy2',
        'gensim',
        'pyemd',
        'nltk'
    ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"]
)
