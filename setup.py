import setuptools

with open('README.md', 'r') as handle:
    long_description = handle.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name='synbio-tooling',
    version='0.1.1',
    author='Jonathan Calles',
    author_email='callesjonathan@gmail.com',
    description='Tools for doing synthetic biology',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jecalles/synbio',
    packages=setuptools.find_packages(),
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
    install_requires=requirements,
    include_package_data=True
)
