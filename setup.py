from setuptools import setup, find_packages


def load_requirements(filename='requirements.txt'):
    with open(filename, 'r') as file:
        return file.read().splitlines()

setup(
    name='mecc',
    version='0.2.0',
    author='Luke Li',
    author_email='zhongwei.li@mavs.uta.edu',
    description='A elliptic curve package using montgomery arithmetic.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/r5by/montgomery-ecc',
    packages=find_packages(),
    install_requires=[
        # change this local version should use no-cashed pip install cmd:
        #   pip install -e . --no-cache-dir
        'py-mont @ file:///Users/r5by_/git/py-mont#egg=py-mont-0.2.0'
    ] + load_requirements(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    python_requires='>=3.6',
)
