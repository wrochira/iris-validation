import setuptools

with open('README.md', 'r') as infile:
    readme = infile.read()

setuptools.setup(name='iris-validation',
                 version='0.0.0',
                 author='William Rochira',
                 author_email='william.rochira@hotmail.co.uk',
                 description='A package for interactive all-in-one graphical validation of 3D protein model iterations',
                 long_description=readme,
                 long_description_content_type='text/markdown',
                 package_data={ 'iris_validation': [ 'interface/template/*',
                                                     'interface/template/*/*',
                                                     'metrics/*/data/*' ]},
                 include_package_data=True,
                 url='https://github.com/wrochira/iris-validation',
                 packages=setuptools.find_packages(),
                 classifiers=[ 'Programming Language :: Python :: 2.7',
                               'License :: OSI Approved :: MIT License',
                               'Operating System :: OS Independent' ],
                 python_requires='>=2.7, <3')
