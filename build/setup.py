from setuptools import setup
setup(
	name='chonk',
	version='0.0.1',
	url='https://github.com/dantaki/chonk',
	author=['Danny Antaki','Dan Averbuj'],
	author_email='dantaki@ucsd.edu',
	packages=['chonk'],
	package_dir={'chonk': 'chonk/'},
	include_package_data=True,
	scripts= ['chonk/chonk']
)
