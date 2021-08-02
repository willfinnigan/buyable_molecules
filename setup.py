from setuptools import setup, find_packages
from buyable_molecules import __version__

setup(
  name = 'buyable_molecules',
  packages = find_packages(),
  include_package_data=True,
  version = __version__,
  author = 'William Finnigan',
  author_email = 'wjafinnigan@gmail.com',
)

