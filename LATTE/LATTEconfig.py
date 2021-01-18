import os


class LATTEconfig(object):
	"""Provide configuration for LATTE"""

	_config_dir = None

	def __init__(self, config_dir=None):
		if config_dir is None:
			config_dir = os.path.join(os.path.expanduser('~'), '.LATTE')

		if os.path.isdir(config_dir):
			self._config_dir = config_dir
			return
		else:
			# if it doesn't exist, make a new directory
			try:
				os.mkdir(config_dir)
				self._config_dir = config_dir
			except OSError as ose:
				msg = "LATTE config directory does not exist, and cannot be created: {}".format(config_dir)
				raise ValueError(msg) from ose

	@property
	def output_path(self):
		"""LATTE's output path."""
		with open("{}/_config.txt".format(self._config_dir), 'r') as f:
			return str(f.readlines()[-1])

	@output_path.setter
	def output_path(self, output_path):
		with open("{}/_config.txt".format(self._config_dir), 'w') as f:
			f.write(str(output_path))

	@property
	def inited(self):
		"""return ``True`` if configuration exists."""
		return os.path.exists("{}/_config.txt".format(self._config_dir))

