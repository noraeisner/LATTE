
'''
Unit test for LATTE configuration logic.
'''
import os
import shutil
import unittest

from LATTE.LATTEconfig import LATTEconfig


class TestLATTEconfig(unittest.TestCase):

	config_dir_for_test = os.path.join(os.path.dirname(__file__), '_LATTE_test_config')

	def setUp(self):
		# Note: we do not remove config dir at tearDown(),
		# so that in case tests fail, one can inspect the content
		if os.path.exists(self.config_dir_for_test):
			shutil.rmtree(self.config_dir_for_test)

	def test_typical(self):
		# Test initial output path assignment, starting with no config dir exist
		config = LATTEconfig(self.config_dir_for_test)
		assert config.inited is False

		output_path_a = 'latteDataA'
		config.output_path = output_path_a
		assert config.output_path == output_path_a
		assert config.inited is True

		# Ensure the setting is persisted by using another instance
		config2 = LATTEconfig(self.config_dir_for_test)
		assert config2.output_path == output_path_a

        #
		# Now test updating the path
		#
		output_path_b = 'latteDataB'
		config2.output_path = output_path_b
		assert config2.output_path == output_path_b

		# ensure the setting is persisted
		assert config.output_path == output_path_b


	def test_config_dir_cannot_be_created(self):
		config_dir_for_fail = os.path.join(os.path.dirname(__file__), '_LATTE_test_config_fail')
		with open(config_dir_for_fail, 'w') as f:
			f.write("Config dir failure case: the path is an existing file")

		with self.assertRaises(ValueError) as context:
			config = LATTEconfig(config_dir_for_fail)
		# clean up
		os.remove(config_dir_for_fail)



if __name__ == '__main__':

	unittest.main()

