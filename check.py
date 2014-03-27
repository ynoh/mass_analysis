"""
to backup codes and texfile
reference: http://code.activestate.com/recipes/191017-backup-your-files/
"""

import os, sys

def check(data_dir, str_to_check, filetype):
	#backup_list = [".py"]

	for dir, subdirs, files in os.walk(data_dir):
		command = "grep %s %s" % (str_to_check, dir + os.sep + filetype)
		os.system(command)

if __name__ == "__main__":
	print sys.argv

	if len(sys.argv) != 4:
		sys.exit("Usage: type [tree top directory to check] [string to check] [in what types of file]" )	

	tree_top = os.path.abspath(sys.argv[1])

	print "checking if :\"%s\" is directory" % (tree_top)

	if os.path.isdir(tree_top):
		print "directories under %s will be checked" % (sys.argv[1])
		check(sys.argv[1], sys.argv[2], sys.argv[3])
