#!/usr/bin/python

import os
import sys
import getopt

host = 'FO_Hebert@10.37.129.2'

mount_points = ['/Users/FO_Hebert/Documents', '/home/hebertfo/Documents_OSX']

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        opts, args = getopt.getopt(argv[1:], 'u', ['unmount'])
    except getopt.error, msg:
        print msg
        print 'for help use --help'
        sys.exit(2)
    mount_op = mount
    for option, a in opts:
        if option in ('-u', '--unmount'):
            mount_op = unmount
    mount_op()    
    return 0

def mount():
    argv = ['/usr/bin/sshfs', None, None]
    for mount in mount_points:
        argv[1] = '%s:%s' % (host, mount[0])
        argv[2] = mount[1]
        cmd = ' '.join(argv)
        print cmd
        os.system(cmd)

def unmount():
    argv = ['/usr/bin/fusermount', '-u', None]
    for mount in mount_points:
        argv[2] = mount[1]
        cmd = ' '.join(argv)
        print cmd
        os.system(cmd)

if __name__ == '__main__':
    sys.exit(main())
