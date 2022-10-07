from __future__ import print_function
import hashlib
import os

class Hash(object):

    def __init__(self, sample, minCut = '1', subCut = None, branches = None, splitFilesChunkSize = -1, debug = False, inputPath=None):
        
        debug = debug or ('XBBDEBUG_HASH' in os.environ)
        # basic cut(+part) hash
        self.hashKey = '%s_%s' % (sample, minCut)
        if splitFilesChunkSize > 0:
            self.hashKey += '_split%d' % (splitFilesChunkSize)
        if inputPath:
            self.hashKey += '_from:%s'%inputPath.strip().strip('/')
        self.hash = hashlib.sha224(self.hashKey.encode('utf-8')).hexdigest()

        # including subcut
        if subCut:
            if debug:
                print ('DEBUG: hash function debug:')
                print (' > \x1b[32mKEY:', self.hashKey, '\x1b[0m')
            self.hashKey = '%s_[%s]' % (self.hash, subCut)
            self.hash = hashlib.sha224(self.hashKey.encode('utf-8')).hexdigest()

        # including branchnames
        if branches:
            if debug:
                print ('DEBUG: hash function debug:')
                print (' > \x1b[32mKEY:', self.hashKey, '\x1b[0m')
            branchNames = ','.join(sorted(branches))
            self.hashKey = '%s_<%s>' % (self.hash, branchNames)
            self.hash = hashlib.sha224(self.hashKey.encode('utf-8')).hexdigest()

        if debug:
            print ('DEBUG: hash function debug:')
            print (' > \x1b[32mKEY:', self.hashKey, '\x1b[0m')
            print (' > \x1b[33mHASH:', self.hash, '\x1b[0m')


    def get(self):
        return self.hash
