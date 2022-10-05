from __future__ import print_function
import base64
import zlib

# ------------------------------------------------------------------------------
# FileList 
# 
# compresses file lists to work around argument size limits when passing lists
# of many files (with absolute path)
# ------------------------------------------------------------------------------
class FileList(object):
    def __init__(self):
        pass

    @staticmethod
    def compress(fileList):
        if type(fileList) == list:
            fileListString = ';'.join(fileList)
        else:
            fileListString = fileList
        compressedList = 'base64:' + base64.b64encode(zlib.compress(fileListString, 9))
        return compressedList
    
    @staticmethod
    def decompress(fileList):
        if fileList.startswith('base64:'):
            decompressedList = zlib.decompress(base64.b64decode(fileList[7:])).split(';')
        else:
            print ("\x1b[31mERROR: invalid compressed file list format: {fileList}\x1b[0m".format(fileList=fileList))
            raise Exception("InvalidCompressedFileListFormat")
        return decompressedList

