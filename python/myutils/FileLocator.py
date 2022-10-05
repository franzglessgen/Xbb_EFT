from __future__ import print_function
import os
import shutil
import sys
import subprocess
import ROOT
import hashlib
import fnmatch
import glob
from time import sleep

# TODO: make it possible to have several Xrootd clients active at the same time
# e.g. self.clients = {'t3dcachedb03.psi.ch': <client object>, 'eoscms.cern.ch': <client object>} ...

class FileLocator(object):

    def __init__(self, config=None, xrootdRedirector=None, usePythonXrootD=True, useDirectoryListingCache=False):
        self.config = config
        self.debug = 'XBBDEBUG' in os.environ
        self.timeBetweenAttempts = 1
        self.useDirectoryListingCache = useDirectoryListingCache
        self.dirListingCache = {}
        self.failureCounters = {'pyxrd': 0, 'xrdfs': 0}

        self.disableXrootdForDirlist = False
        if self.config is not None:
            if self.config.has_option('Configuration', 'disableXrootdForDirlist'):
                self.disableXrootdForDirlist = eval(self.config.get('Configuration', 'disableXrootdForDirlist'))

        try:
            self.xrootdRedirectors = [x.strip() for x in self.config.get('Configuration', 'xrootdRedirectors').split(',') if len(x.strip())>0]
            #print('self.xrootdRedirectors', self.xrootdRedirectors) #['root://t3dcachedb03.psi.ch:1094/']
            if len(self.xrootdRedirectors) < 1:
                print("WARNING: empty list of xrootd redirectors given!")
                
        except Exception as e:
            if config:
                print (e)
                print("WARNING: no xrootd redirector given!")
            self.xrootdRedirectors = None

        # needed for cases there are /store/ path accessible locally but with additional prefix path
        try:
            self.pnfsStoragePath = self.config.get('Configuration', 'pnfsStoragePath').strip()
            #print('self.pnfsStoragePath', self.pnfsStoragePath) #.
        except:
            if config:
                print("WARNING: no pnfs storage path given!")
            self.pnfsStoragePath = None

        # if redirector is given as argument, use it as the main redirector
        if xrootdRedirector:
            #print('xrootdRedirector')
            if self.xrootdRedirectors:
                self.xrootdRedirectors = [xrootdRedirector]+self.xrootdRedirectors
            else:
                self.xrootdRedirectors = [xrootdRedirector]
        #print('self.xrootdRedirectors', self.xrootdRedirectors) #['root://t3dcachedb03.psi.ch:1094/']

        # use python bindings for xrootd (can be disabled setting the 'usePythonXrootD' option to False
        self.usePythonXrootD = eval(self.config.get('Configuration', 'usePythonXrootD')) if self.config and self.config.has_option('Configuration', 'usePythonXrootD') else True
        self.client = None
        self.server = None
        if self.usePythonXrootD or usePythonXrootD:
            try:
                import XRootD.client
                self.server = self.xrootdRedirectors[0].strip('/')
                self.client = XRootD.client.FileSystem(self.server)
                if self.debug:
                    print('DEBUG: initialized xrootd client, server:', self.server)
                    print('DEBUG: client:', self.client)
            except Exception as e:
                if self.debug:
                    print('DEBUG: xrootd:', e)
                    print('DEBUG: xrootd could not be initialized, trying to use xrdfs as fallback. To use the faster Python bindings upgrade CMSSW to version 9.')

        # prefixes to distinguish remote file paths from local ones
        self.storagePathPrefix = '/store/'
        self.pnfsPrefix = '/pnfs/'
        self.remotePrefixes = [self.storagePathPrefix, self.pnfsPrefix]
        if self.config and self.config.has_option('Configuration', 'remotePrefixes'):
            self.remotePrefixes = self.config.get('Configuration', 'remotePrefixes').split(',')
        #print('self.remotePrefixes', self.remotePrefixes)  #self.remotePrefixes ['/store/', '/pnfs/']

        # TODO: use XrootD python bindings
        self.remoteStatDirectory = 'xrdfs {server} stat -q IsDir {path}'
        self.remoteStatFile = 'xrdfs {server} stat {path}'
        self.remoteMkdir = 'xrdfs {server} mkdir {path}'
        self.remoteRm = 'xrdfs {server} rm {path}'
        self.remoteCp = 'xrdcp -d 1 -f --posc --nopbar {source} {target}'
        self.makedirsMinLevel = 5   # don't even try to create/access the 5 lowest levels in the path

    # special Xbb function: get filename after prep step from original file name
    def getFilenameAfterPrep(self, inputFile):
        inputFile = inputFile.replace('\r','')
        # this allows to have duplicate tree names, e.g. to import 2 trees which are both called tree_1.root
        if self.config and self.config.has_option('Configuration', 'AllowDuplicateTrees') and eval(self.config.get('Configuration', 'AllowDuplicateTrees')):
            subfolder = '_'.join(inputFile.split('/')[-4:-1])
        else:
            subfolder = inputFile.split('/')[-4]
        filename = inputFile.split('/')[-1]
        filename = filename.split('_')[0]+'_'+subfolder+'_'+filename.split('_')[1]
        hash = hashlib.sha224(filename).hexdigest()
        return filename.replace('.root','')+'_'+str(hash)+'.root'

    def getFilePath(self, basePath, sampleIdentifier, originalFileName):
        return "{path}/{subfolder}/{filename}".format(path=basePath, subfolder=sampleIdentifier, filename=self.getFilenameAfterPrep(originalFileName))

    # check if path is relative to /store/
    def isStoragePath(self, path):
        return path.startswith(self.storagePathPrefix)

    # check if a path is a PNFS directory (with or without redirector)
    def isPnfs(self, path):
        return self.pnfsPrefix in path

    # check if the path is remote (= has to be accessed via xrootd)
    def isRemotePath(self, path):
        return self.isPnfs(path) or self.isStoragePath(path) or '://' in path

    def isValidRootFile(self, path):
        if self.debug:
            print("DEBUG: check validity of ", path)
        f = ROOT.TFile.Open(path, 'read')
        if f:
            isValid = not (f.IsZombie() or f.GetNkeys() == 0 or f.TestBit(ROOT.TFile.kRecovered))
            try:
                f.Close()
            except:
                pass
            return isValid
        else:
            return False

    def runCommand(self, command, forceDebug=False):
        if self.debug or forceDebug:
            print ("RUN: \x1b[32m",command,"\x1b[0m")
            return subprocess.call([command], shell=True)
        else:
            with open(os.devnull, 'w') as fp:
                result = subprocess.call([command], shell=True, stdout=fp, stderr=fp)
            return result 

    def getRemoteFileserver(self, path=None):
        if path and '://' in path:
            server = path.split('://')[1].split('/')[0].split(':')[0].strip()
        else:
            server = self.xrootdRedirectors[0].split('://')[1].split('/')[0].split(':')[0].strip()
        return server

    def getRedirector(self, path):
        if '://' in path:
            redirector = path.split('://')[0] + '://' + path.split('://')[1].split('/')[0]
        else:
            redirector = None
        return redirector

    def getHostnameFromRedirector(self, redirector):
        return redirector.replace('root://','').split(':')[0].strip()

    def remoteDirectoryExists(self, path):
        statCommand = self.remoteStatDirectory.format(server=self.getRemoteFileserver(path), path=self.getLocalFileName(path))
        result = self.runCommand(statCommand)
        return result==0

    def remoteFileExists(self, path):
        if self.client:
            pathOnServer = self.removeRedirector(path.strip())
            response = self.client.stat(pathOnServer)
            existing = response[0].ok
            if self.debug:
                print('DEBUG: remoteFileExists("' + pathOnServer + '")')
                print('DEBUG: remoteFileExists() returned', response[0])
        else:
            statCommand = self.remoteStatFile.format(server=self.getRemoteFileserver(path), path=self.getLocalFileName(path))
            result = self.runCommand(statCommand)
            existing = result==0
        return existing

    def remoteFileRm(self, path):
        command = self.remoteRm.format(server=self.getRemoteFileserver(path), path=self.getLocalFileName(path))
        result = self.runCommand(command)
        return result==0

    def remoteCopy(self, source, target, force=False):
        if self.client:
            status,_ = self.client.copy(source, target, force=force)
            if self.debug:
                print('DEBUG:', status.message, 'copy from', source, 'to', target, " force=", force)
            result = 0 if status.ok else 5
        else:
            command = self.remoteCp.format(source=source, target=target)
            result = self.runCommand(command)
        return result==0

    def cp(self, source, target, force=False):
        if self.isRemotePath(source) or self.isRemotePath(target):
            return self.remoteCopy(source, target, force=force)
        else:
            shutil.copyfile(source, target)
            return True

    def rm(self, path):
        if self.isRemotePath(path):
            self.remoteFileRm(path)
        else:
            os.remove(path)

    # attempts > 1 can be given to check existence of files multiple times
    #  this is to workaround a storage problem where files are randomly reported as not existent
    def exists(self, path, attempts=1):
        attemptsLeft = attempts
        found = False
        while attemptsLeft > 0:
            attemptsLeft -= 1
            found = self.remoteFileExists(path) if self.isRemotePath(path) else os.path.exists(path)
            if found:
                return True
            if attemptsLeft > 0:
                if attemptsLeft < attempts-1 and self.timeBetweenAttempts and self.timeBetweenAttempts<10:
                    self.timeBetweenAttempts *= 2
                for xrdenv in ['X509_VOMS_DIR','X509_VOMS_DIR','X509_USER_PROXY']:
                    if xrdenv in os.environ:
                        print('DEBUG:', xrdenv, '=', os.environ[xrdenv])
                    else:
                        print('DEBUG: not found:', xrdenv)
                print('INFO: file was not found:'+path+', trying %d more times...'%attemptsLeft)
                if self.timeBetweenAttempts:
                    reconnect = False
                    if self.client:
                        self.client = None
                        # don't reconnect and try to use fallback solution using xrdcp/xrdfs
                        reconnect = False
                        # try to reconnect
                        #reconnect = True
                    print('INFO: wait %r seconds before trying again'%self.timeBetweenAttempts)
                    sleep(self.timeBetweenAttempts)
                    if reconnect:
                        try:
                            try:
                                reload(XRootD.client)
                            except Exception as e:
                                print('INFO: > ', e)
                                try:
                                    import XRootD.client 
                                    print('INFO: > reloaded XrootD') 
                                except Exception as e2:
                                    print('INFO: >> ', e2)
                            self.client = XRootD.client.FileSystem(self.server)
                            print('INFO: > reconnected file system! ==>', self.client)
                        except Exception as e3:
                            print('INFO: >>> ', e3)

        return found

    def mkdir(self, path):
        status = False
        if not self.isRemotePath(path):
            try:
                os.mkdir(path)
                status = True
                if self.debug:
                    print('DEBUG: created the local directory:', path)
            except:
                pass
        else:
            command = self.remoteMkdir.format(server=self.getRemoteFileserver(), path=self.getLocalFileName(path))
            status = self.runCommand(command)
        return status

    # create folder and parent folders 
    def makedirs(self, path):
        status = True
        # use os functionality for local files
        if not self.isRemotePath(path):
            try:
                os.makedirs(path)
            except:
                status = False
        else:
            if not self.remoteDirectoryExists(path):
                localName = self.getLocalFileName(path.strip('/'))
                directories = localName.split('/')
                i = len(directories)
                while i>self.makedirsMinLevel:
                    parentDirectory = self.getRemoteFileName('/'.join(directories[:i-1]))
                    if self.remoteDirectoryExists(parentDirectory):
                        break
                    i -= 1
                for j in range(i, len(directories) + 1):
                    remoteName = self.getRemoteFileName('/'.join(directories[:j]))
                    self.mkdir(remoteName)
        return status

    def fileExists(self, fileName):
        return os.path.isfile(self.getLocalFileName(fileName))

    def directoryExists(self, fileName):
        return os.path.isdir(self.getLocalFileName(fileName))

    def getRemoteFileName(self, path):
        # only xrootd protocol supported at the moment
        return self.getXrootdFileName(path)

    def removeRedirector(self, fileName):
        if '://' in fileName:
            fileName = '/'.join(fileName.split('://')[1].split('/')[1:])
        return fileName

    def addRedirector(self, redirector=None, fileName=None):
        if redirector is None:
            redirector = self.xrootdRedirectors[0]
        if fileName is None:
            fileName = ""
        return redirector + '/' + fileName

    # ------------------------------------------------------------------------------
    # get file name WITH redirector
    # ------------------------------------------------------------------------------
    def getXrootdFileName(self, rawFileName):
        xrootdFileName = rawFileName.strip()
        if self.isStoragePath(xrootdFileName):
            if self.pnfsStoragePath:
                xrootdFileName = self.pnfsStoragePath + xrootdFileName.strip()
            if self.xrootdRedirectors:
                xrootdFileName = self.addRedirector(fileName=self.removeRedirector(xrootdFileName), redirector=self.xrootdRedirectors[0])
            return xrootdFileName.strip()
        elif self.isPnfs(xrootdFileName):
            # replace already existing redirectors with primary one
            if self.xrootdRedirectors:
                xrootdFileName = self.addRedirector(fileName=self.removeRedirector(xrootdFileName), redirector=self.xrootdRedirectors[0])
            return xrootdFileName.strip()
        else:
            return xrootdFileName.strip()

    # ------------------------------------------------------------------------------
    # get file name WITHOUT redirector, but with complete PNFS path in case of a 
    # /store/ path.
    # ------------------------------------------------------------------------------
    def getLocalFileName(self, rawFileName):
        if rawFileName:
            localFileName = self.removeRedirector(rawFileName.strip())
            if localFileName.startswith('/store/') and self.pnfsStoragePath:
                localFileName = self.pnfsStoragePath + localFileName.strip()
            return localFileName
        else:
            print ("\x1b[31mERROR: invalid file name\x1b[0m")
            raise Exception("InvalidFileName")

    def getXrootdRedirector(self):
        return self.xrootdRedirectors[0] if self.xrootdRedirectors and len(self.xrootdRedirectors)>0 else None

    # filesystem should be mounted read-only as normal nfs mount for e.g. glob operations, but since this very often fails on T3 worker nodes, use this 
    # as fallback when applicable! (can NOT fully replace glob, only list ALL files in a directory so use with care!)
    def lsRemote(self, path, allowCachedRead=True):
        if self.client:
            if allowCachedRead and self.useDirectoryListingCache and path in self.dirListingCache:
                listing = self.dirListingCache[path]
            else:
                # when python bindings are available, use them for directory listing
                status, listing = self.client.dirlist(path)
                self.dirListingCache[path] = listing
            if listing is None:
                return None 
            fileList = sorted([(x.name if x.name.startswith('/') else path + '/' + x.name) for x in listing])
        else:
            # fallback for the fallback...
            serverName = self.xrootdRedirectors[0].replace('root://','').split(':')[0].strip()
            lsRemoteCommand = ["xrdfs {serverName} ls {path}".format(serverName=serverName, path=path)]
            print("DEBUG: command=", lsRemoteCommand)
            counter = 0
            status = 0
            # since every call to xrdfs has a random chance to produce [Auth failed] errors on T3, try few times
            while status < 1 and counter < 5:
                counter += 1
                try:
                    out = subprocess.check_output(
                        lsRemoteCommand, stderr=subprocess.STDOUT, shell=True,
                        universal_newlines=True)
                    status = 1
                except subprocess.CalledProcessError as exc:
                    print("ERROR: xrdfs ls failed:", exc.returncode, exc.output)
            
            if counter > 4:
                raise Exception("RemoteFSError")

            fileList = []
            for fileName in sorted(out.strip("\n").split("\n")):
                if fileName.startswith('/'):
                    fileList.append(fileName)
                else:
                    fileList.append(path + '/' + fileName)
        return fileList

    # like glob but allows wildcards only in the last level
    def glob_XrootdPyBindings(self, path):
        path = self.removeRedirector(path.strip())
        basePath = '/'.join(path.split('/')[:-1])
        if '*' in basePath:
            raise Exception("NotImplemented")
        status, fileList = self.client.dirlist(basePath)
        if "ERROR" in str(status):
            print("\x1b[31mERROR: directory listing with python xrootd bindings failed:", status, "\x1b[0m")
            self.failureCounters['pyxrd'] += 1
            raise Exception("DirectoryListingFailed")

        if not path.startswith('/'):
            fileList = sorted(fnmatch.filter(fileList, path))
        fileList = [basePath + '/' + x if not x.startswith('/') else x for x in fileList]
        if path.startswith('/'):
            fileList = sorted(fnmatch.filter(fileList, path))

        return fileList

    def glob_Xrdfs(self, path):
        redirector = self.getRedirector(path)
        if redirector is None:
            redirector = self.xrootdRedirectors[0]
        path = self.removeRedirector(path.strip())
        basePath = '/'.join(path.split('/')[:-1])
        if '*' in basePath:
            raise Exception("NotImplemented")
        serverName = self.getHostnameFromRedirector(redirector)
        lsRemoteCommand = ["xrdfs {serverName} ls {path}".format(serverName=serverName, path=basePath)]
        print("DEBUG: dirlist fallback #1: command=", lsRemoteCommand)
        try:
            out = subprocess.check_output(
                lsRemoteCommand, stderr=subprocess.STDOUT, shell=True,
                universal_newlines=True)
            success = True

            fileList = out.strip("\n").split("\n")
            if not path.startswith('/'):
                fileList = sorted(fnmatch.filter(fileList, path))
            fileList = [basePath + '/' + x if not x.startswith('/') else x for x in fileList]
            if path.startswith('/'):
                fileList = sorted(fnmatch.filter(fileList, path))

        except Exception as exc:
            print("\x1b[31mERROR: dirlist fallback #1 (xrdfs ls) failed:", exc.returncode, exc.output, "\x1b[0m")
            self.failureCounters['xrdfs'] += 1
            raise exc

        return fileList

    def glob(self, path):
        if self.isRemotePath(path):
            if self.client:
                return self.glob_XrootdPyBindings(path)
            else:
                return self.glob_Xrdfs(path)
        else:
            return glob.glob(path)

    def glob_with_fallback(self, path):
        if self.isRemotePath(path):
            success = False

            # try xrootd bindings first
            if self.client and self.failureCounters['pyxrd'] < 3 and not self.disableXrootdForDirlist:
                try:
                    fileList = self.glob_XrootdPyBindings(path)
                    success  = True
                except Exception as e:
                    print("ERROR:", e)

            # try xrdfs .. ls
            if not success and self.failureCounters['xrdfs'] < 3 and not self.disableXrootdForDirlist:
                try:
                    fileList = self.glob_Xrdfs(path)
                    success = True
                except Exception as e:
                    print("ERROR:", e)

            # try ls on local fs
            if not success:
                path = self.removeRedirector(path)
                basePath = '/'.join(path.split('/')[:-1])
                if os.path.isdir(basePath):
                    try:
                        fileList = glob.glob(path)
                        success = True
                    except Exception as e:
                        print("\x1b[31mERROR: dirlist fallback #2 (local glob) failed:", e, "\x1b[0m")
                else:
                    print("\x1b[31mERROR: dirlist fallback #2 (local glob) failed: not a directory:", basePath, "\x1b[0m")

            if not success:
                print("\x1b[31mERROR: none of the methods to obtain the directory listing succeeded.\x1b[0m")
                raise Exception("DirectoryListingFailed")

            return fileList
        else:
            return glob.glob(path)

    def get_numbered_file_list(self, mask, start, end):
        if mask.count('*') > 1:
            print("INFO: filename mask contains more than one *!")
            raise Exception("NotImplemented")
        elif mask.count('*') < 1:
            print("WARNING: filename mask does not contain any * => use exists() instead!")
            # since there is only 1-file, check only once
            start = 1
            end = 1
        existingFiles = []
        for i in range(start, end+1):
            fileName = mask.replace('*','%d'%i)
            if self.exists(fileName):
                existingFiles.append(fileName)
        return existingFiles






